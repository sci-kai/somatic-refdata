"""Workflow profiles: turn a dataset into the samplesheet, params and run of one workflow.

A profile is a small YAML file (profiles/<workflow>.yaml, or any path):

  name: variantinterpretation
  pipeline: cio-abcd/variantinterpretation   # nextflow run target (repo or local path)
  revision: main                             # -r (left out for a local path)
  nf_profile: conda                          # -profile
  accepts:                                   # optional: the dataset must fit the workflow
    callers: [muTect2]
  requires: [reference]                      # optional: dataset resources that must exist
  paths: absolute                            # absolute or relative (to the samplesheet)
  samplesheet:
    file: samplesheet.csv
    rows: vcf                                # vcf: one row per VCF; pair: per pair x
                                             # technology; bam: per alignment
    missing: n/a                             # value of empty fields
    columns:                                 # column: template
      sample: "{pair}_{technology}_{caller}"
      vcf: "{vcf}"
  params:                                    # param: template (left out if empty)
    input: "{samplesheet}"
    fasta: "{reference}"
  files:                                     # extra files: relative path -> template
    known.csv: |
      ...

Templates are Python format strings. Row variables: pair, center, replicate, technology,
tech (lower case), library_type (wgs / capture), aligner, caller, vcf, tumor_bam,
tumor_bai, normal_bam, normal_bai, tumor_id, normal_id, target_bed, role, bam, bai.
Dataset variables (also in rows): dataset, dataset_dir, samplesheet, workdir, pairs,
first_pair, all_bams (true/false: every pair has tumor and normal BAMs), reference,
reference_fai, target_bed_wes, truth_snv, truth_indel,
truth_regions, gencode_gtf, giab_<name>. A variable without value is empty. Rendered
params 'true' / 'false' are written as booleans.
"""

import csv
import io
import json
import logging
import os
import string
from importlib import resources
from pathlib import Path

import yaml

from .dataset import ConfigError, merge

log = logging.getLogger("seqc2")

PROFILE_DEFAULTS = {
    "name": None,
    "description": "",
    "pipeline": None,
    "revision": None,
    "nf_profile": None,
    "nextflow_args": [],
    "accepts": {},
    "requires": [],
    "paths": "absolute",
    "samplesheet": {
        "file": "samplesheet.csv",
        "rows": "vcf",
        "missing": "n/a",
        "separator": ",",
        "columns": {},
    },
    "params": {},
    "files": {},
}
LIBRARY_TYPE = {"WGS": "wgs", "WES": "capture"}


def shipped_profiles():
    """Names of the profiles shipped with the toolkit (profiles/*.yaml in the package)."""
    folder = resources.files("seqc2_toolkit").joinpath("profiles")
    return sorted(
        p.name.removesuffix(".yaml") for p in folder.iterdir() if p.name.endswith(".yaml")
    )


def load_profile(name_or_path):
    """A profile by name (shipped) or path, merged with the defaults; checked for typos."""
    path = Path(name_or_path)
    if path.suffix in (".yaml", ".yml") and path.exists():
        text, where = path.read_text(), str(path)
    else:
        resource = resources.files("seqc2_toolkit").joinpath("profiles", f"{name_or_path}.yaml")
        if not resource.is_file():
            raise ConfigError(
                f"profile '{name_or_path}' not found; shipped profiles: "
                f"{', '.join(shipped_profiles())}, "
                "or give the path of a profile YAML"
            )
        text, where = resource.read_text(), f"profile {name_or_path}"
    profile = merge(PROFILE_DEFAULTS, yaml.safe_load(text) or {}, where)
    profile["name"] = profile["name"] or path.stem
    if profile["samplesheet"]["rows"] not in ("vcf", "pair", "bam"):
        raise ConfigError(f"{where}: samplesheet.rows must be vcf, pair or bam")
    if profile["paths"] not in ("absolute", "relative"):
        raise ConfigError(f"{where}: paths must be absolute or relative")
    return profile


class Values(dict):
    """Template variables; an unknown or empty variable renders as '' and is noted."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.empty = []

    def __missing__(self, key):
        self.empty.append(key)
        return ""


def render(template, values):
    """(text, empty variables used) of a template with the given variables."""
    values = Values({k: v for k, v in values.items() if v not in (None, "")})
    text = string.Formatter().vformat(str(template), (), values)
    return text, values.empty


def rows(facts, unit_kind, path):
    """Row variables: one dict per VCF, pair x technology, or alignment of the dataset."""
    tumor_id, normal_id = facts["vcf_samples"]["tumor"], facts["vcf_samples"]["normal"]
    result = []
    for unit in facts["units"]:
        base = {
            "pair": unit["pair"],
            "center": unit["center"],
            "replicate": unit["replicate"],
            "technology": unit["technology"],
            "tech": unit["technology"].lower(),
            "library_type": LIBRARY_TYPE[unit["technology"]],
            "target_bed": path(unit["target_bed"]),
            "tumor_id": tumor_id,
            "normal_id": normal_id,
            "tumor_bam": path(unit["bams"].get("T")),
            "tumor_bai": path(unit["bams"].get("T") and unit["bams"]["T"] + ".bai"),
            "normal_bam": path(unit["bams"].get("N")),
            "normal_bai": path(unit["bams"].get("N") and unit["bams"]["N"] + ".bai"),
        }
        if unit_kind == "vcf":
            for vcf in unit["vcfs"]:
                result.append(
                    {
                        **base,
                        "aligner": vcf["aligner"],
                        "caller": vcf["caller"],
                        "vcf": path(vcf["path"]),
                    }
                )
        elif unit_kind == "pair":
            first = unit["vcfs"][0] if unit["vcfs"] else {}
            result.append(
                {
                    **base,
                    "aligner": first.get("aligner"),
                    "caller": first.get("caller"),
                    "vcf": path(first.get("path")),
                }
            )
        else:
            for role, bam in unit["bams"].items():
                result.append(
                    {
                        **base,
                        "role": role,
                        "aligner": "bwa",
                        "bam": path(bam),
                        "bai": path(bam + ".bai"),
                    }
                )
    return result


def check(profile, facts):
    """Stop if the dataset doesn't fit the profile (accepts, requires)."""
    problems = []
    for field, allowed in profile["accepts"].items():
        if field not in ("technologies", "callers", "aligners"):
            raise ConfigError(
                f"profile {profile['name']}: accepts.{field} unknown "
                "(technologies, callers, aligners)"
            )
        key = field.removesuffix("s").replace("technologie", "technology")
        values = (
            {u["technology"] for u in facts["units"]}
            if key == "technology"
            else {v[key] for u in facts["units"] for v in u["vcfs"]}
        )
        bad = sorted(values - set(allowed))
        if bad:
            problems.append(
                f"{field}: {', '.join(bad)} not supported by {profile['name']} "
                f"(supported: {', '.join(allowed)})"
            )
    available = set(facts["resources"])
    if any(u["bams"] for u in facts["units"]):
        available.add("bam")
    if any(u["vcfs"] for u in facts["units"]):
        available.add("vcf")
    for need in profile["requires"]:
        if need not in available:
            problems.append(
                f"requires {need}: not in the dataset "
                "(set it in the data section of the dataset config)"
            )
    if profile["samplesheet"]["rows"] == "bam" and not any(u["bams"] for u in facts["units"]):
        problems.append("samplesheet rows are alignments, but the dataset has no BAMs (data.bam)")
    if problems:
        raise ConfigError(
            f"the dataset doesn't fit profile {profile['name']}:\n  " + "\n  ".join(problems)
        )


def write(profile, dataset_dir, workdir=None):
    """Write samplesheet, params.yaml and extra files of a profile; return their paths.

    workdir: folder of the workflow files (default: <dataset>/<profile name>/).
    """
    dataset_dir = Path(dataset_dir).resolve()
    facts = json.loads((dataset_dir / "dataset.json").read_text())
    check(profile, facts)
    workdir = Path(workdir).resolve() if workdir else dataset_dir / profile["name"]
    workdir.mkdir(parents=True, exist_ok=True)
    sheet = profile["samplesheet"]
    sheet_path = workdir / sheet["file"]

    def path(rel):
        """A dataset path as written into the workflow files (absolute or relative)."""
        if not rel:
            return None
        full = dataset_dir / rel
        if profile["paths"] == "relative":
            return os.path.relpath(full, workdir)
        return str(full)

    pairs = [u["pair"] for u in facts["units"]]
    dataset_values = {
        "dataset": facts["name"],
        "dataset_dir": str(dataset_dir),
        "workdir": str(workdir),
        "samplesheet": (
            path(str(sheet_path.relative_to(dataset_dir)))
            if sheet_path.is_relative_to(dataset_dir)
            else str(sheet_path)
        ),
        "pairs": ",".join(dict.fromkeys(pairs)),
        "first_pair": pairs[0] if pairs else None,
        # 'true' if every unit has tumor and normal BAMs (e.g. for params that need them)
        "all_bams": str(
            bool(facts["units"]) and all(len(u["bams"]) == 2 for u in facts["units"])
        ).lower(),
        **{k: path(v) for k, v in facts["resources"].items()},
    }
    # samplesheet: one row per unit of the chosen kind
    out = io.StringIO()
    writer = csv.writer(out, delimiter=sheet["separator"], lineterminator="\n")
    writer.writerow(list(sheet["columns"]))
    for row in rows(facts, sheet["rows"], path):
        values = {**dataset_values, **row}
        writer.writerow(
            [render(t, values)[0] or sheet["missing"] for t in sheet["columns"].values()]
        )
    sheet_path.write_text(out.getvalue())
    written = [sheet_path]
    # params: a param with an empty variable is left out (e.g. no gencode_gtf in the dataset)
    params = {}
    for key, template in profile["params"].items():
        if isinstance(template, (bool, int, float)):
            params[key] = template
            continue
        text, empty = render(template, dataset_values)
        if empty:
            log.info("param %s left out: %s not in the dataset", key, ", ".join(empty))
            continue
        # 'true' / 'false' become booleans (a string 'false' would be true in nextflow)
        params[key] = {"true": True, "false": False}.get(text, text)
    params_path = workdir / "params.yaml"
    params_path.write_text(
        f"# Parameters of {profile['name']} for the SEQC2 dataset {facts['name']} (seqc2-toolkit)\n"
        + yaml.safe_dump(params, sort_keys=False)
    )
    written.append(params_path)
    # extra files (e.g. a known-variant sheet); an empty variable is an error here
    for rel, template in profile["files"].items():
        text, empty = render(template, dataset_values)
        if empty:
            raise ConfigError(
                f"profile {profile['name']}: file {rel} needs {', '.join(empty)}, "
                "which the dataset lacks"
            )
        target = workdir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
        written.append(target)
    return written


def resolve_params(workdir, rundir):
    """params.yaml of the workflow folder with absolute paths, written into the run folder.

    Relative values that name an existing file or folder (relative to the workflow folder,
    paths: relative) are made absolute: nextflow resolves them against the launch folder.
    """
    workdir, rundir = Path(workdir), Path(rundir)
    params = yaml.safe_load((workdir / "params.yaml").read_text()) or {}
    for key, value in params.items():
        if isinstance(value, str) and not Path(value).is_absolute() and (workdir / value).exists():
            params[key] = str((workdir / value).resolve())
    rundir.mkdir(parents=True, exist_ok=True)
    target = rundir / "params.yaml"
    target.write_text(yaml.safe_dump(params, sort_keys=False))
    return target


def run_command(
    profile,
    params_file,
    outdir,
    pipeline=None,
    revision=None,
    nf_profile=None,
    configs=(),
    extra=(),
):
    """The nextflow command of a profile (list), with optional overrides."""
    pipeline = pipeline or profile["pipeline"]
    if not pipeline:
        raise ConfigError(f"profile {profile['name']}: no pipeline; give --pipeline")
    cmd = ["nextflow", "run", pipeline]
    revision = revision or profile["revision"]
    # a local pipeline folder has no revision
    if revision and not Path(pipeline).exists():
        cmd += ["-r", revision]
    nf_profile = nf_profile or profile["nf_profile"]
    if nf_profile:
        cmd += ["-profile", nf_profile]
    cmd += ["-params-file", str(params_file), "--outdir", str(outdir)]
    for config in configs:
        cmd += ["-c", str(Path(config).resolve())]
    return cmd + list(profile["nextflow_args"]) + list(extra)
