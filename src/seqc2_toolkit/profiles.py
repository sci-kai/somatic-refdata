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
  params_remote:                             # optional: replaces params of a remote dataset
    callable_regions: false                  # (access: remote); null removes a param
  files_remote: {}                           # the same for extra files (null: not written)
  files:                                     # extra files: relative path -> template
    known.csv: |
      ...

Templates are Python format strings. Row variables: pair, center, replicate, technology,
tech (lower case), library_type (wgs / capture), aligner, caller, vcf, tumor_bam,
tumor_bai, normal_bam, normal_bai, tumor_id, normal_id, target_bed, role, bam, bai;
HG008 only: series, platform, build, tumor_sample, normal_sample, tumor_coverage,
normal_coverage (bai: .crai of a CRAM).
Dataset variables (also in rows): dataset, source (seqc2 / hg008), source_name
(SEQC2 / HG008), dataset_dir, samplesheet, workdir, pairs, first_pair, all_bams
(true/false: every pair has tumor and normal BAMs), reference, reference_fai,
target_bed_wes, truth_snv, truth_indel, truth_regions, gencode_gtf, giab_<name>, access
(local / remote), contigs (comma-separated contigs the records are, or for remote datasets
should be, restricted to); HG008 only: reference_gzi, truth_sv, truth_sv_regions,
truth_sv_regions_clonal. A variable without value is empty. Rendered params 'true' /
'false' are written as booleans. accepts may also restrict sources (seqc2, hg008).

Remote datasets (access: remote, SEQC2 only) hold URLs instead of paths; they are written
as given. Next to params.yaml, params.config holds the same params as a Nextflow config
(for -c or includeConfig); paths inside the workflow folder start with config_base there.
"""

import csv
import io
import json
import logging
import os
import re
import string
from importlib import resources
from pathlib import Path

import yaml

from .dataset import ConfigError, merge

log = logging.getLogger("seqc2")

PROFILE_DEFAULTS = {
    "name": None,
    "description": "",
    # the one workflow a profile is written for (shown by 'seqc2 profiles'): homepage,
    # folder of its dataset configs, notes on revisions or required pipeline params
    "workflow": {"url": None, "configs": None, "notes": None},
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
    "params_remote": {},
    "files": {},
    "files_remote": {},
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


def index_of(alignment):
    """Index path of a dataset alignment (.bam.bai or .cram.crai), None without one."""
    if not alignment:
        return None
    return alignment + (".crai" if alignment.endswith(".cram") else ".bai")


def vcf_samples(vcf, tumor_id, normal_id):
    """tumor_id / normal_id of a VCF: its own sample columns (HG008) or the dataset's."""
    if "tumor_id" in vcf:
        return {"tumor_id": vcf["tumor_id"], "normal_id": vcf.get("normal_id")}
    return {"tumor_id": tumor_id, "normal_id": normal_id}


def rows(facts, unit_kind, path):
    """Row variables: one dict per VCF, pair x technology, or alignment of the dataset."""
    tumor_id, normal_id = facts["vcf_samples"]["tumor"], facts["vcf_samples"]["normal"]
    result = []
    for unit in facts["units"]:
        samples, coverage = unit.get("samples") or {}, unit.get("coverage") or {}
        base = {
            "pair": unit["pair"],
            "center": unit["center"],
            "replicate": unit["replicate"],
            "technology": unit["technology"],
            "tech": unit["technology"].lower(),
            "library_type": unit.get("library_type") or LIBRARY_TYPE[unit["technology"]],
            "series": unit.get("series"),
            "platform": unit.get("platform"),
            "build": unit.get("build"),
            "tumor_sample": samples.get("T"),
            "normal_sample": samples.get("N"),
            "tumor_coverage": coverage.get("T"),
            "normal_coverage": coverage.get("N"),
            "target_bed": path(unit["target_bed"]),
            "tumor_id": tumor_id,
            "normal_id": normal_id,
            "tumor_bam": path(unit["bams"].get("T")),
            "tumor_bai": path(index_of(unit["bams"].get("T"))),
            "normal_bam": path(unit["bams"].get("N")),
            "normal_bai": path(index_of(unit["bams"].get("N"))),
        }
        if unit_kind == "vcf":
            for vcf in unit["vcfs"]:
                result.append(
                    {
                        **base,
                        **vcf_samples(vcf, tumor_id, normal_id),
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
                    **vcf_samples(first, tumor_id, normal_id),
                    "aligner": first.get("aligner") or unit.get("aligner"),
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
                        "aligner": unit.get("aligner", "bwa"),
                        "bam": path(bam),
                        "bai": path(index_of(bam)),
                    }
                )
    return result


def check(profile, facts):
    """Stop if the dataset doesn't fit the profile (accepts, requires)."""
    problems = []
    for field, allowed in profile["accepts"].items():
        if field not in ("sources", "technologies", "callers", "aligners"):
            raise ConfigError(
                f"profile {profile['name']}: accepts.{field} unknown "
                "(sources, technologies, callers, aligners)"
            )
        key = field.removesuffix("s").replace("technologie", "technology")
        if key == "source":
            values = {facts.get("source", "seqc2")}
        elif key == "technology":
            values = {u["technology"] for u in facts["units"]}
        else:
            values = {v[key] for u in facts["units"] for v in u["vcfs"]}
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


def is_url(value):
    """True for a URL (scheme://...), e.g. a file of a remote dataset."""
    return isinstance(value, str) and re.match(r"^[A-Za-z][A-Za-z0-9+.-]*://", value) is not None


def write(profile, dataset_dir, workdir=None, config_base=None):
    """Write samplesheet, params.yaml, params.config and extra files; return their paths.

    workdir: folder of the workflow files (default: <dataset>/<profile name>/).
    config_base: Nextflow expression for the workflow folder in params.config, e.g.
    '${projectDir}/assets/seqc2/test' when the files are committed into a pipeline
    (default: the absolute path of the workflow folder).
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
        # files of a remote dataset stay URLs
        if is_url(rel):
            return rel
        full = dataset_dir / rel
        if profile["paths"] == "relative":
            return os.path.relpath(full, workdir)
        return str(full)

    pairs = [u["pair"] for u in facts["units"]]
    dataset_values = {
        "dataset": facts["name"],
        # seqc2 / hg008, and the name used in sample names (SEQC2 / HG008)
        "source": facts.get("source", "seqc2"),
        "source_name": facts.get("source_name", "SEQC2"),
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
        "access": facts.get("access", "local"),
        "contigs": ",".join(facts.get("contigs") or []),
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
    templates = dict(profile["params"])
    # a remote dataset may need other params (e.g. no BAM-based features); null removes one
    if dataset_values["access"] == "remote":
        templates.update(profile["params_remote"])
    for key, template in templates.items():
        if template is None:
            log.info("param %s left out: null in params_remote", key)
            continue
        if isinstance(template, bool | int | float):
            params[key] = template
            continue
        text, empty = render(template, dataset_values)
        if empty:
            log.info("param %s left out: %s not in the dataset", key, ", ".join(empty))
            continue
        # 'true' / 'false' become booleans (a string 'false' would be true in nextflow)
        params[key] = {"true": True, "false": False}.get(text, text)
    params_path = workdir / "params.yaml"
    header = (
        f"# Parameters of {profile['name']} for the {dataset_values['source_name']} dataset "
        f"{facts['name']} (seqc2-toolkit)\n"
    )
    yaml_params = params
    if config_base:
        # files meant to be committed (e.g. into a pipeline): no path of this machine;
        # files of the workflow folder are written relative to it
        yaml_params = {k: local_relpath(v, workdir) or v for k, v in params.items()}
        header += "# paths relative to this folder\n"
    params_path.write_text(header + yaml.safe_dump(yaml_params, sort_keys=False))
    written.append(params_path)
    # extra files (e.g. a known-variant sheet); an empty variable is an error here
    files = dict(profile["files"])
    if dataset_values["access"] == "remote":
        files.update(profile["files_remote"])
    for rel, template in files.items():
        if template is None:
            log.info("file %s not written: null in files_remote", rel)
            continue
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
    # params as Nextflow config, after the extra files (a param may name one of them)
    config_path = workdir / "params.config"
    config_path.write_text(
        f"// Parameters of {profile['name']} for the SEQC2 dataset {facts['name']} "
        "(seqc2-toolkit; do not edit, regenerate with seqc2 samplesheet)\n"
        + nextflow_params(params, workdir, config_base)
    )
    written.append(config_path)
    return written


def groovy_string(text):
    """A single-quoted Groovy string literal (no interpolation)."""
    return "'" + text.replace("\\", "\\\\").replace("'", "\\'") + "'"


def local_relpath(value, workdir):
    """Path of a local file relative to the workflow folder, or None (URL, other value).

    A file is local if value is an existing path relative to the workflow folder
    (paths: relative) or an absolute path below it.
    """
    workdir = Path(workdir).resolve()
    if not isinstance(value, str) or is_url(value):
        return None
    if Path(value).is_absolute():
        full = Path(value).resolve()
        return os.path.relpath(full, workdir) if full.is_relative_to(workdir) else None
    return value if (workdir / value).exists() else None


def nextflow_params(params, workdir, config_base=None):
    """params as a Nextflow config block; files in the workflow folder start with config_base.

    A string value is a local file if it is a path relative to the workflow folder that
    exists (paths: relative) or an absolute path below the workflow folder. With
    config_base (a Nextflow expression such as ${projectDir}/assets/x) such a value becomes
    "<config_base>/<relative path>", else its absolute path. URLs and other values are
    written as they are.
    """
    workdir = Path(workdir).resolve()
    lines = ["params {"]
    width = max((len(k) for k in params), default=0)
    for key, value in params.items():
        if isinstance(value, bool):
            lines.append(f"    {key:<{width}} = {'true' if value else 'false'}")
            continue
        if isinstance(value, int | float):
            lines.append(f"    {key:<{width}} = {value}")
            continue
        value = str(value)
        rel = local_relpath(value, workdir)
        if rel is None:
            text = groovy_string(value)
        elif config_base:
            # interpolated Groovy string: only the base expression is evaluated
            escaped = rel.replace("\\", "\\\\").replace('"', '\\"').replace("$", "\\$")
            text = f'"{config_base}/{escaped}"'
        else:
            text = groovy_string(str((workdir / rel).resolve()))
        lines.append(f"    {key:<{width}} = {text}")
    lines.append("}")
    return "\n".join(lines) + "\n"


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
