"""Dataset config (YAML) and the selection of reference data files it describes.

source: seqc2 (default; SEQC2 HCC1395/HCC1395BL, below) or hg008 (GIAB HG008, see hg008.py).
A SEQC2 dataset config has three optional sections; everything left out takes the default:

  access: local           local: download and build the dataset folder (seqc2 fetch);
                          remote: no download, the workflow reads the files from the
                          SEQC2 server (URLs; whole files only, see build.remote_facts)

  select:                 which tumor/normal pairs, technologies, aligners, callers
    technologies: [WGS]   WGS and/or WES
    centers: [IL]         sequencing centers, or 'all'
    replicates: [1]       replicate numbers, or 'all'
    pairs: null           explicit pairs (e.g. [IL_1, FD_2]); replaces centers/replicates
    aligners: [bwa]       bwa, bowtie, novo
    callers: [muTect2]    muTect2, strelka, somaticSniper
  data:                   which kinds of files the dataset contains
    vcf: true             the somatic VCFs of the selection
    bam: false            alignments: false, true (every selected pair) or a list of pairs
    reference: true       GRCh38.d1.vd1 FASTA (+ .fai)
    target_bed: true      exome target BED (used for WES)
    truth: v1.2.1         SEQC2 somatic truth set (release), or false
    gencode: false        GENCODE basic GTF (release number, e.g. 46), or false
    giab: []              GIAB stratifications (homopolymer, tandem_repeat, ...)
  subset:                 make a small dataset (default: everything genome-wide)
    contigs: null         keep only these contigs in every file, e.g. [chr22]
    bam_windows: null     BAMs: read only these windows (e.g. chr22:23000000-24000000)
    bam_targets: true     BAMs: within the windows, only the WES targets (+ padding)
    bam_padding: 100      bp around the targets
    bam_subsample: null   fraction of read pairs to keep (fixed seed 42)
  vcf:                    optional changes of the VCFs (local datasets; see vcfprep.py)
    sample_order: null    sample columns in this order, e.g. [NORMAL, TUMOR]
    format_dp: false      add FORMAT/DP = sum(FORMAT/AD) to records without FORMAT/DP

Default (empty config): WGS of pair IL_1 (Illumina center, replicate 1), bwa + MuTect2
VCF, reference and truth set; genome-wide.
"""

import logging
from pathlib import Path

import yaml

log = logging.getLogger("seqc2")

SOURCE_NAMES = ("seqc2", "hg008")
ACCESS = ("local", "remote")

DEFAULTS = {
    "source": "seqc2",
    "name": None,
    "access": "local",
    "outdir": None,
    "cache": None,
    "select": {
        "technologies": ["WGS"],
        "centers": ["IL"],
        "replicates": [1],
        "pairs": None,
        "aligners": ["bwa"],
        "callers": ["muTect2"],
    },
    "data": {
        "vcf": True,
        "bam": False,
        "reference": True,
        "target_bed": True,
        "truth": "v1.2.1",
        "gencode": False,
        "giab": [],
    },
    "subset": {
        "contigs": None,
        "bam_windows": None,
        "bam_targets": True,
        "bam_padding": 100,
        "bam_subsample": None,
    },
    "vcf": {
        "sample_order": None,
        "format_dp": False,
    },
}


class ConfigError(Exception):
    """An invalid dataset config or profile (the message says what to change)."""


def as_list(value):
    """A single value or a list as a list (None stays None)."""
    if value is None or isinstance(value, list):
        return value
    return [value]


def load(path, overrides=None):
    """The dataset config of a YAML file, merged with the defaults; checked for unknown keys.

    Relative 'outdir' and 'cache' are resolved against the folder of the config file.
    """
    path = Path(path)
    user = yaml.safe_load(path.read_text()) or {}
    if not isinstance(user, dict):
        raise ConfigError(f"{path}: expected a mapping with the sections select, data, subset")
    config = merge(defaults_of(user.get("source", "seqc2"), str(path)), user, str(path))
    # paths in the config file are relative to its folder
    for key in ("outdir", "cache"):
        if config[key] is not None:
            config[key] = str((path.parent / config[key]).resolve())
    # command-line values replace them; they are relative to the working directory
    for key, value in (overrides or {}).items():
        if value is not None:
            config[key] = str(Path(value).resolve()) if key in ("outdir", "cache") else value
    if config["access"] not in ACCESS:
        raise ConfigError(f"{path}: access must be one of {', '.join(ACCESS)}")
    # remote datasets (URLs instead of downloaded files) exist for SEQC2 only
    if config["access"] == "remote" and config["source"] != "seqc2":
        raise ConfigError(
            f"{path}: access 'remote' is only available for source seqc2 "
            f"(source is {config['source']}); use access: local"
        )
    # the dataset name defaults to the config file name
    config["name"] = config["name"] or path.stem
    config["config_file"] = str(path.resolve())
    return config


def defaults_of(source, where="config"):
    """The config defaults of a source: seqc2 (HCC1395) or hg008 (GIAB HG008)."""
    if source == "seqc2":
        return DEFAULTS
    if source == "hg008":
        from . import hg008

        return hg008.DEFAULTS
    raise ConfigError(f"{where}: source {source} unknown; choose from {', '.join(SOURCE_NAMES)}")


def override_select(config, technologies=None, pairs=None, centers=None, replicates=None):
    """Replace select values of a config with command-line values (comma-separated or 'all').

    Lets the user choose the samples of a dataset without editing its config, e.g.
    --pairs IL_1,FD_1 or --centers all --replicates all. Giving centers or replicates
    clears the pairs of the config (pairs would replace them). The changes are logged.
    """
    values = {
        "technologies": technologies,
        "pairs": pairs,
        "centers": centers,
        "replicates": replicates,
    }
    sel = config["select"]
    if (centers or replicates) and not pairs and sel["pairs"]:
        log.info("select.pairs %s of the config replaced by centers/replicates", sel["pairs"])
        sel["pairs"] = None
    for key, text in values.items():
        if text is None:
            continue
        # 'all' stays a word; everything else becomes a list of the comma-separated items
        new = "all" if text.strip() == "all" else [v.strip() for v in text.split(",") if v.strip()]
        if not new:
            raise ConfigError(f"--{key}: no value given (comma-separated list or 'all')")
        log.info("select.%s: %s (command line; config: %s)", key, new, sel[key])
        sel[key] = new
    return config


def merge(defaults, user, where):
    """defaults updated with user values; unknown keys are an error (typos must not pass)."""
    unknown = sorted(set(user) - set(defaults))
    if unknown:
        raise ConfigError(
            f"{where}: unknown key(s) {', '.join(unknown)}; allowed: {', '.join(defaults)}"
        )
    merged = {}
    for key, default in defaults.items():
        if key not in user:
            merged[key] = default
        elif isinstance(default, dict):
            if not isinstance(user[key], dict):
                raise ConfigError(f"{where}: '{key}' must be a mapping")
            # an empty default is a free mapping (e.g. params, columns): any keys allowed
            merged[key] = merge(default, user[key], f"{where}: {key}") if default else user[key]
        else:
            merged[key] = user[key]
    return merged


def pair_order(pair):
    """Sort key of a pair name (center, replicate number)."""
    center, replicate = pair.split("_")
    return center, int(replicate)


def choose(name, requested, available):
    """Requested values ('all' = every available one), checked against the available ones."""
    if requested == "all" or requested == ["all"]:
        return sorted(available)
    requested = [str(v) if name != "replicates" else int(v) for v in as_list(requested)]
    unknown = [v for v in requested if v not in available]
    if unknown:
        raise ConfigError(
            f"select.{name}: {', '.join(map(str, unknown))} not available; "
            f"choose from: {', '.join(map(str, sorted(available)))}"
        )
    return requested


def select(config, records):
    """The selected units: one per pair x technology, with its VCFs and BAMs.

    Every unit is a dict: pair, center, replicate, technology, vcfs (list of catalog
    records, one per aligner x caller), bams ({'T': record, 'N': record} or {}),
    bais (same for the indices). Missing files of a requested combination are an error.
    HG008 configs (source: hg008) are selected by hg008.select.
    """
    if config.get("source", "seqc2") == "hg008":
        from . import hg008

        return hg008.select(config, records)
    sel, data = config["select"], config["data"]
    vcfs = [r for r in records if r["kind"] == "vcf"]
    technologies = choose("technologies", sel["technologies"], {r["technology"] for r in records})
    aligners = choose("aligners", sel["aligners"], {r["aligner"] for r in vcfs})
    callers = choose("callers", sel["callers"], {r["caller"] for r in vcfs})
    # pairs: explicit list (kept in the given order), or centers x replicates (sorted)
    if sel["pairs"]:
        pairs = list(dict.fromkeys(choose("pairs", sel["pairs"], {r["pair"] for r in records})))
    else:
        centers = choose("centers", sel["centers"], {r["center"] for r in records})
        replicates = choose("replicates", sel["replicates"], {r["replicate"] for r in records})
        pairs = sorted(
            {r["pair"] for r in records if r["center"] in centers and r["replicate"] in replicates},
            key=pair_order,
        )
    if not pairs:
        raise ConfigError("select: no pair matches the centers and replicates")
    # BAMs: none, all selected pairs, or a subset of them
    bam_pairs = set()
    if data["bam"] is True:
        bam_pairs = set(pairs)
    elif data["bam"]:
        bam_pairs = set(as_list(data["bam"]))
        unknown = sorted(bam_pairs - set(pairs))
        if unknown:
            raise ConfigError(
                f"data.bam: {', '.join(unknown)} not among the selected pairs ({', '.join(pairs)})"
            )

    units, missing = [], []
    for pair in pairs:
        for technology in technologies:
            have = [r for r in records if r["pair"] == pair and r["technology"] == technology]
            if not have:
                # e.g. NS (WGS only) with WES: the pair doesn't exist for this technology
                missing.append(f"{pair} {technology} (pair not sequenced with this technology)")
                continue
            unit_vcfs = []
            for aligner in aligners:
                for caller in callers:
                    found = [
                        r
                        for r in have
                        if r["kind"] == "vcf" and r["aligner"] == aligner and r["caller"] == caller
                    ]
                    if data["vcf"] and not found:
                        missing.append(f"{pair} {technology} VCF {aligner}/{caller}")
                    unit_vcfs += found
            bams, bais = {}, {}
            if pair in bam_pairs:
                for role in ("T", "N"):
                    bam = [r for r in have if r["kind"] == "bam" and r["role"] == role]
                    bai = [r for r in have if r["kind"] == "bai" and r["role"] == role]
                    if not bam or not bai:
                        missing.append(f"{pair} {technology} BAM {role}")
                        continue
                    bams[role], bais[role] = bam[0], bai[0]
            units.append(
                {
                    "pair": pair,
                    "center": have[0]["center"],
                    "replicate": have[0]["replicate"],
                    "technology": technology,
                    "vcfs": unit_vcfs if data["vcf"] else [],
                    "bams": bams,
                    "bais": bais,
                }
            )
    if missing:
        raise ConfigError(
            "the selection asks for files that SEQC2 doesn't provide:\n  "
            + "\n  ".join(missing)
            + "\nnarrow the selection (e.g. other centers, technologies or callers)"
        )
    return units


def describe(units):
    """One line per unit: what the selection contains."""
    lines = []
    for unit in units:
        vcfs = ", ".join(f"{v['aligner']}/{v['caller']}" for v in unit["vcfs"]) or "no VCF"
        bams = "tumor + normal BAM" if unit["bams"] else "no BAM"
        lines.append(f"{unit['pair']} {unit['technology']}: {vcfs}; {bams}")
    return "\n".join(lines)
