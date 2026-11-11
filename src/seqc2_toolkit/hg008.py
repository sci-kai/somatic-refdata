"""GIAB cancer reference HG008: catalog, dataset config and selection (sources_hg008.yaml).

HG008-T (pancreatic cancer cell line) with two normals of the same donor, HG008-N-D
(duodenal tissue) and HG008-N-P (pancreatic tissue), sequenced by many centers with
short- and long-read platforms. A dataset unit is one tumor/normal pair (T-ND or T-NP)
sequenced in one series (platform + center + run, e.g. illumina-bcm, hifi-bcm, ont-ul-ucsc).

Catalog records (from the md5 list of every series and caller folder):

  kind      aln (BAM or CRAM) or vcf
  series    series id (sources_hg008.yaml)
  sample    T, N-D, N-P (alignments)
  build     GRCh38 (GRCh38-GIABv3) or GRCh37 (hs37d5)
  coverage  mean coverage from the file name (alignments)
  format    bam or cram; index: {url, md5} of the .bai / .crai
  caller    muTect2, strelka2-snv, strelka2-indel, deepsomatic (VCFs)
  url, md5

Dataset config (source: hg008); everything left out takes the default:

  select:
    series: [illumina-nygc]  sequencing series, or all
    normals: [N-D]           N-D and/or N-P: pairs T-ND, T-NP
    build: GRCh38            GRCh38 (GIABv3) or GRCh37 (hs37d5)
    callers: [muTect2]       caller VCFs, for the selected series that have them
  data:
    vcf: true                the caller VCFs of the selection
    bam: false               false, true (every pair) or a list of pairs (T-ND, T-NP)
    reference: true          GIAB mapping reference of the build (+ .fai, .gzi)
    truth: v0.4              small-variant somatic draft benchmark, or false
    truth_tier: all          benchmark regions: all, nogermlineinterference,
                             nogermlinewithin50bp
    truth_sv: false          SV somatic draft benchmark (v0.6, GRCh38 only), or false
    gencode: false           GENCODE basic GTF release (GRCh38 only), or false
    giab: []                 GIAB stratifications of the build
  subset: as for SEQC2 (contigs, bam_windows, bam_subsample; no capture targets)
"""

import copy
import logging
import re

from . import catalog
from .dataset import ConfigError, as_list, choose

log = logging.getLogger("seqc2")

PAIR_NAMES = {"N-D": "T-ND", "N-P": "T-NP"}

DEFAULTS = {
    "source": "hg008",
    "name": None,
    "outdir": None,
    "cache": None,
    "select": {
        "series": ["illumina-nygc"],
        "normals": ["N-D"],
        "build": "GRCh38",
        "callers": ["muTect2"],
    },
    "data": {
        "vcf": True,
        "bam": False,
        "reference": True,
        "truth": "v0.4",
        "truth_tier": "all",
        "truth_sv": False,
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
}


def load_sources():
    """The HG008 server description (sources_hg008.yaml)."""
    return catalog.load_sources("hg008")


def strip_path(path):
    """A path of an md5 list relative to its folder ('./x', '/x' and 'x' are the same)."""
    return path.removeprefix("./").lstrip("/")


def parse_md5(lines):
    """{relative path: md5} of an md5 list ('<md5>  <path>' or tab-separated)."""
    result = {}
    for line in lines:
        parts = line.split()
        if len(parts) == 2 and re.fullmatch(r"[0-9a-f]{32}", parts[0]):
            result[strip_path(parts[1])] = parts[0]
    return result


def parse_alignments(series_id, spec, md5s, sources, folder_url):
    """Alignment records of one series from its md5 list (files without index skipped)."""
    pattern = re.compile(sources["alignment_pattern"])
    coverage = re.compile(sources["coverage_pattern"])
    records = []
    for name, md5 in sorted(md5s.items()):
        match = pattern.match(name)
        if not match or "/" in name:
            continue
        if spec.get("match") and spec["match"] not in name:
            continue
        if spec.get("exclude") and spec["exclude"] in name:
            continue
        rest = match["rest"]
        builds = [b for b, token in sources["build_tokens"].items() if re.search(token, rest)]
        if len(builds) != 1:
            continue  # unaligned reads, CHM13 or other references
        index_name = f"{name}.{'bai' if match['format'] == 'bam' else 'crai'}"
        if index_name not in md5s:
            log.warning("%s: %s has no index in the md5 list; skipped", series_id, name)
            continue
        cov = coverage.search(rest)
        records.append(
            {
                "kind": "aln",
                "series": series_id,
                "sample": match["sample"],
                "build": builds[0],
                "coverage": int(cov["coverage"]) if cov else None,
                "format": match["format"],
                "name": name,
                "url": f"{folder_url}/{name}",
                "md5": md5,
                "index": {"url": f"{folder_url}/{index_name}", "md5": md5s[index_name]},
            }
        )
    return records


def build_catalog(cache, refresh=False):
    """All HG008 catalog records: alignments of every series, caller VCFs."""
    sources = load_sources()
    base = sources["base_url"]
    data = f"{base}/{sources['data_folder']}"
    lists = {}

    def md5s(folder_url, name):
        """Parsed md5 list of a folder (fetched once per run, cached on disk)."""
        url = f"{folder_url}/{name}"
        if url not in lists:
            lists[url] = parse_md5(catalog.md5_list(cache, url, refresh))
        return lists[url]

    records = []
    for series_id, spec in sources["series"].items():
        folder_url = f"{data}/{spec['folder']}"
        found = parse_alignments(
            series_id, spec, md5s(folder_url, spec["md5"]), sources, folder_url
        )
        # one alignment per sample and build: otherwise the series must be split
        seen = {}
        for r in found:
            key = (r["sample"], r["build"])
            if key in seen:
                raise RuntimeError(
                    f"sources_hg008.yaml: series {series_id} has two alignments of "
                    f"{r['sample']} {r['build']}: {seen[key]}, {r['name']} (add match/exclude)"
                )
            seen[key] = r["name"]
        records += found
    for caller, spec in sources["callers"].items():
        for build, build_folder in spec["builds"].items():
            folder_url = f"{data}/{spec['folder'].format(build_folder=build_folder)}"
            rel = spec["file"].format(build_folder=build_folder)
            md5 = md5s(folder_url, spec["md5"]).get(rel)
            if md5 is None:
                log.warning("caller %s %s: %s not in the md5 list; skipped", caller, build, rel)
                continue
            records.append(
                {
                    "kind": "vcf",
                    "series": spec["series"],
                    "caller": caller,
                    "normal": spec["normal"],
                    "build": build,
                    "samples": spec["samples"],
                    "name": rel.rsplit("/", 1)[-1],
                    "url": f"{folder_url}/{rel}",
                    "md5": md5,
                }
            )
    return records


def summary(records):
    """Text overview of the HG008 catalog: alignments per series, caller VCFs."""
    sources = load_sources()
    lines = ["HG008 series (alignments per sample: coverage; builds):"]
    for series_id, spec in sources["series"].items():
        alns = [r for r in records if r["kind"] == "aln" and r["series"] == series_id]
        samples = {}
        for r in alns:
            samples.setdefault(r["sample"], r["coverage"])
        builds = sorted({r["build"] for r in alns})
        have = ", ".join(f"{s} {c}x" for s, c in sorted(samples.items()))
        normal = f"; normals from {spec['normal_series']}" if spec.get("normal_series") else ""
        lines.append(
            f"  {series_id:<18} {spec['platform']} ({spec['center']}, {spec['aligner']}): "
            f"{have or 'none'}; {', '.join(builds)}{normal}"
        )
    lines.append("Caller VCFs:")
    for caller in sources["callers"]:
        vcfs = [r for r in records if r["kind"] == "vcf" and r["caller"] == caller]
        if vcfs:
            lines.append(
                f"  {caller:<18} series {vcfs[0]['series']}, T vs {vcfs[0]['normal']}; "
                f"{', '.join(sorted(r['build'] for r in vcfs))}"
            )
    lines.append(
        "Truth: small variants "
        + ", ".join(sources["truth"])
        + " (GRCh38, GRCh37); SVs "
        + ", ".join(sources["truth_sv"])
        + " (GRCh38)"
    )
    return "\n".join(lines)


def check_config(config, sources):
    """Value checks of an HG008 config that don't need the catalog."""
    sel, data = config["select"], config["data"]
    build = sel["build"]
    if build not in sources["references"]:
        raise ConfigError(
            f"select.build: {build} not available; choose from {', '.join(sources['references'])}"
        )
    for key in ("truth", "truth_sv"):
        release = data[key]
        if not release:
            continue
        spec = sources[key].get(release)
        if spec is None:
            raise ConfigError(
                f"data.{key}: release {release} unknown; choose from {', '.join(sources[key])}"
            )
        if build not in spec["builds"]:
            raise ConfigError(
                f"data.{key}: {release} is not available for {build} "
                f"(builds: {', '.join(spec['builds'])}); set data.{key}: false"
            )
    tiers = sources["truth"][data["truth"]]["tiers"] if data["truth"] else []
    if data["truth"] and data["truth_tier"] not in tiers:
        raise ConfigError(
            f"data.truth_tier: {data['truth_tier']} unknown; choose from {', '.join(tiers)}"
        )
    if data["gencode"] and build not in sources["annotation"]["gencode"]["builds"]:
        raise ConfigError(f"data.gencode: GENCODE GTFs are used for GRCh38 only, not {build}")


def select(config, records):
    """The selected units: one per tumor/normal pair x series, with VCFs and alignments.

    Units have the fields of SEQC2 units (pair, center, replicate, technology, vcfs, bams,
    bais) and: series, platform, aligner, build, library_type, samples ({T, N}: sample
    names), coverage ({T, N}), bam_names ({T, N}: file names in the dataset).
    """
    sources = load_sources()
    check_config(config, sources)
    sel, data = config["select"], config["data"]
    build = sel["build"]
    alns = [r for r in records if r["kind"] == "aln" and r["build"] == build]
    vcfs = [r for r in records if r["kind"] == "vcf" and r["build"] == build]
    # series with a tumor alignment can form pairs (ont-ne holds only normals)
    tumor_series = {r["series"] for r in alns if r["sample"] == sources["tumor"]}
    series = choose("series", sel["series"], tumor_series)
    normals = choose("normals", sel["normals"], set(sources["normals"]))
    callers = choose("callers", sel["callers"] or [], set(sources["callers"]))
    pairs = [PAIR_NAMES[n] for n in normals]
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

    def alignment(series_id, sample):
        found = [r for r in alns if r["series"] == series_id and r["sample"] == sample]
        return found[0] if found else None

    units, missing, used_callers = [], [], set()
    for series_id in series:
        spec = sources["series"][series_id]
        normal_series = spec.get("normal_series", series_id)
        for normal in normals:
            pair = PAIR_NAMES[normal]
            tumor_aln = alignment(series_id, sources["tumor"])
            normal_aln = alignment(normal_series, normal)
            if normal_aln is None:
                missing.append(
                    f"{pair} {series_id}: no {normal} alignment ({build}) in series "
                    f"{normal_series}"
                )
                continue
            unit_vcfs = []
            if data["vcf"]:
                for caller in callers:
                    found = [
                        r
                        for r in vcfs
                        if r["caller"] == caller
                        and r["series"] == series_id
                        and r["normal"] == normal
                    ]
                    if found:
                        used_callers.add(caller)
                        unit_vcfs.append({**found[0], "aligner": spec["aligner"]})
            bams, bais = {}, {}
            if pair in bam_pairs:
                for role, rec in (("T", tumor_aln), ("N", normal_aln)):
                    bams[role] = rec
                    bais[role] = rec["index"]
            units.append(
                {
                    "pair": pair,
                    "center": spec["center"],
                    "replicate": None,
                    "technology": series_id,
                    "series": series_id,
                    "platform": spec["platform"],
                    "aligner": spec["aligner"],
                    "build": build,
                    "library_type": "wgs",
                    "samples": {
                        "T": sources["samples"][sources["tumor"]],
                        "N": sources["samples"][normal],
                    },
                    "coverage": {"T": tumor_aln["coverage"], "N": normal_aln["coverage"]},
                    "bam_names": {
                        "T": f"{sources['samples'][sources['tumor']]}_{series_id}",
                        "N": f"{sources['samples'][normal]}_{normal_series}",
                    },
                    "vcfs": unit_vcfs,
                    "bams": bams,
                    "bais": bais,
                }
            )
    if data["vcf"]:
        for caller in sorted(set(callers) - used_callers):
            spec = sources["callers"][caller]
            missing.append(
                f"caller {caller}: VCFs only for series {spec['series']}, tumor vs "
                f"{spec['normal']}, builds {', '.join(spec['builds'])}"
            )
    if missing:
        raise ConfigError(
            "the selection asks for files that HG008 doesn't provide:\n  "
            + "\n  ".join(missing)
            + "\nnarrow the selection (other series, normals or callers; callers: [] or "
            "data.vcf: false for no caller VCFs)"
        )
    return units


def builder_sources(config):
    """The HG008 server description resolved for one config (build, truth tier).

    Returns the layout of sources.yaml (SEQC2) that build.Builder and build.plan read.
    """
    sources = load_sources()
    build = config["select"]["build"]
    data = config["data"]
    ref = sources["references"][build]
    resolved = {
        "name": sources["name"],
        "citation": sources["citation"],
        "base_url": sources["base_url"],
        "reference": {
            "name": ref["name"],
            "fasta": ref["fasta"],
            "fai": ref["fasta"] + ".fai",
            "gzi": ref["fasta"] + ".gzi",
        },
        "target_bed": {},
        "vcf_samples": {"tumor": None, "normal": None},
        "caller_fields": {c: spec["fields"] for c, spec in sources["callers"].items()},
        "truth": {},
        "truth_sv": {},
        "annotation": copy.deepcopy(sources["annotation"]),
    }

    def files(spec, keys, **fmt):
        """{key: {path, md5_list}} of a benchmark release for the build."""
        values = {"build_folder": spec["builds"][build], **fmt}
        md5_list = f"{spec['folder']}/{spec['md5']}"
        return {
            key: {
                "path": f"{spec['folder']}/{spec[src].format(**values)}",
                "md5_list": md5_list,
                "md5_key": spec[src].format(**values),
                **extra,
            }
            for key, (src, extra) in keys.items()
        }

    for release, spec in sources["truth"].items():
        if build in spec["builds"]:
            tier = data.get("truth_tier") or "all"
            resolved["truth"][release] = files(
                spec,
                {
                    "snv": ("vcf", {"types": "snps"}),
                    "indel": ("vcf", {"types": "indels"}),
                    "regions": ("regions", {}),
                },
                tier=tier,
            )
    for release, spec in sources["truth_sv"].items():
        if build in spec["builds"]:
            resolved["truth_sv"][release] = files(
                spec,
                {
                    "vcf": ("vcf", {}),
                    "regions": ("regions", {}),
                    "regions_clonal": ("regions_clonal", {}),
                },
            )
    giab = resolved["annotation"]["giab"]
    giab["base_url"] = giab["base_url"].format(build=build)
    giab["md5"] = giab["md5"].format(build=build)
    giab["files"] = {k: v.format(build=build) for k, v in giab["files"].items()}
    return resolved
