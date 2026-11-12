"""Catalog of the SEQC2 files: built from the server's md5 lists (sources_seqc2.yaml).

Every VCF, BAM and BAM index of the paired tumor/normal series (WGS, WES) becomes one
record with the fields parsed from its file name:

  kind        vcf, bam or bai
  technology  WGS or WES
  center      sequencing center (EA, FD, IL, LL, NC, NS, NV)
  replicate   replicate number (int)
  pair        <center>_<replicate>, e.g. IL_1: one tumor/normal "patient"
  role        T / N (BAMs only; a VCF holds both)
  aligner     bwa, bowtie, novo (BAMs: bwa only)
  caller      muTect2, strelka, somaticSniper (VCFs only)
  url, md5    download location and md5 from the server's list

Files of the other series (FFPE FFG/FFX, library prep LBP, purity SPP, AmpliSeq) are
listed by the server but not part of the catalog yet.
"""

import logging
import re
from importlib import resources
from pathlib import Path

import yaml

from . import net

log = logging.getLogger("somref")


SOURCES = {"seqc2": "sources_seqc2.yaml", "hg008": "sources_hg008.yaml"}


def load_sources(source="seqc2"):
    """The server description of a source (YAML file shipped with the package)."""
    if source not in SOURCES:
        raise ValueError(f"unknown source {source}; choose from {', '.join(SOURCES)}")
    text = resources.files("somatic_refdata").joinpath(SOURCES[source]).read_text()
    return yaml.safe_load(text)


def builder_sources(config):
    """The server description as build.Builder and build.plan read it, for one config."""
    if config.get("source", "seqc2") == "hg008":
        from . import hg008

        return hg008.builder_sources(config)
    return load_sources()


MD5_LINE = re.compile(r"^[0-9a-fA-F]{32}\s+\S")


def is_md5_list(text):
    """True if text has at least one '<md5> <path>' line (not an error page)."""
    return any(MD5_LINE.match(line) for line in text.splitlines())


def md5_list(cache, url, refresh=False):
    """Lines of a server md5 list; cached in the cache folder (used offline as fallback).

    A reply without any md5 line (an error page, the server's throttling reply) counts
    as a failed fetch: it is neither cached nor used.
    """
    path = cache.path(url)
    if refresh or not path.exists():
        text = net.get_text(url)
        if text is not None and not is_md5_list(text):
            log.warning("%s is not an md5 list (%d bytes): %s", url, len(text), text[:80])
            text = None
        if text is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(path.name + ".tmp")
            tmp.write_text(text)
            tmp.replace(path)
        elif not path.exists():
            raise RuntimeError(f"{url} can't be fetched and is not in the cache {cache.root}")
        else:
            log.warning("using the cached copy of %s", url)
    return path.read_text().splitlines()


def parse_md5_list(lines, kind, pattern, folder, base_url):
    """Catalog records of one md5 list ('<md5>  <path>' per line); other files are skipped."""
    regex = re.compile(pattern)
    records, skipped = [], 0
    for line in lines:
        parts = line.split()
        if len(parts) != 2:
            continue
        md5, rel = parts[0], parts[1].removeprefix("./")
        match = regex.match(rel)
        if not match:
            skipped += 1
            continue
        fields = match.groupdict()
        fields["replicate"] = int(fields["replicate"])
        # the folder may depend on the fields (BAMs: data/<technology>/)
        remote = f"{folder.format(**fields)}/{rel}"
        records.append(
            {
                "kind": kind,
                "pair": f"{fields['center']}_{fields['replicate']}",
                **fields,
                "name": Path(rel).name,
                "url": f"{base_url}/{remote}",
                "md5": md5,
            }
        )
    log.debug(
        "%s list: %d files in the catalog, %d of other series skipped", kind, len(records), skipped
    )
    return records


def build(cache, refresh=False, source="seqc2"):
    """All catalog records of a source (SEQC2: VCFs, BAMs, BAM indices of WGS/WES)."""
    if source == "hg008":
        from . import hg008

        return hg008.build_catalog(cache, refresh)
    sources = load_sources()
    base = sources["base_url"]
    records = []
    for kind, spec in sources["md5_lists"].items():
        lines = md5_list(cache, f"{base}/{spec['list']}", refresh)
        records += parse_md5_list(lines, kind, sources["patterns"][kind], spec["folder"], base)
    return records


def summary(records, source="seqc2"):
    """Text overview of the catalog: what exists per technology."""
    if source == "hg008":
        from . import hg008

        return hg008.summary(records)
    lines = []
    for technology in sorted({r["technology"] for r in records}):
        tech = [r for r in records if r["technology"] == technology]
        vcfs = [r for r in tech if r["kind"] == "vcf"]
        bams = [r for r in tech if r["kind"] == "bam"]
        pairs = sorted(
            {r["pair"] for r in tech}, key=lambda p: (p.split("_")[0], int(p.split("_")[1]))
        )
        lines.append(f"{technology}: {len(pairs)} pairs: {', '.join(pairs)}")
        lines.append(
            f"  VCFs: {len(vcfs)}; aligners: {', '.join(sorted({r['aligner'] for r in vcfs}))}; "
            f"callers: {', '.join(sorted({r['caller'] for r in vcfs}))}"
        )
        lines.append(
            f"  BAMs: {len(bams)} (tumor + normal, aligner: "
            f"{', '.join(sorted({r['aligner'] for r in bams}))})"
        )
    return "\n".join(lines)
