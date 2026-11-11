"""Build a dataset folder from the selection: download (cache), subset, index, describe.

Every output file of the dataset has a recipe: its source (URL + server md5 if known) and
the settings that shape it (contigs, BAM windows, ...). seqc2_manifest.json records the
recipe hash and the md5 of every output; an output is only rebuilt if its recipe changed
or the file was modified. Outputs that the current config no longer produces are removed.

Dataset layout (paths relative to the dataset folder):
  reference/GRCh38.d1.vd1[_<contigs>].fa[.gz] (+ .fai, .gzi)
  bed/S07604624_Covered_human_all_v6_plus_UTR_hg38.bed     WES targets, merged BED3
  vcf/<pair>_<technology>.<aligner>.<caller>.vcf.gz (+ .tbi) bgzipped, indexed
  bam/<pair>_<technology>_<T|N>.bam (+ .bam.bai)
  truth/<release files>                                     truth VCFs (+ .tbi), regions
  annotation/gencode.v<release>.basic.annotation.gtf.gz (+ .tbi)
  annotation/giab/GRCh38/<release path>                      GIAB stratifications
  dataset.json      the dataset for workflow profiles (units, files, resources)
  checksums.md5     md5 of every file; README.md: provenance
"""

import gzip
import hashlib
import json
import logging
import os
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

from . import __version__, net
from .cache import file_state, md5sum
from .dataset import ConfigError, as_list

log = logging.getLogger("seqc2")

MANIFEST = "seqc2_manifest.json"
# unit fields of HG008 datasets copied into dataset.json (SEQC2 units don't have them)
UNIT_EXTRAS = ("series", "platform", "aligner", "build", "library_type", "samples", "coverage")
SEQC2_CITATION = (
    "Cite the SEQC2 somatic reference data: Fang et al. 2021, *Nat Biotechnol* 39:1151"
    " (doi:10.1038/s41587-021-00993-6); Xiao et al. 2021, *Nat Biotechnol* 39:1141"
    " (doi:10.1038/s41587-021-00994-5)."
)


def run(cmd, cwd=None, stdout=None):
    """Run a command (list of str), log it and stop with its stderr on failure."""
    cmd = [str(c) for c in cmd]
    log.debug("$ %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=cwd, stdout=stdout, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(cmd)}\n{result.stderr}"
        )


def require_tools(*tools):
    """Stop if a command-line tool needed for the next step is not on PATH."""
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise RuntimeError(
            f"missing tools on PATH: {', '.join(missing)} " "(conda env: envs/seqc2-toolkit.yaml)"
        )


def recipe_hash(recipe):
    """Stable hash of a recipe (dict): changes when a source or a setting changes."""
    return hashlib.sha256(json.dumps(recipe, sort_keys=True).encode()).hexdigest()[:16]


def link(source, dest):
    """Hard link a cached file into the dataset (symlink if on another file system)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    try:
        os.link(source, dest)
    except OSError:
        os.symlink(Path(source).resolve(), dest)
        log.info("symlinked (cache on another file system): %s", dest.name)


class Manifest:
    """The outputs of a dataset folder with their recipes and md5 sums."""

    def __init__(self, root):
        self.root = root
        self.path = root / MANIFEST
        self.entries = json.loads(self.path.read_text())["files"] if self.path.exists() else {}
        self.produced = set()

    def up_to_date(self, rel, recipe, verify):
        """True if rel exists, was built by this recipe and is unchanged."""
        entry, path = self.entries.get(rel), self.root / rel
        if entry is None or not path.exists() or entry["recipe"] != recipe_hash(recipe):
            return False
        if verify or [entry["size"], entry["mtime_ns"]] != list(file_state(path)):
            return md5sum(path) == entry["md5"]
        return True

    def record(self, rel, recipe, source=None):
        """Remember an output (after it was built)."""
        path = self.root / rel
        size, mtime = file_state(path)
        self.entries[rel] = {
            "recipe": recipe_hash(recipe),
            "md5": md5sum(path),
            "size": size,
            "mtime_ns": mtime,
            "source": source or "",
        }
        self.produced.add(rel)
        self.save()

    def keep(self, rel):
        """Mark an up-to-date output as part of the current dataset."""
        self.produced.add(rel)

    def prune(self):
        """Remove outputs of earlier configs that the current config doesn't produce."""
        for rel in sorted(set(self.entries) - self.produced):
            path = self.root / rel
            if path.exists() or path.is_symlink():
                path.unlink()
            log.info("removed (no longer selected): %s", rel)
            del self.entries[rel]
        self.save()

    def save(self):
        """Write the manifest (after every output, so an interrupted run keeps it)."""
        data = {"toolkit_version": __version__, "files": self.entries}
        self.path.write_text(json.dumps(data, indent=1, sort_keys=True))


class Builder:
    """Builds one dataset folder; every step is skipped if its outputs are up to date."""

    def __init__(self, config, sources, cache, outdir, verify=False):
        self.config, self.sources, self.cache = config, sources, cache
        self.subset = config["subset"]
        self.contigs = as_list(self.subset["contigs"])
        self.root = Path(outdir).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.verify = verify
        self.manifest = Manifest(self.root)
        self.base = sources["base_url"]

    def output(self, rel, recipe, build, companions=(), source=None):
        """Build rel (and its companion files, e.g. indices) unless up to date.

        build(path) writes the file at path (its companions next to it).
        """
        rels = [rel, *companions]
        if all(self.manifest.up_to_date(r, recipe, self.verify) for r in rels):
            for r in rels:
                self.manifest.keep(r)
            log.info("up to date: %s", rel)
            return self.root / rel
        path = self.root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        log.info("building: %s", rel)
        build(path)
        for r in rels:
            self.manifest.record(r, recipe, source)
        return path

    # reference -------------------------------------------------------------------

    def contig_order(self):
        """Contig names of the reference in .fai order (sort order of BEDs and GTFs)."""
        fai = self.cache.fetch(f"{self.base}/{self.sources['reference']['fai']}")
        order = [line.split("\t", 1)[0] for line in fai.read_text().splitlines() if line]
        if self.contigs:
            unknown = [c for c in self.contigs if c not in order]
            if unknown:
                raise ConfigError(
                    f"subset.contigs: {', '.join(unknown)} not in "
                    f"{self.sources['reference']['name']}"
                )
        return order

    def reference(self):
        """FASTA + index: contig subset read remotely, or the whole FASTA from the cache."""
        ref = self.sources["reference"]
        url = f"{self.base}/{ref['fasta']}"
        if self.contigs:
            require_tools("samtools", "bgzip")
            rel = f"reference/{ref['name']}_{'_'.join(self.contigs)}.fa.gz"
            recipe = {"source": url, "contigs": self.contigs}

            def build(path):
                plain = path.with_suffix("")
                with tempfile.TemporaryDirectory() as tmp, plain.open("w") as fh:
                    # samtools reads the remote .fai and only the needed blocks of the FASTA;
                    # it caches the remote index in the working directory (tmp)
                    run(["samtools", "faidx", url, *self.contigs], cwd=tmp, stdout=fh)
                run(["bgzip", "-f", plain])
                run(["samtools", "faidx", path])

            return self.output(rel, recipe, build, [rel + ".fai", rel + ".gzi"], url)
        # whole FASTA (3.2 GB) and its index (+ .gzi of a bgzipped FASTA) from the cache
        rel = f"reference/{Path(ref['fasta']).name}"
        recipe = {"source": url}
        companions = [rel + ".fai"] + ([rel + ".gzi"] if ref.get("gzi") else [])

        def build(path):
            link(self.cache.fetch(url), path)
            link(self.cache.fetch(f"{self.base}/{ref['fai']}"), Path(str(path) + ".fai"))
            if ref.get("gzi"):
                link(self.cache.fetch(f"{self.base}/{ref['gzi']}"), Path(str(path) + ".gzi"))

        return self.output(rel, recipe, build, companions, url)

    # regions ---------------------------------------------------------------------

    def target_bed(self, technology, order):
        """Target BED of a technology: BED3, contig subset, sorted in reference order, merged."""
        url = f"{self.base}/{self.sources['target_bed'][technology]}"
        name = Path(url).name.replace(".liftover.to.hg38.bed6.gz", "_hg38.bed")
        rel = f"bed/{name}"
        recipe = {"source": url, "contigs": self.contigs}

        def build(path):
            src = self.cache.fetch(url)
            rank = {c: i for i, c in enumerate(order)}
            intervals = []
            with gzip.open(src, "rt") as fh:
                for line in fh:
                    if line.startswith(("#", "track", "browser")) or not line.strip():
                        continue
                    chrom, start, end = line.split("\t")[:3]
                    if chrom in rank and (not self.contigs or chrom in self.contigs):
                        intervals.append((rank[chrom], int(start), int(end), chrom))
            intervals.sort()
            path.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in merge(intervals)))

        return self.output(rel, recipe, build, source=url)

    # VCFs ------------------------------------------------------------------------

    def vcf(self, record, rel):
        """Somatic VCF: bgzipped (some SEQC2 VCFs are plain gzip), contig subset, indexed."""
        require_tools("bcftools", "tabix")
        recipe = {"source": record["url"], "md5": record["md5"], "contigs": self.contigs}

        def build(path):
            src = self.cache.fetch(record["url"], record["md5"], self.verify)
            subset_vcf(src, path, self.contigs)

        return self.output(rel, recipe, build, [rel + ".tbi"], record["url"])

    # BAMs ------------------------------------------------------------------------

    def bam_regions(self, target_bed):
        """BED of the regions read from remote BAMs (windows, optionally only targets)."""
        windows = []
        for window in as_list(self.subset["bam_windows"]):
            chrom, span = window.split(":")
            start, end = (int(x.replace(",", "")) for x in span.split("-"))
            if self.contigs and chrom not in self.contigs:
                raise ConfigError(f"subset.bam_windows: {window} is not on subset.contigs")
            windows.append((chrom, start - 1, end))
        if not (self.subset["bam_targets"] and target_bed):
            regions = windows
        else:
            # targets +/- padding, clipped to the windows
            pad = int(self.subset["bam_padding"])
            regions = []
            for line in Path(target_bed).read_text().splitlines():
                chrom, start, end = line.split("\t")[:3]
                for wchrom, wstart, wend in windows:
                    s, e = max(int(start) - pad, wstart), min(int(end) + pad, wend)
                    if chrom == wchrom and s < e:
                        regions.append((chrom, s, e))
        rel = "bam/bam_regions.bed"
        text = "".join(
            f"{c}\t{s}\t{e}\n" for c, s, e in merge([(0, s, e, c) for c, s, e in sorted(regions)])
        )
        recipe = {"regions": text}
        return self.output(rel, recipe, lambda path: path.write_text(text))

    def bam(self, record, index, name, regions):
        """Alignment: a slice of the remote BAM/CRAM (regions, written as BAM) or the whole
        BAM/CRAM from the cache. name: dataset path without extension (bam/<name>).

        CRAMs are decoded with the reference of the source (remote FASTA, read by blocks).
        """
        cram = record.get("format") == "cram"
        ref_url = f"{self.base}/{self.sources['reference']['fasta']}"
        recipe = {"source": record["url"], "md5": record["md5"]}
        if regions is not None:
            require_tools("samtools")
            rel, index_ext = f"{name}.bam", ".bai"
            recipe["regions"] = md5sum(regions)
            recipe["subsample"] = self.subset["bam_subsample"]

            def build(path):
                fraction = self.subset["bam_subsample"]
                subsample = ["--subsample", fraction, "--subsample-seed", 42] if fraction else []
                reference = ["-T", ref_url] if cram else []
                with tempfile.TemporaryDirectory() as tmp:
                    # -M: multi-region iterator (each read once); -X: index given explicitly
                    # (SEQC2 names it .bai, not .bam.bai); samtools caches it in cwd (tmp)
                    run(
                        [
                            "samtools",
                            "view",
                            "-b",
                            "-M",
                            "-L",
                            regions,
                            *subsample,
                            *reference,
                            "-o",
                            path,
                            "-X",
                            record["url"],
                            index["url"],
                        ],
                        cwd=tmp,
                    )
                run(["samtools", "index", path])

        else:
            rel, index_ext = (f"{name}.cram", ".crai") if cram else (f"{name}.bam", ".bai")

            def build(path):
                link(self.cache.fetch(record["url"], record["md5"], self.verify), path)
                link(
                    self.cache.fetch(index["url"], index["md5"], self.verify),
                    Path(str(path) + index_ext),
                )

        return self.output(rel, recipe, build, [rel + index_ext], record["url"])

    # truth and annotation --------------------------------------------------------

    def truth(self, release, kind="truth"):
        """Truth set (benchmark) files: VCFs (subset, indexed) and region BEDs (subset).

        kind: truth (SEQC2 / HG008 small variants) or truth_sv (HG008 SVs). An entry of
        the release is a path (SEQC2) or {path, md5_list, md5_key, types} (HG008: md5 from
        the benchmark's list; types: snps / indels keeps only these from a combined VCF).
        """
        files = self.sources[kind].get(release)
        if files is None:
            raise ConfigError(
                f"data.{kind}: release {release} unknown; "
                f"choose from {', '.join(self.sources[kind])}"
            )
        lists = {}
        result = {}
        for key, entry in files.items():
            entry = {"path": entry} if isinstance(entry, str) else entry
            url = f"{self.base}/{entry['path']}"
            expected = None
            if entry.get("md5_list"):
                list_url = f"{self.base}/{entry['md5_list']}"
                if list_url not in lists:
                    lists[list_url] = parse_md5_text(self.cache.fetch(list_url).read_text())
                expected = lists[list_url].get(entry["md5_key"])
                if expected is None:
                    raise RuntimeError(f"{entry['md5_key']} is not in {list_url}")
            name = Path(url).name
            types = entry.get("types")
            if types:
                name = name.replace(".vcf.gz", f"_{types}.vcf.gz")
            rel = f"{kind}/{name}"
            # md5 / types only when set: SEQC2 recipes (and their outputs) stay unchanged
            recipe = {"source": url, "contigs": self.contigs}
            recipe.update({k: v for k, v in (("md5", expected), ("types", types)) if v})
            if name.endswith(".vcf.gz"):

                def build(path, url=url, expected=expected, types=types):
                    require_tools("bcftools", "tabix")
                    src = self.cache.fetch(url, expected, self.verify)
                    subset_vcf(src, path, self.contigs, types)

                result[key] = self.output(rel, recipe, build, [rel + ".tbi"], url)
            else:

                def build(path, url=url, expected=expected):
                    src = self.cache.fetch(url, expected, self.verify)
                    lines = src.read_text().splitlines(keepends=True)
                    path.write_text(
                        "".join(
                            ln
                            for ln in lines
                            if not self.contigs or ln.split("\t", 1)[0] in self.contigs
                        )
                    )

                result[key] = self.output(rel, recipe, build, source=url)
        return result

    def gencode(self, release, order):
        """GENCODE basic GTF (VEP --gtf): contig subset, sorted, bgzipped, tabix-indexed."""
        require_tools("bgzip", "tabix")
        url = self.sources["annotation"]["gencode"]["url"].format(release=release)
        rel = f"annotation/{Path(url).name}"
        recipe = {"source": url, "contigs": self.contigs}

        def build(path):
            src = self.cache.fetch(url)
            rank = {c: i for i, c in enumerate(order)}
            records = []
            with gzip.open(src, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    fields = line.split("\t", 5)
                    if fields[0] in rank and (not self.contigs or fields[0] in self.contigs):
                        records.append((rank[fields[0]], int(fields[3]), int(fields[4]), line))
            records.sort(key=lambda r: r[:3])
            plain = path.with_suffix("")
            plain.write_text("".join(r[3] for r in records))
            run(["bgzip", "-f", plain])
            run(["tabix", "-f", "-p", "gff", path])

        return self.output(rel, recipe, build, [rel + ".tbi"], url)

    def giab(self, names):
        """GIAB stratifications (md5-checked), contig subset, in the release layout."""
        spec = self.sources["annotation"]["giab"]
        unknown = [n for n in names if n not in spec["files"]]
        if unknown:
            raise ConfigError(
                f"data.giab: {', '.join(unknown)} unknown; choose from {', '.join(spec['files'])}"
            )
        base = spec["base_url"]
        md5 = {}
        md5_url = f"{base}/{spec['md5']}"
        for line in self.cache.fetch(md5_url).read_text().splitlines():
            parts = line.split()
            if len(parts) == 2:
                md5[parts[1].removeprefix("./")] = parts[0]
        # local folder named after the build: no '@' in file names
        folder = f"annotation/giab/{Path(base).name.split('@')[0]}"
        result = {}
        for name in names:
            rel_src = spec["files"][name]
            url = f"{base}/{rel_src}"
            rel = f"{folder}/{rel_src}"
            expected = md5.get(rel_src)
            recipe = {"source": url, "md5": expected, "contigs": self.contigs}

            def build(path, url=url, expected=expected):
                src = self.cache.fetch(url, expected, self.verify)
                # mtime=0: the same content gives the same gzip file (stable checksums)
                with gzip.open(src, "rt") as fin, open(path, "wb") as raw:
                    with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fout:
                        for line in fin:
                            if line.startswith("#"):
                                continue
                            if not self.contigs or line.split("\t", 1)[0] in self.contigs:
                                fout.write(line.encode())

            result[name] = self.output(rel, recipe, build, source=url)
        return result

    # whole dataset -----------------------------------------------------------------

    def build(self, units):
        """Build every output of the selection; write dataset.json, checksums and README."""
        data = self.config["data"]
        rel = lambda p: str(Path(p).relative_to(self.root))  # noqa: E731
        order = self.contig_order()
        facts = {"name": self.config["name"], "resources": {}, "units": []}
        res = facts["resources"]
        if data["reference"]:
            fasta = self.reference()
            res["reference"] = rel(fasta)
            res["reference_fai"] = rel(fasta) + ".fai"
            if Path(str(fasta) + ".gzi").exists():
                res["reference_gzi"] = rel(fasta) + ".gzi"
        technologies = sorted({u["technology"] for u in units})
        targets = {}
        if data.get("target_bed"):
            for technology in technologies:
                if technology in self.sources["target_bed"]:
                    targets[technology] = self.target_bed(technology, order)
                    res[f"target_bed_{technology.lower()}"] = rel(targets[technology])
        regions = None
        if any(u["bams"] for u in units) and self.subset["bam_windows"]:
            regions = self.bam_regions(targets.get("WES"))
        for unit in units:
            tech = unit["technology"].lower()
            entry = {k: unit[k] for k in ("pair", "center", "replicate", "technology")}
            # HG008 units: series, platform, aligner, build, library type, samples, coverage
            entry.update({k: unit[k] for k in UNIT_EXTRAS if k in unit})
            entry["target_bed"] = (
                rel(targets[unit["technology"]]) if unit["technology"] in targets else None
            )
            entry["vcfs"] = []
            for record in unit["vcfs"]:
                out = self.vcf(
                    record,
                    f"vcf/{unit['pair']}_{tech}.{record['aligner']}.{record['caller']}.vcf.gz",
                )
                vcf_entry = {
                    "aligner": record["aligner"],
                    "caller": record["caller"],
                    "path": rel(out),
                }
                # sample columns of this VCF (HG008: differ per caller; SEQC2: dataset-wide)
                if record.get("samples"):
                    vcf_entry["tumor_id"] = record["samples"]["tumor"]
                    vcf_entry["normal_id"] = record["samples"]["normal"]
                entry["vcfs"].append(vcf_entry)
            entry["bams"] = {}
            for role, record in unit["bams"].items():
                name = unit.get("bam_names", {}).get(role, f"{unit['pair']}_{tech}_{role}")
                out = self.bam(record, unit["bais"][role], f"bam/{name}", regions)
                entry["bams"][role] = rel(out)
            facts["units"].append(entry)
        if data["truth"]:
            res.update({f"truth_{k}": rel(v) for k, v in self.truth(data["truth"]).items()})
        if data.get("truth_sv"):
            res.update(
                {
                    f"truth_sv_{k}".replace("truth_sv_vcf", "truth_sv"): rel(v)
                    for k, v in self.truth(data["truth_sv"], "truth_sv").items()
                }
            )
        if data["gencode"]:
            res["gencode_gtf"] = rel(self.gencode(data["gencode"], order))
        if data["giab"]:
            res.update({f"giab_{k}": rel(v) for k, v in self.giab(as_list(data["giab"])).items()})
        self.manifest.prune()
        facts["source"] = self.config.get("source", "seqc2")
        facts["source_name"] = self.sources.get("name", "SEQC2")
        facts["vcf_samples"] = self.sources["vcf_samples"]
        facts["caller_fields"] = self.sources["caller_fields"]
        facts["config"] = self.settings()
        facts["toolkit_version"] = __version__
        (self.root / "dataset.json").write_text(json.dumps(facts, indent=1))
        self.write_provenance()
        return facts

    def write_provenance(self):
        """checksums.md5 (md5 of every output) and README.md (sources, settings)."""
        entries = self.manifest.entries
        (self.root / "checksums.md5").write_text(
            "".join(f"{entries[r]['md5']}  {r}\n" for r in sorted(entries))
        )
        sources = sorted(
            {(r.split("/")[0], e["source"]) for r, e in entries.items() if e["source"]}
        )
        config = json.dumps(self.settings(), indent=1)
        lines = [
            f"# {self.sources.get('name', 'SEQC2')} dataset `{self.config['name']}`",
            "",
            "Created by [seqc2-toolkit](https://github.com/core-unit-bioinformatics/seqc2-toolkit)"
            f" {__version__} on {datetime.now().date().isoformat()}. Do not edit by hand;"
            " change the dataset config and run `seqc2 fetch` again.",
            "",
            "## Settings",
            "",
            "```json",
            config,
            "```",
            "",
            "## Sources",
            "",
            "| Folder | Source |",
            "|---|---|",
            *[f"| `{folder}/` | {url} |" for folder, url in sources],
            "",
            self.sources.get("citation", SEQC2_CITATION),
        ]
        (self.root / "README.md").write_text("\n".join(lines) + "\n")

    def settings(self):
        """The config sections that shape the dataset (written to dataset.json, README)."""
        return {
            "source": self.config.get("source", "seqc2"),
            **{k: self.config[k] for k in ("select", "data", "subset")},
        }


def merge(intervals):
    """Merge sorted (rank, start, end, chrom) intervals: overlapping/adjacent ones joined."""
    merged = []
    for _rank, start, end, chrom in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    return merged


def subset_vcf(src, dest, contigs, types=None):
    """bgzipped, indexed copy of a VCF, restricted to contigs (all if None) and variant
    types (bcftools -v: snps, indels; all if None).

    Streams with -t (targets), not -r: some SEQC2 VCFs are plain gzip without index.
    """
    cmd = ["bcftools", "view", "--no-version", "-Oz", "-o", dest]
    if contigs:
        cmd += ["-t", ",".join(contigs)]
    if types:
        cmd += ["-v", types]
    run([*cmd, src])
    run(["tabix", "-f", "-p", "vcf", dest])


def parse_md5_text(text):
    """{relative path: md5} of an md5 list ('<md5>  <path>'; './x', '/x', 'x' the same)."""
    result = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            result[parts[1].removeprefix("./").lstrip("/")] = parts[0]
    return result


def plan(config, units, cache, sources):
    """What fetch would get: per remote source, its size and whether it is cached already.

    Returns dicts with item, url, size (bytes or None), how (download / remote read) and
    cached (ok, unverified, missing, or '-' for remote reads).
    """
    data, subset = config["data"], config["subset"]
    base = sources["base_url"]
    contigs = as_list(subset["contigs"])
    rows = []
    if data["reference"]:
        url = f"{base}/{sources['reference']['fasta']}"
        rows.append(
            ("reference FASTA", url, None, "remote read (contigs)" if contigs else "download")
        )
    if data.get("target_bed") and any(u["technology"] in sources["target_bed"] for u in units):
        for technology in sorted({u["technology"] for u in units} & set(sources["target_bed"])):
            rows.append(
                (
                    f"{technology} target BED",
                    f"{base}/{sources['target_bed'][technology]}",
                    None,
                    "download",
                )
            )
    for unit in units:
        label = f"{unit['pair']} {unit['technology']}"
        for record in unit["vcfs"]:
            rows.append(
                (
                    f"{label} VCF {record['aligner']}/{record['caller']}",
                    record["url"],
                    record["md5"],
                    "download",
                )
            )
        for role, record in unit["bams"].items():
            how = "remote read (bam_windows)" if subset["bam_windows"] else "download"
            rows.append((f"{label} BAM {role}", record["url"], record["md5"], how))
            if not subset["bam_windows"]:
                index = unit["bais"][role]
                rows.append((f"{label} BAM index {role}", index["url"], index["md5"], "download"))
    for kind in ("truth", "truth_sv"):
        release = data.get(kind)
        if not release or release not in sources.get(kind, {}):
            continue
        seen = set()
        for key, entry in sources[kind][release].items():
            path = entry if isinstance(entry, str) else entry["path"]
            if path in seen:
                continue  # HG008: SNVs and indels come from one VCF
            seen.add(path)
            rows.append((f"{kind} {release} {key}", f"{base}/{path}", None, "download"))
    if data["gencode"]:
        url = sources["annotation"]["gencode"]["url"].format(release=data["gencode"])
        rows.append((f"GENCODE {data['gencode']} GTF", url, None, "download"))
    giab = sources["annotation"]["giab"]
    for name in as_list(data["giab"]) or []:
        if name in giab["files"]:
            rows.append(
                (f"GIAB {name}", f"{giab['base_url']}/{giab['files'][name]}", None, "download")
            )
    result, seen = [], set()
    for item, url, md5, how in rows:
        # HG008: the tumor alignment is shared by the pairs T-ND and T-NP
        if url in seen:
            continue
        seen.add(url)
        remote = how != "download"
        size = None if remote else net.head(url)[0]
        cached = "-" if remote else cache.status(url, md5)
        result.append({"item": item, "url": url, "size": size, "how": how, "cached": cached})
    return result
