"""Build a dataset folder from the selection: download (cache), subset, index, describe.

Every output file of the dataset has a recipe: its source (URL + server md5 if known) and
the settings that shape it (contigs, BAM windows, ...). somref_manifest.json records the
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

from . import __version__, net, vcfprep
from .cache import file_state, md5sum
from .catalog import is_md5_list
from .dataset import ConfigError, as_list

log = logging.getLogger("somref")

MANIFEST = "somref_manifest.json"
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
            f"missing tools on PATH: {', '.join(missing)} " "(conda env: envs/somatic-refdata.yaml)"
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
        # atomic: an interrupted write must not leave a broken manifest
        tmp = self.path.with_name(self.path.name + ".tmp")
        tmp.write_text(json.dumps(data, indent=1, sort_keys=True))
        tmp.replace(self.path)


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
        # remove the old outputs first: they may be hard links or symlinks into the cache
        # (whole BAMs, reference), which writing in place would change for every dataset
        for r in rels:
            old = self.root / r
            if old.exists() or old.is_symlink():
                old.unlink()
        log.info("building: %s", rel)
        build(path)
        for r in rels:
            self.manifest.record(r, recipe, source)
        return path

    def md5_list(self, url):
        """{path: md5} of a server md5 list fetched through the cache; an error page or
        throttling reply instead of a list is removed from the cache and is an error."""
        path = self.cache.fetch(url)
        text = path.read_text(errors="replace")
        if not is_md5_list(text):
            path.unlink()
            raise RuntimeError(
                f"{url} is not an md5 list ({len(text)} bytes: {text[:80]!r}); "
                "the server may be busy: run the command again"
            )
        return parse_md5_text(text)

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
            # plain FASTA (data.reference_plain) for tools that can't read bgzipped FASTA
            plain_ref = self.config["data"].get("reference_plain", False)
            rel = f"reference/{ref['name']}_{'_'.join(self.contigs)}.fa"
            rel += "" if plain_ref else ".gz"
            recipe = {"source": url, "contigs": self.contigs}
            # only a set option enters the recipe: datasets without it keep their recipe hash
            if plain_ref:
                recipe["plain"] = True

            def build(path):
                plain = path if plain_ref else path.with_suffix("")
                with tempfile.TemporaryDirectory() as tmp, plain.open("w") as fh:
                    # samtools reads the remote .fai and only the needed blocks of the FASTA;
                    # it caches the remote index in the working directory (tmp)
                    run(["samtools", "faidx", url, *self.contigs], cwd=tmp, stdout=fh)
                if not plain_ref:
                    run(["bgzip", "-f", plain])
                run(["samtools", "faidx", path])

            companions = [rel + ".fai"] + ([] if plain_ref else [rel + ".gzi"])
            return self.output(rel, recipe, build, companions, url)
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
        """Somatic VCF: bgzipped (some SEQC2 VCFs are plain gzip), contig subset, indexed.

        With the vcf section of the config (sample_order, format_dp) the records are
        changed afterwards (vcfprep.py).
        """
        require_tools("bcftools", "tabix", "bgzip")
        recipe = {"source": record["url"], "md5": record["md5"], "contigs": self.contigs}
        # the vcf section exists for SEQC2 configs only
        changes = {k: v for k, v in self.config.get("vcf", {}).items() if v}
        # only a set option enters the recipe: datasets without it keep their recipe hash
        if changes:
            recipe["vcf"] = changes

        def build(path):
            src = self.cache.fetch(record["url"], record["md5"], self.verify)
            subset_vcf(src, path, self.contigs)
            if changes:
                # plain text next to the output, bgzipped and indexed again in place
                plain = path.with_suffix("")
                vcfprep.transform_file(
                    path, plain, changes.get("sample_order"), changes.get("format_dp", False)
                )
                run(["bgzip", "-f", plain])
                run(["tabix", "-f", "-p", "vcf", path])

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
                    lists[list_url] = self.md5_list(list_url)
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
        md5_url = f"{base}/{spec['md5']}"
        md5 = self.md5_list(md5_url)
        # local folder named after the build: no '@' in file names
        folder = f"annotation/giab/{Path(base).name.split('@')[0]}"
        result = {}
        for name in names:
            rel_src = spec["files"][name]
            url = f"{base}/{rel_src}"
            rel = f"{folder}/{rel_src}"
            expected = md5.get(rel_src)
            if expected is None:
                raise RuntimeError(f"{rel_src} is not in {md5_url}")
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
        check_bam_subset(self.config, units)  # before any remote work
        rel = lambda p: str(Path(p).relative_to(self.root))  # noqa: E731
        order = self.contig_order()
        # contigs: the records of every file are restricted to them already (local build)
        facts = {
            "name": self.config["name"],
            "access": "local",
            "contigs": self.contigs,
            "resources": {},
            "units": [],
        }
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
            "Created by [somatic-refdata](https://github.com/sci-kai/somatic-refdata)"
            f" {__version__} on {datetime.now().date().isoformat()}. Do not edit by hand;"
            " change the dataset config and run `somref fetch` again.",
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
            "access": self.config.get("access", "local"),
            **{k: self.config[k] for k in ("select", "data", "subset")},
            # the vcf section (sample order, FORMAT/DP) exists for SEQC2 configs only
            **({"vcf": self.config["vcf"]} if "vcf" in self.config else {}),
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


def check_bam_subset(config, units):
    """Stop if alignments would be linked whole although the config asks for a subset:
    contigs and subsampling are applied to BAMs/CRAMs only within subset.bam_windows."""
    subset = config["subset"]
    if not any(u["bams"] for u in units) or subset["bam_windows"]:
        return
    asked = [k for k in ("contigs", "bam_subsample") if subset[k]]
    if asked:
        raise ConfigError(
            f"subset.{' / subset.'.join(asked)} with data.bam but no subset.bam_windows: "
            "whole alignments would be used unchanged; give subset.bam_windows on the "
            "contigs (e.g. chr22:1-50818468) or leave the subset out"
        )


def plan(config, units, cache, sources):
    """What fetch would get: per remote source, its size and whether it is cached already.

    Returns dicts with item, url, size (bytes or None), how (download / remote read) and
    cached (ok, unverified, missing, or '-' for remote reads).
    """
    data, subset = config["data"], config["subset"]
    check_bam_subset(config, units)
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


# remote datasets (access: remote) ------------------------------------------------------


def remote_reference(sources, contigs):
    """(fasta URL, fai URL, name) of the reference for a remote dataset.

    Without subset.contigs: the SEQC2 GRCh38.d1.vd1 FASTA (3.2 GB, read whole by the
    workflow). With contigs: a smaller reference of remote_references that holds exactly
    these contigs; there is none for other contig sets (error with the choices).
    """
    if not contigs:
        ref = sources["reference"]
        base = sources["base_url"]
        return f"{base}/{ref['fasta']}", f"{base}/{ref['fai']}", ref["name"]
    for name, ref in sources.get("remote_references", {}).items():
        if sorted(ref["contigs"]) == sorted(contigs):
            return ref["fasta"], ref["fai"], name
    choices = "; ".join(
        f"{name}: {', '.join(ref['contigs'])}"
        for name, ref in sources.get("remote_references", {}).items()
    )
    raise ConfigError(
        f"subset.contigs {', '.join(contigs)}: no remote reference with exactly these contigs "
        f"(available: {choices or 'none'}); use one of them, leave subset.contigs out "
        "(whole SEQC2 reference, 3.2 GB) or access: local (the toolkit cuts the contigs)"
    )


def check_remote(config):
    """Stop if a remote dataset asks for something that needs local processing."""
    subset = config["subset"]
    problems = []
    if subset["bam_windows"] or subset["bam_subsample"]:
        problems.append(
            "subset.bam_windows / bam_subsample: BAM slices are cut by somref fetch (access: local)"
        )
    if any(config["vcf"].values()):
        problems.append(
            "vcf.sample_order / vcf.format_dp: the VCFs are changed by somref fetch (access: local)"
        )
    if config["data"]["bam"]:
        log.warning(
            "data.bam with access: remote: the workflow downloads every BAM whole (10-160 GB each)"
        )
    if problems:
        raise ConfigError("access: remote can't do this:\n  " + "\n  ".join(problems))


def local_target_bed(cache, url, fai_url, contigs, root):
    """Target BED of a remote dataset written into its folder: bed/<kit>.bed.gz (relative).

    For workflows that need a target BED with a .bed/.bed.gz name (the remote file is a
    lifted-over BED6 named *.bed6.gz): BED3 of the contigs of the reference (fai_url; only
    subset.contigs if given), sorted in reference order, overlapping intervals merged;
    gzip with mtime 0, so the same input gives the same file (stable checksums).
    """
    order = [ln.split("\t", 1)[0] for ln in cache.fetch(fai_url).read_text().splitlines() if ln]
    rank = {c: i for i, c in enumerate(order) if not contigs or c in contigs}
    intervals = []
    with gzip.open(cache.fetch(url), "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            chrom, start, end = line.split("\t")[:3]
            # only contigs of the reference (and of subset.contigs)
            if chrom in rank:
                intervals.append((rank[chrom], int(start), int(end), chrom))
    intervals.sort()
    rel = f"bed/{Path(url).name.replace('.liftover.to.hg38.bed6.gz', '_hg38.bed.gz')}"
    path = Path(root) / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as out:
        out.write("".join(f"{c}\t{s}\t{e}\n" for c, s, e in merge(intervals)).encode())
    log.info("target BED written: %s (%d intervals before merging)", rel, len(intervals))
    return rel


def remote_facts(config, units, sources, cache=None, root=None):
    """dataset.json content of a remote dataset: every file is a URL on its server.

    Nothing is changed: the VCFs are the original files of the SEQC2 server (plain gzip,
    no index, sample columns TUMOR and NORMAL, all contigs). subset.contigs only chooses
    the reference (remote_reference) and is passed on as 'contigs': the workflow has to
    restrict the records itself (e.g. bcftools view --targets). With data.target_bed:
    local, the target BED is downloaded (cache) and written into the dataset folder root
    (local_target_bed); everything else stays a URL.
    """
    check_remote(config)
    data = config["data"]
    contigs = as_list(config["subset"]["contigs"])
    facts = {"name": config["name"], "access": "remote", "contigs": contigs, "resources": {}}
    res = facts["resources"]
    if data["reference"]:
        fasta, fai, ref_name = remote_reference(sources, contigs)
        res["reference"], res["reference_fai"] = fasta, fai
        facts["reference_name"] = ref_name
    technologies = sorted({u["technology"] for u in units})
    targets = {}
    if data["target_bed"]:
        for technology in technologies:
            if technology in sources["target_bed"]:
                # original BED6 (gzip) of the kit, lifted over to hg38: all contigs
                url = f"{sources['base_url']}/{sources['target_bed'][technology]}"
                if data["target_bed"] == "local":
                    if cache is None or root is None or not data["reference"]:
                        raise ConfigError(
                            "data.target_bed: local needs the reference (its contigs) and a "
                            "dataset folder"
                        )
                    url = local_target_bed(cache, url, res["reference_fai"], contigs, root)
                targets[technology] = url
                res[f"target_bed_{technology.lower()}"] = targets[technology]
    facts["units"] = []
    for unit in units:
        entry = {k: unit[k] for k in ("pair", "center", "replicate", "technology")}
        entry["target_bed"] = targets.get(unit["technology"])
        entry["vcfs"] = [
            {"aligner": r["aligner"], "caller": r["caller"], "path": r["url"], "md5": r["md5"]}
            for r in unit["vcfs"]
        ]
        entry["bams"] = {role: r["url"] for role, r in unit["bams"].items()}
        # the server's index names (SEQC2: <name>.bai, not <name>.bam.bai)
        entry["bais"] = {role: unit["bais"][role]["url"] for role in unit["bams"]}
        facts["units"].append(entry)
    if data["truth"]:
        files = sources["truth"].get(data["truth"])
        if files is None:
            raise ConfigError(
                f"data.truth: release {data['truth']} unknown; "
                f"choose from {', '.join(sources['truth'])}"
            )
        res.update({f"truth_{k}": f"{sources['base_url']}/{v}" for k, v in files.items()})
    if data["gencode"]:
        # plain gzip, not tabix-indexed: only for workflows that read it whole
        url = sources["annotation"]["gencode"]["url"].format(release=data["gencode"])
        res["gencode_gtf"] = url
    giab = sources["annotation"]["giab"]
    for name in as_list(data["giab"]) or []:
        if name not in giab["files"]:
            raise ConfigError(f"data.giab: {name} unknown; choose from {', '.join(giab['files'])}")
        res[f"giab_{name}"] = f"{giab['base_url']}/{giab['files'][name]}"
    facts["vcf_samples"] = sources["vcf_samples"]
    facts["caller_fields"] = sources["caller_fields"]
    facts["config"] = {k: config[k] for k in ("access", "select", "data", "subset", "vcf")}
    facts["toolkit_version"] = __version__
    return facts


def write_remote(root, facts):
    """Write dataset.json and README.md of a remote dataset into root; return root."""
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    (root / "dataset.json").write_text(json.dumps(facts, indent=1))
    files = [v["path"] for u in facts["units"] for v in u["vcfs"]]
    files += [p for u in facts["units"] for p in u["bams"].values()]
    files += sorted(set(facts["resources"].values()))
    lines = [
        f"# SEQC2 dataset `{facts['name']}` (remote)",
        "",
        f"Created by somatic-refdata {__version__} on {datetime.now().date().isoformat()}."
        " Nothing was downloaded: the workflow reads these files from their servers"
        " (internet access needed). Do not edit by hand; change the dataset config and run"
        " `somref fetch` again.",
        "",
        "Contigs the workflow should restrict the records to: "
        + (", ".join(facts["contigs"]) if facts["contigs"] else "all"),
        "",
        "## Settings",
        "",
        "```json",
        json.dumps(facts["config"], indent=1),
        "```",
        "",
        "## Files",
        "",
        *[f"- {url}" for url in files],
        "",
        "Cite the SEQC2 somatic reference data: Fang et al. 2021, *Nat Biotechnol* 39:1151"
        " (doi:10.1038/s41587-021-00993-6); Xiao et al. 2021, *Nat Biotechnol* 39:1141"
        " (doi:10.1038/s41587-021-00994-5).",
    ]
    (root / "README.md").write_text("\n".join(lines) + "\n")
    return root
