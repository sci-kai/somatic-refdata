"""Remote datasets (access: remote): URLs, reference choice, restrictions."""

import pytest
from test_selection import config, records

from seqc2_toolkit import build, catalog, dataset


def test_remote_facts(tmp_path):
    cfg = config(
        tmp_path,
        "access: remote\nselect:\n  technologies: [WES]\n  pairs: [IL_1]\n"
        "subset:\n  contigs: [chr21]\n",
    )
    sources = catalog.load_sources()
    facts = build.remote_facts(cfg, dataset.select(cfg, records()), sources)
    assert facts["access"] == "remote" and facts["contigs"] == ["chr21"]
    # the chr21 reference of remote_references, the original VCFs of the catalog
    chr21 = sources["remote_references"]["nfcore_chr21"]
    assert facts["resources"]["reference"] == chr21["fasta"]
    assert facts["units"][0]["vcfs"][0]["path"].endswith("WES/WES_IL_1.bwa.muTect2.vcf.gz")
    assert facts["resources"]["truth_snv"].startswith(sources["base_url"])
    root = build.write_remote(tmp_path / "ds", facts)
    assert (root / "dataset.json").exists() and "chr21" in (root / "README.md").read_text()


def test_remote_whole_reference_without_contigs(tmp_path):
    cfg = config(tmp_path, "access: remote\n")
    facts = build.remote_facts(cfg, dataset.select(cfg, records()), catalog.load_sources())
    assert facts["resources"]["reference"].endswith("GRCh38.d1.vd1.fa")
    assert facts["contigs"] is None


def test_remote_errors(tmp_path):
    sources = catalog.load_sources()
    cfg = config(tmp_path, "access: remote\nsubset:\n  contigs: [chr22]\n")
    with pytest.raises(dataset.ConfigError, match="no remote reference"):
        build.remote_facts(cfg, dataset.select(cfg, records()), sources)
    cfg = config(
        tmp_path,
        "access: remote\ndata:\n  bam: true\nsubset:\n  bam_windows: ['chr21:1-100']\n",
    )
    with pytest.raises(dataset.ConfigError, match="BAM slices"):
        build.remote_facts(cfg, dataset.select(cfg, records()), sources)


def test_remote_local_target_bed(tmp_path):
    """data.target_bed: local writes a merged BED3 of the reference contigs (no network)."""
    import gzip

    from seqc2_toolkit.cache import Cache

    sources = catalog.load_sources()
    cache = Cache(tmp_path / "cache")
    kit = f"{sources['base_url']}/{sources['target_bed']['WES']}"
    fai = sources["remote_references"]["nfcore_chr21"]["fai"]
    # files already in the cache: fetch() only records them
    cache.path(kit).parent.mkdir(parents=True)
    with gzip.open(cache.path(kit), "wt") as fh:
        fh.write("chr21\t100\t200\tx\t0\t+\nchr1\t5\t9\tx\t0\t+\nchr21\t150\t300\tx\t0\t+\n")
    cache.path(fai).parent.mkdir(parents=True)
    cache.path(fai).write_text("chr21\t46709983\t7\t60\t61\n")
    cfg = config(
        tmp_path,
        "access: remote\nselect:\n  technologies: [WES]\n  pairs: [IL_1]\n"
        "data:\n  target_bed: local\nsubset:\n  contigs: [chr21]\n",
    )
    facts = build.remote_facts(cfg, dataset.select(cfg, records()), sources, cache, tmp_path / "ds")
    rel = facts["resources"]["target_bed_wes"]
    assert rel == "bed/S07604624_Covered_human_all_v6_plus_UTR_hg38.bed.gz"
    assert facts["units"][0]["target_bed"] == rel
    # chr1 dropped (not in the chr21 reference), overlapping chr21 intervals merged
    with gzip.open(tmp_path / "ds" / rel, "rt") as fh:
        assert fh.read() == "chr21\t100\t300\n"


def test_target_bed_values(tmp_path):
    with pytest.raises(dataset.ConfigError, match="target_bed must be"):
        config(tmp_path, "data:\n  target_bed: yes-please\n")
    # local datasets build their BED anyway
    assert config(tmp_path, "data:\n  target_bed: local\n")["data"]["target_bed"] is True
