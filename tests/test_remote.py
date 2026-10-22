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
