"""Optional VCF changes of local datasets: sample order, FORMAT/DP from AD."""

import gzip

import pytest
from test_selection import config, records

from seqc2_toolkit import build, catalog, dataset, vcfprep

HEADER = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL",
]


def test_order_and_dp():
    lines = HEADER + [
        "chr1\t5\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:7,3\t0/0:9,0",
        # existing DP is kept; a missing AD value gives DP '.'
        "chr1\t6\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP\t0/1:1,1:5\t0/0:2,0:4",
        "chr1\t7\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:.,3\t0/0",
    ]
    out, stats = vcfprep.transform_lines(lines, ["NORMAL", "TUMOR"], True)
    assert out[1].startswith("##FORMAT=<ID=DP,")
    assert out[2].startswith("##seqc2_toolkit_vcf=") and "sample order NORMAL TUMOR" in out[2]
    assert out[3].endswith("NORMAL\tTUMOR")
    assert out[4].split("\t")[8:] == ["GT:AD:DP", "0/0:9,0:9", "0/1:7,3:10"]
    assert out[5].split("\t")[8:] == ["GT:AD:DP", "0/0:2,0:4", "0/1:1,1:5"]
    assert out[6].split("\t")[8:] == ["GT:AD:DP", "0/0:.:.", "0/1:.,3:."]
    assert stats == {"records": 3, "dp_added": 2, "dp_existing": 1, "ad_missing": 0}


def test_only_order_keeps_format():
    lines = HEADER + ["chr1\t5\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:7,3\t0/0:9,0"]
    out, _ = vcfprep.transform_lines(lines, ["NORMAL", "TUMOR"], False)
    assert not any(line.startswith("##FORMAT=<ID=DP") for line in out)
    assert out[-1].split("\t")[8:] == ["GT:AD", "0/0:9,0", "0/1:7,3"]


def test_wrong_sample_names():
    with pytest.raises(RuntimeError, match="don't match vcf.sample_order"):
        vcfprep.transform_lines(HEADER, ["N", "T"], False)


def test_transform_file(tmp_path):
    src = tmp_path / "in.vcf.gz"
    with gzip.open(src, "wt") as fh:
        fh.write("\n".join(HEADER + ["chr1\t5\t.\tA\tG\t.\tPASS\t.\tAD\t7,3\t9,0"]) + "\n")
    stats = vcfprep.transform_file(src, tmp_path / "out.vcf", None, True)
    assert stats["dp_added"] == 1
    assert (tmp_path / "out.vcf").read_text().splitlines()[-1].endswith("AD:DP\t7,3:10\t9,0:9")


def test_vcf_options_need_local_access(tmp_path):
    cfg = config(tmp_path, "access: remote\nvcf:\n  format_dp: true\n")
    with pytest.raises(dataset.ConfigError, match="access: local"):
        build.remote_facts(cfg, dataset.select(cfg, records()), catalog.load_sources())
