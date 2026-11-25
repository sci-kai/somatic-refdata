"""Optional VCF changes of local datasets: sample order, FORMAT/DP from AD."""

import gzip
import shutil

import pytest
from test_selection import config, records

from somatic_refdata import build, catalog, dataset, vcfprep

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
    assert out[2].startswith("##somatic_refdata_vcf=") and "sample order NORMAL TUMOR" in out[2]
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


# header repair and sort check (SEQC2 Strelka VCFs: SNVs + indels concatenated) ------------

STRELKA_LIKE = [
    "##fileformat=VCFv4.1",
    "##contig=<ID=chr1>",
    "##contig=<ID=chr2>",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR",
    # SNV block
    "chr1\t10\t.\tA\tG\t.\tPASS\tDP=5\tDP\t5\t6",
    "chr2\t20\t.\tA\tG\t.\tPASS\tDP=7\tDP\t7\t8",
    # indel block: chr1 again, keys not in the header
    "chr1\t15\t.\tAT\tA\t.\tPASS\tQSI=7;RU=T;OVERLAP;SOR=1.5\tDP:TAR:DP50\t9:1,2:3.5\t9:3,4:.",
]


def write_vcf(path, lines):
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def test_scan_vcf_finds_undefined_keys_and_concatenated_blocks(tmp_path):
    found = vcfprep.scan_vcf(write_vcf(tmp_path / "x.vcf.gz", STRELKA_LIKE))
    assert found["unsorted"] and found["records"] == 3
    types = {vcfprep.key_of(m): m.split(",")[1:3] for m in found["missing"]}
    assert types == {
        "INFO/QSI": ["Number=1", "Type=Integer"],
        "INFO/RU": ["Number=1", "Type=String"],
        "INFO/OVERLAP": ["Number=0", "Type=Flag"],
        "INFO/SOR": ["Number=1", "Type=Float"],
        "FORMAT/TAR": ["Number=.", "Type=Integer"],
        "FORMAT/DP50": ["Number=1", "Type=Float"],
    }
    text = vcfprep.header_with(found["header"], found["missing"]).splitlines()
    assert text[-1].startswith("#CHROM") and text[-2].startswith("##FORMAT=<ID=DP50,")


def test_scan_vcf_of_a_valid_vcf_finds_nothing(tmp_path):
    lines = STRELKA_LIKE[:6] + [
        "chr1\t10\t.\tA\tG\t.\tPASS\tDP=5\tDP\t5\t6",
        "chr1\t10\t.\tA\tT\t.\tPASS\tDP=5\tDP\t5\t6",
        "chr2\t3\t.\tA\tG\t.\tPASS\t.\tDP\t7\t8",
    ]
    found = vcfprep.scan_vcf(write_vcf(tmp_path / "x.vcf.gz", lines))
    assert not found["unsorted"] and found["missing"] == [] and found["records"] == 3
    # positions going back within a chromosome
    lines[-1] = "chr1\t3\t.\tA\tG\t.\tPASS\t.\tDP\t7\t8"
    assert vcfprep.scan_vcf(write_vcf(tmp_path / "y.vcf.gz", lines))["unsorted"]


def test_mixed_flag_and_value_is_a_string():
    table = {}
    vcfprep._note_value(table, "X", None)
    vcfprep._note_value(table, "X", "3")
    assert vcfprep._definition("INFO", "X", table["X"]).startswith(
        "##INFO=<ID=X,Number=.,Type=String"
    )


@pytest.mark.skipif(
    not all(shutil.which(t) for t in ("bcftools", "tabix", "bgzip")),
    reason="needs bcftools, tabix and bgzip",
)
def test_subset_vcf_repairs_and_sorts_only_when_needed(tmp_path):
    src = write_vcf(tmp_path / "src.vcf.gz", STRELKA_LIKE)
    out = tmp_path / "out"
    out.mkdir()
    build.subset_vcf(src, out / "x.vcf.gz", None)
    # no temporary folder left behind
    assert sorted(p.name for p in out.iterdir()) == ["x.vcf.gz", "x.vcf.gz.tbi"]
    with gzip.open(out / "x.vcf.gz", "rt") as fh:
        lines = fh.read().splitlines()
    body = [line.split("\t")[:2] for line in lines if not line.startswith("#")]
    assert body == [["chr1", "10"], ["chr1", "15"], ["chr2", "20"]]
    assert any(line.startswith("##INFO=<ID=QSI,Number=1,Type=Integer") for line in lines)
    # a valid source is not changed: no added definitions, records as they were
    valid = write_vcf(tmp_path / "valid.vcf.gz", STRELKA_LIKE[:7])
    build.subset_vcf(valid, out / "v.vcf.gz", None)
    with gzip.open(out / "v.vcf.gz", "rt") as fh:
        text = fh.read()
    assert "somatic-refdata" not in text and text.endswith(STRELKA_LIKE[6] + "\n")
