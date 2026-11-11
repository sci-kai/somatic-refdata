"""GIAB HG008: catalog from md5 lists, config defaults and checks, selection, sources."""

import json

import pytest

from seqc2_toolkit import catalog, dataset, hg008, profiles
from seqc2_toolkit.dataset import ConfigError

# md5 lists of the server (folder name -> lines), in the formats seen on the server
MD5_LISTS = {
    "NYGC_Illumina-WGS_20231023/checksums.md5": [
        "a1  HG008-T_Illumina_161x_GRCh38-GIABv3.bam",
        "a2  HG008-T_Illumina_161x_GRCh38-GIABv3.bam.bai",
        "a3  HG008-N-D_Illumina_118x_GRCh38-GIABv3.bam",
        "a4  HG008-N-D_Illumina_118x_GRCh38-GIABv3.bam.bai",
        "a5  HG008-T_Illumina_161x_GRCh37.bam",
        "a6  HG008-T_Illumina_161x_GRCh37.bam.bai",
        "a7  HG008-N-D_Illumina_118x_GRCh37.bam",
        "a8  HG008-N-D_Illumina_118x_GRCh37.bam.bai",
        "a9  HG008-T_Illumina_161x_CHM13v2.0.bam",  # CHM13: not in the catalog
        "aa  HG008-T_Illumina_161x_CHM13v2.0.bam.bai",
        "ab  HG008-T_fastqs/HVWWFDSX7_L002_R1_001.fastq.gz",
    ],
    "BCM_Revio_20240313/checksum.md5": [
        "b1  ./HG008-T_PacBio-HiFi-Revio_20240313_106x_GRCh38-GIABv3.bam",
        "b2  ./HG008-T_PacBio-HiFi-Revio_20240313_106x_GRCh38-GIABv3.bam.bai",
        "b3  ./HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam",
        "b4  ./HG008-N-D_PacBio-HiFi-Revio_20240313_68x_GRCh38-GIABv3.bam.bai",
        "b5  ./HG008-T_PacBio-Revio_m84039_240113_032943_s4.hifi_reads.bc2005.bam",  # unaligned
    ],
    "UCSC_ONT-UL_20231207/checksum.md5": [
        "c1  HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado0.4.3_sup4.2.0_5mCG_5hmCG_54x_UCSC_"
        "20231031.bam",
        "c2  HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado0.4.3_sup4.2.0_5mCG_5hmCG_54x_UCSC_"
        "20231031.bam.bai",
        "c3  HG008-T_UCSC_20231031_ONT-UL-R10.4.1_1_dorado-v0.4.3_sup4.2.0-5mCG-5hmCG.bam",
    ],
    "Northeastern_ONT-std_20240422/checksum.md5": [
        "d1  HG008-N-D_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_94x.bam",
        "d2  HG008-N-D_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_94x.bam.bai",
        "d3  HG008-N-D_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_13x_gt25kb.bam",
        "d4  HG008-N-D_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_13x_gt25kb.bam.bai",
        "d5  HG008-N-P_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_41x.bam",
        "d6  HG008-N-P_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_41x.bam.bai",
    ],
    "Ultima_UG100-WGS_20240530/checksum.md5": [
        "e1\tHG008-T_Ultima_306x_GRCh38-GIABv3.cram",
        "e2\tHG008-T_Ultima_306x_GRCh38-GIABv3.cram.crai",
        "e3\tHG008-N-D_Ultima_234x_GRCh38-GIABv3.cram",  # no index: skipped
    ],
    "Element_AVITI_20240626/checksum.md5": [
        "f1  HG008-T_Element-StdInsert_102x_GRCh38-GIABv3.bam",
        "f2  HG008-T_Element-StdInsert_102x_GRCh38-GIABv3.bam.bai",
        "f3  HG008-T_Element-LngInsert_31x_GRCh38-GIABv3.bam",
        "f4  HG008-T_Element-LngInsert_31x_GRCh38-GIABv3.bam.bai",
        "f5  HG008-N-P_Element-StdInsert_70x_GRCh38-GIABv3.bam",
        "f6  HG008-N-P_Element-StdInsert_70x_GRCh38-GIABv3.bam.bai",
    ],
    "analysis/NYGC-somatic-pipeline_20240412/GRCh38-GIABv3/checksums.md5": [
        "g1  HG008-T--HG008-N.mutect2.vcf.gz",
        "g2  HG008-T--HG008-N.snv.strelka2.vcf.gz",
        "g3  HG008-T--HG008-N.indel.strelka2.vcf.gz",
    ],
    "analysis/BCM_Revio_20240313/pacbio-wgs-wdl_somatic_20240601/checksum.md5": [
        "h1  /GRCh38-GIABv3/HG008_deepsomatic.PASS.hiphase.vep.vcf.gz",
    ],
}


@pytest.fixture
def records(monkeypatch, hex_md5):
    """The HG008 catalog of the md5 lists above (other folders: empty lists)."""

    def fake_md5_list(cache, url, refresh=False):
        for suffix, lines in MD5_LISTS.items():
            if url.endswith("/Liss_lab/" + suffix):
                return lines
        return []

    monkeypatch.setattr(catalog, "md5_list", fake_md5_list)
    return catalog.build(None, source="hg008")


def config(tmp_path, text):
    """A dataset config from YAML text."""
    path = tmp_path / "ds.yaml"
    path.write_text(text)
    return dataset.load(path)


def find(records, **fields):
    """The records with these field values."""
    return [r for r in records if all(r.get(k) == v for k, v in fields.items())]


def test_parse_md5_formats():
    # md5s must be 32 hex digits (e.g. not an HTML error page)
    assert hg008.parse_md5(["aaaa  ./x.bam", "<title>Object not found!</title>", ""]) == {}
    lines = [f"{'a' * 32}  ./x.bam", f"{'b' * 32}\t/sub/y.vcf.gz", f"{'c' * 32} z.bai"]
    assert hg008.parse_md5(lines) == {
        "x.bam": "a" * 32,
        "sub/y.vcf.gz": "b" * 32,
        "z.bai": "c" * 32,
    }


@pytest.fixture
def hex_md5(monkeypatch):
    """The fixture md5 lists use short ids: accept any word as md5 in the catalog tests."""

    def lenient(lines):
        result = {}
        for line in lines:
            parts = line.split()
            if len(parts) == 2:
                result[hg008.strip_path(parts[1])] = parts[0]
        return result

    monkeypatch.setattr(hg008, "parse_md5", lenient)


def test_catalog_alignments(records):
    alns = find(records, kind="aln")
    # CHM13, unaligned reads, FASTQs, gt25kb subsets and unindexed files are left out
    assert all(r["build"] in ("GRCh38", "GRCh37") for r in alns)
    assert not [r for r in alns if "gt25kb" in r["name"] or "hifi_reads" in r["name"]]
    assert not find(alns, series="ultima", sample="N-D")
    t = find(alns, series="illumina-nygc", sample="T", build="GRCh38")[0]
    assert (t["coverage"], t["format"], t["md5"]) == (161, "bam", "a1")
    assert t["url"].endswith("/Liss_lab/NYGC_Illumina-WGS_20231023/" + t["name"])
    assert t["index"] == {"url": t["url"] + ".bai", "md5": "a2"}
    hifi = find(alns, series="hifi-bcm", sample="T")[0]
    assert hifi["md5"] == "b1" and "/./" not in hifi["url"]
    ont = find(alns, series="ont-ul-ucsc", sample="T")[0]
    assert (ont["build"], ont["coverage"]) == ("GRCh38", 54)  # GRCh38_GIABv3 token
    cram = find(alns, series="ultima", sample="T")[0]
    assert cram["format"] == "cram" and cram["index"]["url"].endswith(".cram.crai")
    # two libraries in one folder: two series
    std = find(alns, series="element-2406")
    lng = find(alns, series="element-long-2406")
    assert {r["coverage"] for r in std} == {102, 70} and [r["coverage"] for r in lng] == [31]


def test_catalog_vcfs(records):
    vcfs = find(records, kind="vcf")
    assert {r["caller"] for r in vcfs} == {
        "muTect2",
        "strelka2-snv",
        "strelka2-indel",
        "deepsomatic",
    }
    deep = find(vcfs, caller="deepsomatic")[0]
    assert deep["md5"] == "h1"
    assert deep["url"].endswith(
        "pacbio-wgs-wdl_somatic_20240601/GRCh38-GIABv3/HG008_deepsomatic.PASS.hiphase.vep.vcf.gz"
    )
    assert deep["samples"] == {"tumor": "HG008.tumor", "normal": None}
    mutect = find(vcfs, caller="muTect2")[0]
    assert (mutect["series"], mutect["normal"], mutect["build"]) == (
        "illumina-nygc",
        "N-D",
        "GRCh38",
    )


def test_catalog_two_alignments_of_one_sample_is_an_error(monkeypatch, hex_md5):
    lines = [
        "x1  HG008-T_Illumina_161x_GRCh38-GIABv3.bam",
        "x2  HG008-T_Illumina_161x_GRCh38-GIABv3.bam.bai",
        "x3  HG008-T_Illumina_80x_GRCh38-GIABv3.bam",
        "x4  HG008-T_Illumina_80x_GRCh38-GIABv3.bam.bai",
    ]
    monkeypatch.setattr(
        catalog,
        "md5_list",
        lambda cache, url, refresh=False: lines if "NYGC_Illumina" in url else [],
    )
    with pytest.raises(RuntimeError, match="two alignments of T GRCh38"):
        catalog.build(None, source="hg008")


def test_load_defaults_and_unknown_source(tmp_path):
    cfg = config(tmp_path, "source: hg008\n")
    assert cfg["select"] == hg008.DEFAULTS["select"]
    assert cfg["data"]["truth"] == "v0.4" and cfg["data"]["truth_sv"] is False
    assert config(tmp_path, "{}\n")["source"] == "seqc2"
    with pytest.raises(ConfigError, match="source giab unknown"):
        config(tmp_path, "source: giab\n")
    # SEQC2 keys are not HG008 keys
    with pytest.raises(ConfigError, match="unknown key.*centers"):
        config(tmp_path, "source: hg008\nselect: {centers: [IL]}\n")


def test_select_default(tmp_path, records):
    units = dataset.select(config(tmp_path, "source: hg008\n"), records)
    assert len(units) == 1
    unit = units[0]
    assert (unit["pair"], unit["technology"], unit["build"]) == ("T-ND", "illumina-nygc", "GRCh38")
    assert unit["samples"] == {"T": "HG008-T", "N": "HG008-N-D"}
    assert unit["coverage"] == {"T": 161, "N": 118}
    assert [v["caller"] for v in unit["vcfs"]] == ["muTect2"]
    assert unit["vcfs"][0]["aligner"] == "bwa-mem"
    assert unit["bams"] == {}


def test_select_ont_normal_from_other_series(tmp_path, records):
    cfg = config(
        tmp_path,
        "source: hg008\nselect: {series: [ont-ul-ucsc], normals: [N-D, N-P], callers: []}\n"
        "data: {bam: [T-NP]}\n",
    )
    units = dataset.select(cfg, records)
    assert [u["pair"] for u in units] == ["T-ND", "T-NP"]
    assert units[0]["bams"] == {} and units[0]["coverage"] == {"T": 54, "N": 94}
    np_unit = units[1]
    assert np_unit["bams"]["N"]["series"] == "ont-ne"
    assert np_unit["bam_names"] == {"T": "HG008-T_ont-ul-ucsc", "N": "HG008-N-P_ont-ne"}
    assert np_unit["bais"]["T"]["url"].endswith(".bam.bai")


def test_select_errors(tmp_path, records):
    def select(text):
        return dataset.select(config(tmp_path, "source: hg008\n" + text), records)

    # NYGC sequenced no pancreatic normal
    with pytest.raises(ConfigError, match="T-NP illumina-nygc: no N-P alignment"):
        select("select: {normals: [N-P]}\n")
    # MuTect2 VCFs exist only for the NYGC series
    with pytest.raises(ConfigError, match="caller muTect2: VCFs only for series illumina-nygc"):
        select("select: {series: [hifi-bcm]}\n")
    select("select: {series: [hifi-bcm], callers: [deepsomatic]}\n")
    select("select: {series: [hifi-bcm]}\ndata: {vcf: false}\n")
    # normals-only series can't be selected
    with pytest.raises(ConfigError, match="select.series: ont-ne not available"):
        select("select: {series: [ont-ne], callers: []}\n")
    with pytest.raises(ConfigError, match="truth_sv: v0.6 is not available for GRCh37"):
        select("select: {build: GRCh37}\ndata: {truth_sv: v0.6}\n")
    with pytest.raises(ConfigError, match="GENCODE GTFs are used for GRCh38 only"):
        select("select: {build: GRCh37}\ndata: {gencode: 46}\n")
    with pytest.raises(ConfigError, match="truth_tier: strict unknown"):
        select("data: {truth_tier: strict}\n")
    with pytest.raises(ConfigError, match="select.build: CHM13 not available"):
        select("select: {build: CHM13}\n")
    with pytest.raises(ConfigError, match="data.bam: T-NP not among the selected pairs"):
        select("data: {bam: [T-NP]}\n")


def test_builder_sources_build_and_tier(tmp_path):
    cfg = config(
        tmp_path,
        "source: hg008\nselect: {build: GRCh37}\ndata: {truth_tier: nogermlineinterference}\n",
    )
    src = catalog.builder_sources(cfg)
    assert src["name"] == "HG008"
    assert src["reference"]["fasta"].endswith("GRCh37/hs37d5.fa.gz")
    assert src["reference"]["gzi"] == src["reference"]["fasta"] + ".gzi"
    truth = src["truth"]["v0.4"]
    assert truth["snv"]["types"] == "snps" and truth["indel"]["types"] == "indels"
    assert truth["snv"]["path"] == truth["indel"]["path"]
    assert truth["snv"]["md5_key"] == (
        "GRCh37/HG008-T_somatic_smvar_benchmark_v0.4_GRCh37_tumorvariants.vcf.gz"
    )
    assert truth["regions"]["path"].endswith("_GRCh37_nogermlineinterference.bed")
    assert src["truth_sv"] == {}  # SV benchmark: GRCh38 only
    giab = src["annotation"]["giab"]
    assert giab["base_url"].endswith("/v3.5/GRCh37@all")
    assert giab["files"]["all_difficult"] == "Union/GRCh37_alldifficultregions.bed.gz"
    # SEQC2 configs keep the unchanged server description
    assert catalog.builder_sources(config(tmp_path, "{}\n")) == catalog.load_sources()


def test_profile_rows_and_checks(tmp_path):
    root = tmp_path / "ds"
    root.mkdir()
    facts = {
        "name": "hg008_test",
        "source": "hg008",
        "source_name": "HG008",
        "resources": {"reference": "reference/ref.fa.gz", "truth_sv": "truth_sv/sv.vcf.gz"},
        "vcf_samples": {"tumor": None, "normal": None},
        "units": [
            {
                "pair": "T-ND",
                "center": "Ultima",
                "replicate": None,
                "technology": "ultima",
                "series": "ultima",
                "platform": "Ultima UG100",
                "aligner": "ua",
                "build": "GRCh38",
                "library_type": "wgs",
                "samples": {"T": "HG008-T", "N": "HG008-N-D"},
                "coverage": {"T": 306, "N": 234},
                "target_bed": None,
                "vcfs": [
                    {
                        "aligner": "ua",
                        "caller": "deepsomatic",
                        "path": "vcf/t.vcf.gz",
                        "tumor_id": "HG008.tumor",
                        "normal_id": None,
                    }
                ],
                "bams": {"T": "bam/HG008-T_ultima.cram", "N": "bam/HG008-N-D_ultima.bam"},
            }
        ],
    }
    (root / "dataset.json").write_text(json.dumps(facts))
    profile = {
        **profiles.PROFILE_DEFAULTS,
        "name": "demo",
        "accepts": {"sources": ["hg008"]},
        "samplesheet": {
            **profiles.PROFILE_DEFAULTS["samplesheet"],
            "rows": "pair",
            "columns": {
                "sample": "{source_name}_{pair}_{series}",
                "tumor": "{tumor_id}",
                "normal": "{normal_id}",
                "tumor_bai": "{tumor_bai}",
                "normal_bai": "{normal_bai}",
                "cov": "{tumor_coverage}/{normal_coverage}",
                "type": "{library_type}",
            },
        },
        "params": {"input": "{samplesheet}", "sv_truth": "{truth_sv}"},
        "paths": "relative",
    }
    written = profiles.write(profile, root)
    sheet = written[0].read_text().splitlines()
    assert sheet[1].split(",") == [
        "HG008_T-ND_ultima",
        "HG008.tumor",
        "n/a",
        "../bam/HG008-T_ultima.cram.crai",
        "../bam/HG008-N-D_ultima.bam.bai",
        "306/234",
        "wgs",
    ]
    assert "for the HG008 dataset hg008_test" in written[1].read_text()
    with pytest.raises(ConfigError, match="sources: hg008 not supported"):
        profiles.write({**profile, "accepts": {"sources": ["seqc2"]}}, root)
