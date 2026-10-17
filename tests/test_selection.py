"""Catalog parsing, dataset config and selection."""

import pytest

from seqc2_toolkit import catalog, dataset

BASE = "https://example.org/seqc2"


def records():
    """A small catalog: IL_1/IL_2 WGS+WES, NS_1 WGS only; bwa/novo x muTect2/strelka; BAMs."""
    sources = catalog.load_sources()
    vcf_lines, bam_lines, bai_lines = [], [], []
    for tech, pairs in (("WGS", ["IL_1", "IL_2", "NS_1"]), ("WES", ["IL_1", "IL_2"])):
        for pair in pairs:
            for aligner in ("bwa", "novo"):
                for caller in ("muTect2", "strelka"):
                    vcf_lines.append(f"aaa  {tech}/{tech}_{pair}.{aligner}.{caller}.vcf.gz")
            center, rep = pair.split("_")
            for role in "TN":
                bam_lines.append(f"bbb   {tech}_{center}_{role}_{rep}.bwa.dedup.bam")
                bai_lines.append(f"ccc {tech}_{center}_{role}_{rep}.bwa.dedup.bai")
    # a file of another series: not in the catalog
    vcf_lines.append("ddd  FFG/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz")
    result = []
    for kind, lines in (("vcf", vcf_lines), ("bam", bam_lines), ("bai", bai_lines)):
        spec = sources["md5_lists"][kind]
        result += catalog.parse_md5_list(
            lines, kind, sources["patterns"][kind], spec["folder"], BASE
        )
    return result


def config(tmp_path, text):
    """A dataset config from YAML text."""
    path = tmp_path / "ds.yaml"
    path.write_text(text)
    return dataset.load(path)


def test_parse_fields_and_urls():
    recs = records()
    vcf = next(r for r in recs if r["kind"] == "vcf")
    assert vcf["url"] == f"{BASE}/analysis/SNVs/vcfs/WGS/WGS_IL_1.bwa.muTect2.vcf.gz"
    assert (vcf["pair"], vcf["technology"], vcf["aligner"], vcf["caller"]) == (
        "IL_1",
        "WGS",
        "bwa",
        "muTect2",
    )
    bam = next(r for r in recs if r["kind"] == "bam" and r["technology"] == "WES")
    # BAMs are listed by name only; their folder is data/<technology>
    assert bam["url"] == f"{BASE}/data/WES/WES_IL_T_1.bwa.dedup.bam"
    assert not any("FFG" in r["url"] for r in recs)


def test_default_selection(tmp_path):
    units = dataset.select(config(tmp_path, "{}"), records())
    assert [(u["pair"], u["technology"]) for u in units] == [("IL_1", "WGS")]
    assert [(v["aligner"], v["caller"]) for v in units[0]["vcfs"]] == [("bwa", "muTect2")]
    assert units[0]["bams"] == {}


def test_all_centers_replicates_and_bams(tmp_path):
    cfg = config(
        tmp_path,
        "select: {technologies: [WGS, WES], centers: all, replicates: all,"
        " callers: [muTect2, strelka]}\n"
        "data: {bam: [IL_2]}\n",
    )
    with pytest.raises(dataset.ConfigError, match="NS_1 WES"):
        dataset.select(cfg, records())
    cfg["select"]["centers"] = ["IL"]
    units = dataset.select(cfg, records())
    assert [(u["pair"], u["technology"]) for u in units] == [
        ("IL_1", "WGS"),
        ("IL_1", "WES"),
        ("IL_2", "WGS"),
        ("IL_2", "WES"),
    ]
    assert len(units[0]["vcfs"]) == 2
    assert units[0]["bams"] == {} and set(units[2]["bams"]) == {"T", "N"}


def test_explicit_pairs_replace_centers_and_keep_order(tmp_path):
    units = dataset.select(
        config(tmp_path, "select: {pairs: [NS_1, IL_2, IL_1], centers: [IL]}\n"), records()
    )
    assert [u["pair"] for u in units] == ["NS_1", "IL_2", "IL_1"]


def test_unknown_key_and_value(tmp_path):
    with pytest.raises(dataset.ConfigError, match="unknown key"):
        config(tmp_path, "select: {caller: [muTect2]}\n")
    with pytest.raises(dataset.ConfigError, match="choose from: muTect2, strelka"):
        dataset.select(config(tmp_path, "select: {callers: [mutect]}\n"), records())
    with pytest.raises(dataset.ConfigError, match="not among the selected pairs"):
        dataset.select(config(tmp_path, "data: {bam: [IL_2]}\n"), records())


def test_relative_outdir(tmp_path):
    cfg = config(tmp_path, "outdir: data/small\n")
    assert cfg["outdir"] == str((tmp_path / "data" / "small").resolve())
    assert cfg["name"] == "ds"


def test_command_line_paths_relative_to_working_directory(tmp_path, monkeypatch):
    (tmp_path / "configs").mkdir()
    path = tmp_path / "configs" / "ds.yaml"
    path.write_text("outdir: from_config\n")
    monkeypatch.chdir(tmp_path)
    assert dataset.load(path)["outdir"] == str(tmp_path / "configs" / "from_config")
    assert dataset.load(path, {"outdir": "cli"})["outdir"] == str(tmp_path / "cli")
