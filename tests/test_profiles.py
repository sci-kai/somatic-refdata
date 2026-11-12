"""Workflow profiles: templates, samplesheets, params, checks."""

import json

import pytest
import yaml

from somatic_refdata import profiles
from somatic_refdata.dataset import ConfigError


def make_dataset(tmp_path, bams=True):
    """A dataset folder with dataset.json (two technologies of IL_1)."""
    root = tmp_path / "ds"
    root.mkdir()
    units = []
    for tech in ("WGS", "WES"):
        units.append(
            {
                "pair": "IL_1",
                "center": "IL",
                "replicate": 1,
                "technology": tech,
                "target_bed": "bed/wes.bed" if tech == "WES" else None,
                "vcfs": [
                    {
                        "aligner": "bwa",
                        "caller": "muTect2",
                        "path": f"vcf/IL_1_{tech.lower()}.vcf.gz",
                    }
                ],
                "bams": (
                    {"T": f"bam/IL_1_{tech}_T.bam", "N": f"bam/IL_1_{tech}_N.bam"} if bams else {}
                ),
            }
        )
    facts = {
        "name": "ds",
        "units": units,
        "resources": {"reference": "reference/ref.fa.gz", "truth_snv": "truth/snv.vcf.gz"},
        "vcf_samples": {"tumor": "TUMOR", "normal": "NORMAL"},
    }
    (root / "dataset.json").write_text(json.dumps(facts))
    return root


def profile(tmp_path, text):
    path = tmp_path / "wf.yaml"
    path.write_text(text)
    return profiles.load_profile(path)


def test_render_missing_variables():
    assert profiles.render("{a}-{b}", {"a": 1, "b": None}) == ("1-", ["b"])


def test_samplesheet_params_files(tmp_path):
    root = make_dataset(tmp_path)
    prof = profile(
        tmp_path,
        "name: wf\npaths: relative\n"
        "samplesheet:\n  columns:\n    sample: 'S_{pair}'\n    tech: '{tech}'\n"
        "    bed: '{target_bed}'\n    vcf: '{vcf}'\n"
        "params:\n  input: '{samplesheet}'\n  fasta: '{reference}'\n"
        "  gtf: '{gencode_gtf}'\n  bams: '{all_bams}'\n"
        "files:\n  extra.txt: 'truth={truth_snv}'\n",
    )
    written = profiles.write(prof, root)
    sheet = (root / "wf" / "samplesheet.csv").read_text().splitlines()
    assert sheet == [
        "sample,tech,bed,vcf",
        "S_IL_1,wgs,n/a,../vcf/IL_1_wgs.vcf.gz",
        "S_IL_1,wes,../bed/wes.bed,../vcf/IL_1_wes.vcf.gz",
    ]
    params = yaml.safe_load((root / "wf" / "params.yaml").read_text())
    # gtf is left out (not in the dataset); 'true' is a boolean
    assert params == {"input": "samplesheet.csv", "fasta": "../reference/ref.fa.gz", "bams": True}
    assert (root / "wf" / "extra.txt").read_text() == "truth=../truth/snv.vcf.gz"
    # params.config is written last: a param may name an extra file
    assert [p.name for p in written] == [
        "samplesheet.csv",
        "params.yaml",
        "extra.txt",
        "params.config",
    ]
    config = (root / "wf" / "params.config").read_text()
    # local files as absolute paths, booleans as Groovy booleans; the reference file
    # doesn't exist in this test, so its relative value is written as it is
    assert f"input = '{root / 'wf' / 'samplesheet.csv'}'" in config
    assert "fasta = '../reference/ref.fa.gz'" in config
    assert "bams  = true" in config


def test_resolve_params_absolute(tmp_path):
    root = make_dataset(tmp_path)
    (root / "reference").mkdir()
    (root / "reference" / "ref.fa.gz").write_text(">x\n")
    prof = profile(
        tmp_path, "name: wf\npaths: relative\nparams:\n  fasta: '{reference}'\n  x: '7'\n"
    )
    profiles.write(prof, root)
    params = yaml.safe_load(profiles.resolve_params(root / "wf", tmp_path / "run").read_text())
    assert params == {"fasta": str((root / "reference" / "ref.fa.gz").resolve()), "x": "7"}


def test_accepts_and_requires(tmp_path):
    root = make_dataset(tmp_path, bams=False)
    with pytest.raises(ConfigError, match="strelka"):
        profiles.write(profile(tmp_path, "accepts: {callers: [strelka]}\n"), root)
    with pytest.raises(ConfigError, match="requires gencode_gtf"):
        profiles.write(profile(tmp_path, "requires: [gencode_gtf]\n"), root)
    with pytest.raises(ConfigError, match="no BAMs"):
        profiles.write(profile(tmp_path, "samplesheet: {rows: bam}\n"), root)


def test_shipped_profiles_load():
    for name in profiles.shipped_profiles():
        prof = profiles.load_profile(name)
        assert prof["pipeline"] and prof["samplesheet"]["columns"]


def test_run_command(tmp_path):
    prof = profiles.load_profile("variantinterpretation")
    cmd = profiles.run_command(prof, tmp_path / "params.yaml", tmp_path / "out", extra=["-resume"])
    assert cmd[:6] == [
        "nextflow",
        "run",
        "cio-abcd/variantinterpretation",
        "-r",
        "v1.2.0",
        "-profile",
    ]
    # a local pipeline folder has no revision
    cmd = profiles.run_command(prof, tmp_path / "p.yaml", tmp_path / "o", pipeline=str(tmp_path))
    assert "-r" not in cmd


def make_remote_dataset(tmp_path):
    """A remote dataset folder: URLs instead of paths."""
    root = tmp_path / "remote"
    root.mkdir()
    facts = {
        "name": "remote",
        "access": "remote",
        "contigs": ["chr21"],
        "units": [
            {
                "pair": "IL_1",
                "center": "IL",
                "replicate": 1,
                "technology": "WES",
                "target_bed": None,
                "vcfs": [
                    {
                        "aligner": "bwa",
                        "caller": "muTect2",
                        "path": "https://example.org/WES_IL_1.vcf.gz",
                        "md5": "x",
                    }
                ],
                "bams": {},
            }
        ],
        "resources": {"reference": "https://example.org/chr21.fa"},
        "vcf_samples": {"tumor": "TUMOR", "normal": "NORMAL"},
    }
    (root / "dataset.json").write_text(json.dumps(facts))
    return root


def test_remote_urls_params_remote_and_config_base(tmp_path):
    root = make_remote_dataset(tmp_path)
    prof = profile(
        tmp_path,
        "name: wf\npaths: relative\n"
        "samplesheet:\n  columns:\n    sample: 'S_{pair}'\n    vcf: '{vcf}'\n"
        "params:\n  input: '{samplesheet}'\n  fasta: '{reference}'\n  vep: true\n"
        "  regions: '{contigs}'\n  bams: '{all_bams}'\n"
        "params_remote:\n  vep: false\n  bams: null\n"
        "files:\n  a.txt: 'a'\n  b.txt: 'b'\nfiles_remote:\n  b.txt: null\n",
    )
    workdir = tmp_path / "pipeline" / "assets" / "test"
    profiles.write(prof, root, workdir, config_base="${projectDir}/assets/test")
    # URLs are written as they are, also with paths: relative
    assert (workdir / "samplesheet.csv").read_text().splitlines()[1] == (
        "S_IL_1,https://example.org/WES_IL_1.vcf.gz"
    )
    params = yaml.safe_load((workdir / "params.yaml").read_text())
    # with config_base, files of the workflow folder are relative to it
    assert params == {
        "input": "samplesheet.csv",
        "fasta": "https://example.org/chr21.fa",
        "vep": False,
        "regions": "chr21",
    }
    config = (workdir / "params.config").read_text().splitlines()
    assert '    input   = "${projectDir}/assets/test/samplesheet.csv"' in config
    assert "    fasta   = 'https://example.org/chr21.fa'" in config
    assert "    vep     = false" in config
    # files_remote: b.txt is not written for a remote dataset
    assert (workdir / "a.txt").exists() and not (workdir / "b.txt").exists()
    # no absolute path of this machine in a file meant to be committed
    assert str(tmp_path) not in "\n".join(config)


def test_nextflow_params_quoting(tmp_path):
    (tmp_path / "a b$c.csv").write_text("")
    text = profiles.nextflow_params(
        {"x": "it's", "n": 3, "f": "a b$c.csv", "u": "s3://b/k"}, tmp_path, "${base}"
    )
    assert "x = 'it\\'s'" in text
    assert "n = 3" in text
    # only the base expression is interpolated
    assert 'f = "${base}/a b\\$c.csv"' in text
    assert "u = 's3://b/k'" in text


def test_shipped_profiles_name_their_workflow():
    """Every shipped profile is workflow-specific: it names its pipeline and configs."""
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    for name in profiles.shipped_profiles():
        profile = profiles.load_profile(name)
        assert profile["pipeline"] and profile["workflow"]["url"].startswith("https://")
        assert (root / profile["workflow"]["configs"]).is_dir()


def test_profile_workflow_keys_are_checked(tmp_path):
    path = tmp_path / "mine.yaml"
    path.write_text("workflow: {homepage: https://example.org}\n")
    with pytest.raises(ConfigError, match="unknown key.*homepage"):
        profiles.load_profile(path)
    path.write_text("samplesheet: {columns: {vcf: '{vcf}'}}\n")
    assert profiles.load_profile(path)["workflow"]["url"] is None
