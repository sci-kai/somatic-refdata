"""Workflow profiles: templates, samplesheets, params, checks."""

import json

import pytest
import yaml

from seqc2_toolkit import profiles
from seqc2_toolkit.dataset import ConfigError


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
    assert len(written) == 3


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
