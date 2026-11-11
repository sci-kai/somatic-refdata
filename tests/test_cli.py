"""Command-line parsing."""

from seqc2_toolkit import cli


def test_nextflow_args_after_double_dash():
    argv = ["run", "ds.yaml", "--profile", "wf", "--dry-run", "--", "-with-report", "r.html"]
    seqc2_args, extra = cli.split_nextflow_args(argv)
    assert extra == ["-with-report", "r.html"]
    args = cli.parser().parse_args(seqc2_args)
    assert (args.config, args.profile, args.dry_run) == ("ds.yaml", "wf", True)


def test_human_sizes():
    assert [cli.human(s) for s in (None, 999, 1500, 3.2e9)] == ["?", "999 B", "1.5 KB", "3.2 GB"]


def test_with_defaults_completes_sections_of_older_datasets():
    from seqc2_toolkit import cli, dataset

    # a dataset.json of an older toolkit: data without reference_plain, no vcf section
    old_data = {k: v for k, v in dataset.DEFAULTS["data"].items() if k != "reference_plain"}
    built = {"data": old_data}
    assert cli.with_defaults(built, "data") == dataset.DEFAULTS["data"]
    assert cli.with_defaults(built, "vcf") == dataset.DEFAULTS["vcf"]
    # a value of the built dataset wins over the default
    built["data"] = {**old_data, "bam": not dataset.DEFAULTS["data"]["bam"]}
    assert cli.with_defaults(built, "data")["bam"] is not dataset.DEFAULTS["data"]["bam"]
