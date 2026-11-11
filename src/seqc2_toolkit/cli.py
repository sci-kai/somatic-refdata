"""Command line: seqc2 <command> ...

  catalog      what SEQC2 provides (technologies, pairs, aligners, callers, BAMs);
               --source hg008: GIAB HG008 series, caller VCFs, benchmarks
  plan         what a dataset config selects: files, sizes, what is cached already
  fetch        download (only what is missing) and build the dataset folder
  samplesheet  write the samplesheet, params and extra files of a workflow profile
  run          samplesheet + nextflow run of a workflow profile
  verify       recompute the md5 of every cached source and dataset file
  profiles     list the shipped workflow profiles

Examples
  seqc2 catalog
  seqc2 catalog --source hg008
  seqc2 plan datasets/default.yaml
  seqc2 plan datasets/hg008_small.yaml
  seqc2 fetch datasets/variantcomp_small.yaml --outdir data/small
  seqc2 samplesheet datasets/variantcomp_small.yaml --profile variantcomp
  seqc2 run datasets/variantcomp_small.yaml --profile variantcomp \\
      --pipeline ../workflow-nxf-somatic-variantcomp --nf-profile conda
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

from . import __version__, build, catalog, dataset, net, profiles
from .cache import Cache, default_cache_dir, md5sum

log = logging.getLogger("seqc2")


def human(size):
    """File size for humans (None: unknown)."""
    if size is None:
        return "?"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1000 or unit == "TB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1000


def context(args, need_config=True):
    """(config or None, cache, catalog records) for a command."""
    config = None
    if need_config:
        config = dataset.load(args.config, {"outdir": args.outdir, "cache": args.cache})
    cache_dir = (config and config["cache"]) or args.cache or default_cache_dir()
    cache = Cache(cache_dir)
    source = config["source"] if config else getattr(args, "source", "seqc2")
    records = catalog.build(cache, refresh=getattr(args, "refresh", False), source=source)
    return config, cache, records


def outdir_of(config):
    """Dataset folder: --outdir / config outdir, else ./<dataset name>."""
    return Path(config["outdir"] or config["name"]).resolve()


def cmd_catalog(args):
    """Print the catalog summary (or every file with --files)."""
    _, _, records = context(args, need_config=False)
    if args.files:
        for r in records:
            if args.source == "hg008":
                # kind, series, sample (alignments) or caller (VCFs), build, coverage, url
                what = r.get("sample") or r.get("caller")
                cov = f"{r['coverage']}x" if r.get("coverage") else "-"
                fields = [r["kind"], r["series"], what, r["build"], cov, r["url"]]
            else:
                extra = r.get("caller") or r.get("role")
                fields = [r["kind"], r["technology"], r["pair"], r["aligner"], extra, r["url"]]
            print("\t".join(fields))
    else:
        print(catalog.summary(records, args.source))


def cmd_plan(args):
    """Print the selection and the remote files with size and cache status."""
    config, cache, records = context(args)
    units = dataset.select(config, records)
    print(f"Dataset {config['name']} -> {outdir_of(config)}\n")
    print(dataset.describe(units) + "\n")
    rows = build.plan(config, units, cache, catalog.builder_sources(config))
    width = max(len(r["item"]) for r in rows)
    todo = 0
    for r in rows:
        print(f"{r['item']:<{width}}  {human(r['size']):>9}  {r['how']:<26} cached: {r['cached']}")
        if r["cached"] == "missing" and r["size"]:
            todo += r["size"]
    print(f"\nto download: {human(todo)} (cache: {cache.root})")


def cmd_fetch(args):
    """Download what is missing and build the dataset folder."""
    config, cache, records = context(args)
    units = dataset.select(config, records)
    builder = build.Builder(
        config, catalog.builder_sources(config), cache, outdir_of(config), args.verify
    )
    builder.build(units)
    print(f"dataset ready: {builder.root}")


def write_profile(args, config):
    """Write the workflow files of --profile for the dataset; return (profile, workdir)."""
    profile = profiles.load_profile(args.profile)
    root = outdir_of(config)
    if not (root / "dataset.json").exists():
        raise dataset.ConfigError(f"{root}: no dataset.json; run 'seqc2 fetch {args.config}' first")
    written = profiles.write(profile, root, args.workdir)
    for path in written:
        print(f"wrote {path}")
    return profile, written[0].parent


def cmd_samplesheet(args):
    """Write samplesheet, params and extra files of a profile."""
    config = dataset.load(args.config, {"outdir": args.outdir, "cache": args.cache})
    write_profile(args, config)


def cmd_run(args):
    """Write the workflow files and run nextflow."""
    config = dataset.load(args.config, {"outdir": args.outdir, "cache": args.cache})
    profile, workdir = write_profile(args, config)
    rundir = Path(args.rundir or f"{config['name']}_{profile['name']}_run").resolve()
    params = profiles.resolve_params(workdir, rundir)
    extra = args.nextflow_args
    cmd = profiles.run_command(
        profile,
        params,
        rundir / "results",
        args.pipeline,
        args.revision,
        args.nf_profile,
        args.nextflow_config or [],
        ([] if args.no_resume else ["-resume"]) + extra,
    )
    print("$ " + " ".join(cmd) + f"   (in {rundir})")
    if args.dry_run:
        return
    sys.exit(subprocess.run(cmd, cwd=rundir).returncode)


def cmd_verify(args):
    """Recompute the md5 of every dataset file (manifest) and report mismatches."""
    config = dataset.load(args.config, {"outdir": args.outdir, "cache": args.cache})
    root = outdir_of(config)
    manifest = build.Manifest(root)
    bad = 0
    for rel, entry in sorted(manifest.entries.items()):
        path = root / rel
        ok = path.exists() and md5sum(path) == entry["md5"]
        bad += not ok
        print(f"{'ok' if ok else 'CHANGED/MISSING'}  {rel}")
    print(f"{len(manifest.entries) - bad} ok, {bad} changed or missing (run fetch to rebuild)")
    sys.exit(1 if bad else 0)


def cmd_profiles(_args):
    """List the shipped workflow profiles."""
    for name in profiles.shipped_profiles():
        profile = profiles.load_profile(name)
        print(f"{name}: {profile['description']}")


def parser():
    """The argument parser with all subcommands."""
    top = argparse.ArgumentParser(
        prog="seqc2",
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(__doc__.split("\n")[1:]),
    )
    top.add_argument("--version", action="version", version=f"seqc2-toolkit {__version__}")
    top.add_argument("-v", "--verbose", action="store_true", help="debug messages")
    top.add_argument("--ca-bundle", help="CA certificate bundle for HTTPS")
    top.add_argument("--cache", help="download cache (default: $SEQC2_CACHE or ~/.cache/seqc2)")
    sub = top.add_subparsers(dest="command", required=True)

    def with_config(p):
        p.add_argument("config", help="dataset config (YAML)")
        p.add_argument("--outdir", help="dataset folder (default: config outdir or ./<name>)")

    p = sub.add_parser("catalog", help="what SEQC2 (or HG008) provides")
    p.add_argument(
        "--source",
        choices=dataset.SOURCE_NAMES,
        default="seqc2",
        help="reference data: seqc2 (HCC1395, default) or hg008 (GIAB HG008)",
    )
    p.add_argument("--refresh", action="store_true", help="fetch the server's md5 lists again")
    p.add_argument("--files", action="store_true", help="list every file")
    p.set_defaults(func=cmd_catalog)
    p = sub.add_parser("plan", help="what a dataset config selects and what is cached")
    with_config(p)
    p.set_defaults(func=cmd_plan)
    p = sub.add_parser("fetch", help="download what is missing and build the dataset")
    with_config(p)
    p.add_argument("--verify", action="store_true", help="hash every cached and dataset file again")
    p.set_defaults(func=cmd_fetch)
    for name, func, text in (
        ("samplesheet", cmd_samplesheet, "write the workflow files of a profile"),
        ("run", cmd_run, "write the workflow files and run nextflow"),
    ):
        p = sub.add_parser(name, help=text)
        with_config(p)
        p.add_argument(
            "--profile", required=True, help="workflow profile: shipped name or YAML path"
        )
        p.add_argument(
            "--workdir", help="folder of the workflow files (default: <dataset>/<profile>)"
        )
        p.set_defaults(func=func)
        if name == "run":
            p.add_argument(
                "--pipeline", help="nextflow pipeline (repo or local path) instead of the profile's"
            )
            p.add_argument("--revision", help="pipeline revision (-r)")
            p.add_argument("--nf-profile", help="nextflow -profile (e.g. conda, singularity)")
            p.add_argument(
                "-c", "--nextflow-config", action="append", help="extra nextflow config file"
            )
            p.add_argument("--rundir", help="run folder (default: ./<dataset>_<profile>_run)")
            p.add_argument("--no-resume", action="store_true")
            p.add_argument("--dry-run", action="store_true", help="only print the nextflow command")
    p = sub.add_parser("verify", help="recompute the md5 of every dataset file")
    with_config(p)
    p.set_defaults(func=cmd_verify)
    p = sub.add_parser("profiles", help="list the shipped workflow profiles")
    p.set_defaults(func=cmd_profiles)
    return top


def split_nextflow_args(argv):
    """(seqc2 arguments, nextflow arguments): everything after '--' goes to nextflow."""
    if "--" not in argv:
        return argv, []
    cut = argv.index("--")
    return argv[:cut], argv[cut + 1 :]


def main(argv=None):
    """Entry point of the seqc2 command."""
    argv, extra = split_nextflow_args(sys.argv[1:] if argv is None else list(argv))
    args = parser().parse_args(argv)
    args.nextflow_args = extra
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    net.setup_tls(args.ca_bundle)
    try:
        args.func(args)
    except (dataset.ConfigError, RuntimeError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
