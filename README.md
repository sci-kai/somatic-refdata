# somatic-refdata

Select, download and prepare public **somatic reference data** for testing and validating
bioinformatics workflows, and write the samplesheets and parameters that the workflows
need.

Two reference datasets are supported:

- **SEQC2** (`source: seqc2`, default): breast cancer cell line HCC1395 and its matched
  normal HCC1395BL; WGS and WES of several centers and replicates, somatic VCFs of 3
  aligners x 3 callers, BAMs, the SEQC2 somatic truth set.
- **GIAB HG008** (`source: hg008`): pancreatic cancer cell line HG008-T and normal
  duodenal / pancreatic tissue; Illumina, Element, PacBio Onso, PacBio HiFi, ONT and
  Ultima alignments, caller VCFs, the GIAB somatic draft benchmarks (small variants, SVs).

What the toolkit does:

- **Select** what a test needs: technologies, centers / series, replicates, aligners,
  callers, VCFs and/or alignments, contig subsets and slices of BAMs/CRAMs read remotely.
- **Download only what is missing**: one shared cache, md5-checked against the server's
  lists; derived files are rebuilt only when their recipe changes.
- **Or download nothing** (`access: remote`, SEQC2): the samplesheet points to the files
  on the server and the workflow reads them online (quick smoke tests).
- **Write workflow inputs** with small profiles: samplesheet, params and extra files for
  any (nf-core style) Nextflow workflow; optionally run it.

The core toolkit knows nothing about a specific workflow. Everything workflow-specific is in
a *profile* and its dataset configs, kept apart and marked as such (see
[below](#workflow-specific-profiles-and-configs)).

## Install

```bash
git clone https://github.com/sci-kai/somatic-refdata.git
cd somatic-refdata
conda env create -f envs/somatic-refdata.yaml   # python, pyyaml, samtools, bcftools, nextflow
conda activate somatic-refdata                   # installs this package (pip -e) as well
somref --help
```

Without conda: `pip install -e .` (Python >= 3.10; building datasets needs samtools,
bcftools, tabix and bgzip on `PATH`).

## Quick start

```bash
somref catalog                                  # what SEQC2 provides
somref catalog --source hg008                   # what GIAB HG008 provides
somref plan datasets/default.yaml               # what a config selects, sizes, cache state
somref fetch datasets/default.yaml --outdir data/default      # download + build
somref fetch datasets/hg008_small.yaml --outdir data/hg008_small
somref verify datasets/hg008_small.yaml --outdir data/hg008_small
somref samplesheet <config> --outdir <dataset> --profile <profile or FILE.yaml>
somref run <config> --outdir <dataset> --profile <profile> --nf-profile conda
```

The download cache is `$SOMREF_CACHE` or `~/.cache/somref` (`--cache`); genome-wide
datasets need several GB, BAMs 10-420 GB each, so use contig subsets and BAM windows for
small tests.

## Documentation

| Topic | |
|---|---|
| [docs/sources.md](docs/sources.md) | what SEQC2 and HG008 provide: technologies, pairs/series, VCFs, truth sets, what is not covered |
| [docs/datasets.md](docs/datasets.md) | dataset config reference (SEQC2, HG008), remote datasets, sample choice on the command line, sizes, cache and checksums, dataset folder layout |
| [docs/profiles.md](docs/profiles.md) | workflow profiles: format, template variables, params.yaml / params.config, `somref run` |
| [docs/testing.md](docs/testing.md) | three test tiers for workflows: `test`, `test_full`, `test_full_offline` |
| [CHANGELOG.md](CHANGELOG.md) | changes per version |

## Workflow-specific profiles and configs

> **These parts are specific to one pipeline each** and not needed for anything else.
> They are examples of the profile mechanism and inputs of the pipelines' test profiles.
> For your own workflow write a profile ([docs/profiles.md](docs/profiles.md)) and pass it
> with `--profile path/to/profile.yaml`.

| Profile | Pipeline (revision) | Sources | Profile file | Dataset configs | Docs |
|---|---|---|---|---|---|
| `variantcomp` | [core-unit-bioinformatics/workflow-nxf-somatic-variantcomp](https://github.com/core-unit-bioinformatics/workflow-nxf-somatic-variantcomp) (`dev`) | SEQC2 | [`profiles/variantcomp.yaml`](src/somatic_refdata/profiles/variantcomp.yaml) | [`datasets/workflows/variantcomp/`](datasets/workflows/variantcomp) | [docs/workflows/variantcomp.md](docs/workflows/variantcomp.md) |
| `variantinterpretation` | [cio-abcd/variantinterpretation](https://github.com/cio-abcd/variantinterpretation) (`v1.2.0`) | SEQC2, HG008 | [`profiles/variantinterpretation.yaml`](src/somatic_refdata/profiles/variantinterpretation.yaml) | [`datasets/workflows/variantinterpretation/`](datasets/workflows/variantinterpretation) | [docs/workflows/variantinterpretation.md](docs/workflows/variantinterpretation.md) |

`somref profiles` lists them with pipeline, homepage, configs and notes. Each profile file
starts with a `WORKFLOW-SPECIFIC PROFILE` header and names its pipeline in the `workflow`
section. General examples are directly in `datasets/`: `default.yaml` (SEQC2 defaults),
`hg008_default.yaml`, `hg008_small.yaml`.

## Repository layout

```
src/somatic_refdata/          the toolkit (command line: somref)
  sources_seqc2.yaml          SEQC2 server description
  sources_hg008.yaml          GIAB HG008 server description
  profiles/                   WORKFLOW-SPECIFIC profiles (one file per pipeline)
datasets/                     general example dataset configs
  workflows/<workflow>/       WORKFLOW-SPECIFIC dataset configs of the shipped profiles
docs/                         documentation; docs/workflows/: workflow-specific pages
envs/somatic-refdata.yaml     conda environment
tests/                        unit tests (no network)
```

## Development

```bash
pip install -e '.[dev]'
pytest                                  # unit tests, no network needed
black src tests && ruff check src tests
```

Contributions are welcome as issues and pull requests. Development follows the CUBI process:
`main` holds releases, `dev` development; work happens on `feat-*` branches from `dev`,
merged linearly (rebase, fast-forward). Versions follow [semver](https://semver.org); see
[CHANGELOG.md](CHANGELOG.md).

## License

MIT, see [LICENSE](LICENSE). The reference data are not part of this repository: they are
downloaded from their providers and are subject to the providers' terms.

## Citation

Cite this toolkit with the metadata in [CITATION.cff](CITATION.cff), and the data you use:

- SEQC2: Fang LT, Zhu B, Zhao Y, et al. Establishing community reference samples, data and
  call sets for benchmarking cancer mutation detection using whole-genome sequencing.
  *Nat Biotechnol* 2021;39:1151-1160 (doi:10.1038/s41587-021-00993-6). Xiao W, Ren L,
  Chen Z, et al. Toward best practice in cancer mutation detection with whole-genome and
  whole-exome sequencing. *Nat Biotechnol* 2021;39:1141-1150
  (doi:10.1038/s41587-021-00994-5).
- HG008: "Development and extensive sequencing of a broadly-consented Genome in a Bottle
  matched tumor-normal pair for somatic benchmarks", *Sci Data*
  (doi:10.1038/s41597-025-05438-2); the draft benchmarks: see their READMEs on the GIAB
  server.
- GIAB genome stratifications: Dwarshuis N, et al. *Nat Commun* 2024;15:9029.
- GENCODE: Frankish A, et al. *Nucleic Acids Res* 2023;51:D942.
