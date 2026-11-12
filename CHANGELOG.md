# somatic-refdata: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Dates follow
ISO 8601 (`YYYY-MM-DD`).

## Unreleased

### `Added`

- GIAB cancer reference HG008 (`source: hg008` in dataset configs, `somref catalog --source hg008`): 14 sequencing series (Illumina, Element, Onso, PacBio HiFi, ONT standard/ultra-long, Ultima CRAM) on GRCh38-GIABv3 and GRCh37 from the md5 list of each folder; tumor/normal pairs T-ND and T-NP (ONT normals from another series); NYGC MuTect2/Strelka2 and BCM DeepSomatic VCFs; somatic small-variant draft benchmark v0.4 (split into SNVs and indels, region tiers) and SV draft benchmark v0.6; GIAB reference with `.gzi`; remote CRAM slices decoded with the reference; example configs `hg008_default`, `hg008_small`, `workflows/variantinterpretation/hg008_small`
- Profiles: variables `source`, `source_name`, HG008 unit fields (`series`, `platform`, `build`, samples, coverage), `truth_sv*`, `reference_gzi`; per-VCF `tumor_id`/`normal_id`; `accepts.sources`
- Remote datasets (`access: remote`, SEQC2): nothing is downloaded; samplesheets and params hold the URLs of the original SEQC2 files, the workflow reads them online. With `subset.contigs` a smaller online reference is used (`remote_references`: GRCh38 chr21 of nf-core/test-datasets, checked identical to GRCh38.d1.vd1); BAM slices are refused
- `--technologies`, `--pairs`, `--centers`, `--replicates` for `plan`, `fetch`, `samplesheet`, `run` and `verify`: choose the samples without editing the config; `samplesheet` stops if the dataset folder was built with another selection
- Workflow profiles: `params.config` (the params as Nextflow config) next to `params.yaml`, with `--config-base` for files committed into a pipeline; `params_remote` / `files_remote` for remote datasets; template variables `access` and `contigs` (also recorded in dataset.json)
- `data.target_bed: local` for remote datasets: the WES target BED is written into the dataset folder (BED3, merged, gzip); VCFs and reference stay online
- `data.reference_plain`: contig subsets of the reference as plain FASTA (for tools that can't read bgzipped FASTA, e.g. GATK without .gzi)
- Dataset config section `vcf` (local datasets): `sample_order` (e.g. NORMAL before TUMOR) and `format_dp` (FORMAT/DP = sum of FORMAT/AD where missing), as needed by variantinterpretation's production filters; recorded in the recipe only when set, so existing datasets are not rebuilt
- Dataset configs for the test tiers of the two workflows: `test.yaml` and `test_full.yaml` (remote), `offline.yaml`
- Profiles: `variantinterpretation` sets `recompress_vcf` for remote datasets; `variantcomp` switches VEP and the known-variant lists off for remote datasets
- Profile `variantcomp`: `igv_report` (IGV views of variantcomp 0.16.0-dev+) is switched on like `callable_regions` and `read_evidence` when every pair has tumor and normal BAMs
- Profiles: optional `workflow` section (homepage, dataset configs, notes); `somref profiles` marks the shipped profiles as workflow-specific and prints pipeline, revision, homepage, configs and notes
- Publication on GitHub: `CITATION.cff`, project URLs and metadata in `pyproject.toml`, GitHub Actions workflow (unit tests, black, ruff on Python 3.10 and 3.12)
- Documentation in `docs/`: `sources.md` (SEQC2, HG008), `datasets.md` (config reference, remote datasets, cache, dataset folders), `profiles.md` (profile format, template variables), `testing.md` (three test tiers), `workflows/variantcomp.md`, `workflows/variantinterpretation.md` (workflow-specific)

### `Changed`

- Renamed from seqc2-toolkit to **somatic-refdata** (SEQC2 is no longer the only source): Python package `somatic_refdata`, command `somref` (was `seqc2`), conda env `somatic-refdata`, cache variable `SOMREF_CACHE` and default cache `~/.cache/somref`, dataset manifest `somref_manifest.json`, VCF header tag `##somatic_refdata_vcf`, SEQC2 server description `sources_seqc2.yaml`. Not backward compatible: dataset folders and caches of seqc2-toolkit are not recognized; build them again with `somref fetch`
- Workflow-specific dataset configs moved from `datasets/<workflow>_<tier>.yaml` to `datasets/workflows/<workflow>/<tier>.yaml` (dataset names unchanged: explicit `name`); shipped profiles start with a `WORKFLOW-SPECIFIC PROFILE` header
- README shortened for a public audience; details moved to `docs/`
- Profile `variantinterpretation`: sample names start with `{source_name}` (unchanged `SEQC2_...` for SEQC2 data)
- Dataset `README.md`, `dataset.json` and `params.yaml` name the source (SEQC2 / HG008); `dataset.json` records `source`

### `Fixed`

Found in the pre-release review; each has a regression test in `tests/test_fixes.py`.

- `somref samplesheet` / `run` crashed for every local HG008 dataset (`KeyError: 'vcf'`) and compared it with the SEQC2 defaults; now with the defaults of the dataset's source
- Rebuilding an output that was a hard link or symlink into the cache (e.g. a whole BAM, then `bam_windows`) wrote through the link and changed the cached file; old outputs are removed before a rebuild
- Downloads cut off early were accepted as complete (files without md5) or restarted from zero (with md5); the expected size (Content-Length / Content-Range) is checked and the download resumed; a partial file of the wrong size is not accepted on HTTP 416
- `somref run --pipeline <relative path>`: the path is made absolute (nextflow runs in the run folder)
- `subset.contigs` / `subset.bam_subsample` were silently ignored for whole BAMs (no `bam_windows`); now an error
- Remote datasets: `tumor_bai` / `normal_bai` / `bai` pointed to `.bam.bai` URLs that don't exist; now the server's index URLs
- `--pairs` / `--technologies` / `--centers` / `--replicates` on an HG008 config crashed; now a clear error
- An error page or throttling reply instead of an md5 list was cached as an (empty) list; md5 lists must hold md5 lines, also for the truth and GIAB lists; a GIAB file missing from its md5 list is an error (was used unchecked)
- Default cache folder `~/.cache/somref` (was still `~/.cache/seqc2`)
- Missing or broken files, bad values and YAML errors give a short error message instead of a traceback (traceback with `-v`); `somref plan` with an empty selection no longer crashes
- The dataset manifest is written atomically; `somref verify` on a folder without a dataset is an error (was "0 ok")
- Module-level config defaults were changed by command-line overrides (library use, tests)
- `params.config` header named every dataset SEQC2; remote datasets ignored the config's `cache`

## 0.1.0-dev

### `Added`

- `somref catalog`: catalog of the SEQC2 paired WGS/WES series from the server's md5 lists (VCFs of 3 aligners x 3 callers, BAMs, indices)
- Dataset configs (`select`, `data`, `subset`) with defaults (WGS of IL_1, bwa + MuTect2 VCF, reference, truth v1.2.1) and strict checks of keys, values and available combinations
- `somref plan`: selected files, sizes and cache state; `somref fetch`: md5-checked shared download cache, dataset folder with bgzipped/indexed VCFs, contig subsets, remote BAM slices, truth set, GENCODE GTF, GIAB stratifications; manifest of recipes so only changed outputs are rebuilt and deselected ones removed; `somref verify`
- Workflow profiles (samplesheet, params and extra files from templates; `accepts` / `requires` checks): `variantcomp`, `variantinterpretation`; `somref samplesheet`, `somref run`
- HTTPS downloads in Python with an explicit CA bundle (no curl), also exported for samtools/bcftools
- Example dataset configs: `default`, `variantcomp_small`, `variantcomp_full`, `variantinterpretation_small`
