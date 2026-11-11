# seqc2-toolkit: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Dates follow
ISO 8601 (`YYYY-MM-DD`).

## Unreleased (0.1.0-dev)

### `Added`

- `seqc2 catalog`: catalog of the SEQC2 paired WGS/WES series from the server's md5 lists (VCFs of 3 aligners x 3 callers, BAMs, indices)
- Dataset configs (`select`, `data`, `subset`) with defaults (WGS of IL_1, bwa + MuTect2 VCF, reference, truth v1.2.1) and strict checks of keys, values and available combinations
- `seqc2 plan`: selected files, sizes and cache state; `seqc2 fetch`: md5-checked shared download cache, dataset folder with bgzipped/indexed VCFs, contig subsets, remote BAM slices, truth set, GENCODE GTF, GIAB stratifications; manifest of recipes so only changed outputs are rebuilt and deselected ones removed; `seqc2 verify`
- Workflow profiles (samplesheet, params and extra files from templates; `accepts` / `requires` checks): `variantcomp`, `variantinterpretation`; `seqc2 samplesheet`, `seqc2 run`
- HTTPS downloads in Python with an explicit CA bundle (no curl), also exported for samtools/bcftools
- Example dataset configs: `default`, `variantcomp_small`, `variantcomp_full`, `variantinterpretation_small`
- GIAB cancer reference HG008 (`source: hg008` in dataset configs, `seqc2 catalog --source hg008`): 14 sequencing series (Illumina, Element, Onso, PacBio HiFi, ONT standard/ultra-long, Ultima CRAM) on GRCh38-GIABv3 and GRCh37 from the md5 list of each folder; tumor/normal pairs T-ND and T-NP (ONT normals from another series); NYGC MuTect2/Strelka2 and BCM DeepSomatic VCFs; somatic small-variant draft benchmark v0.4 (split into SNVs and indels, region tiers) and SV draft benchmark v0.6; GIAB reference with `.gzi`; remote CRAM slices decoded with the reference; example configs `hg008_default`, `hg008_small`, `hg008_variantinterpretation_small`
- Profiles: variables `source`, `source_name`, HG008 unit fields (`series`, `platform`, `build`, samples, coverage), `truth_sv*`, `reference_gzi`; per-VCF `tumor_id`/`normal_id`; `accepts.sources`

### `Changed`

- Profile `variantinterpretation`: sample names start with `{source_name}` (unchanged `SEQC2_...` for SEQC2 data)
- Dataset `README.md`, `dataset.json` and `params.yaml` name the source (SEQC2 / HG008); `dataset.json` records `source`
