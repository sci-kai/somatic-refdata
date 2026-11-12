# Workflow-specific: variantinterpretation

> **Specific to one pipeline.** The profile `variantinterpretation` and the configs in
> `datasets/workflows/variantinterpretation/` write the inputs of
> [cio-abcd/variantinterpretation](https://github.com/cio-abcd/variantinterpretation)
> and of no other workflow. They follow that pipeline's samplesheet and params; when the
> pipeline changes, they change with it.

| | |
|---|---|
| Profile | `src/somatic_refdata/profiles/variantinterpretation.yaml` (`--profile variantinterpretation`) |
| Pipeline / revision | `cio-abcd/variantinterpretation`, `v1.2.0` |
| Sources | SEQC2 and HG008 (`accepts`: caller muTect2; allele fractions from FORMAT/AD) |
| Samplesheet | one row per VCF: `sample` (`<SEQC2\|HG008>_<pair>_<tech>_<aligner>_<caller>`), `vcf`, `group` (`<SEQC2\|HG008>_<pair>`) |
| Params | `input`, `fasta`; remote datasets: `recompress_vcf: true` |
| Not from the toolkit | the Ensembl VEP cache: `--vep_cache` (and `--vep_cache_version`, 115 for v1.2.0) on the command line or in a config (`-c`) |

**Revision note:** remote datasets set `recompress_vcf` (the SEQC2 VCFs are plain gzip).
The pipeline has this param only since its development branch `feat-seqc2-tests`, not in
`v1.2.0`: run remote datasets with a pipeline revision that has it (`--revision`,
`--pipeline`) until it is released. Local datasets work with `v1.2.0`.

## Dataset configs (`datasets/workflows/variantinterpretation/`)

| Config | Nextflow profile | Content |
|---|---|---|
| `test.yaml` | `test` | remote: MuTect2 VCF of WES IL_1, chr21 reference |
| `test_full.yaml` | `test_full` | remote: the same, genome-wide (needs a full VEP cache) |
| `offline.yaml` | `test_full_offline` | local, WES IL_1 genome-wide (more with `--pairs` / `--centers`), truth set |
| `small.yaml` | (CI / laptop) | local: WES IL_1, chr22 |
| `hg008_small.yaml` | (CI / laptop) | local HG008: NYGC MuTect2 VCF (T vs N-D), chr22, GRCh38-GIABv3, GENCODE 46 |

```bash
somref fetch datasets/workflows/variantinterpretation/hg008_small.yaml --outdir data/vi_hg008
somref run datasets/workflows/variantinterpretation/hg008_small.yaml --outdir data/vi_hg008 \
    --profile variantinterpretation --nf-profile singularity -- --vep_cache /path/to/vep_cache
```
