# Workflow-specific: variantcomp

> **Specific to one pipeline.** The profile `variantcomp` and the configs in
> `datasets/workflows/variantcomp/` write the inputs of
> [core-unit-bioinformatics/workflow-nxf-somatic-variantcomp](https://github.com/core-unit-bioinformatics/workflow-nxf-somatic-variantcomp)
> and of no other workflow. They follow that pipeline's samplesheet and params; when the
> pipeline changes, they change with it.

| | |
|---|---|
| Profile | `src/seqc2_toolkit/profiles/variantcomp.yaml` (`--profile variantcomp`) |
| Pipeline / revision | `core-unit-bioinformatics/workflow-nxf-somatic-variantcomp`, `dev` |
| Sources | SEQC2 only (`accepts`: technologies WGS, WES; caller muTect2) |
| Samplesheet | one row per VCF: `sample` (`SEQC2_<pair>`), technology, library type, tumor/normal ids, caller, VCF, BAMs, target BED |
| Params | reference, `truth_technology: wgs`, `callable_regions` / `read_evidence` / `igv_report` (on when every pair has BAMs), GENCODE GTF, GIAB stratifications, `known_variants.csv` (SEQC2 truth as benchmark + a chr22 demo list) |
| Remote datasets | no VEP annotation (the remote GTF is plain gzip), no known-variant list |

The pipeline compares somatic VCFs of the same tumor/normal pair sequenced with different
technologies (SEQC2 WGS vs. WES; WGS is the reference), so each pair needs both.

## Dataset configs (`datasets/workflows/variantcomp/`)

| Config | Nextflow profile | Content |
|---|---|---|
| `test.yaml` | `test` | remote: WGS + WES MuTect2 VCFs of IL_1, chr21 reference, local target BED |
| `test_full.yaml` | `test_full` | remote: the same, genome-wide (3.2 GB reference) |
| `offline.yaml` | `test_full_offline` | local, genome-wide, IL_1 (more with `--pairs` / `--centers`), GENCODE, GIAB |
| `small.yaml` | (CI / laptop) | local: IL_1, FD_1, NV_1, chr22, BAM slices of two 1-Mb windows, about 150 MB |
| `full.yaml` | (HPC) | local: first replicate of every center with WGS + WES, genome-wide, VCFs only |

```bash
seqc2 fetch datasets/workflows/variantcomp/small.yaml --outdir data/variantcomp_small
seqc2 run datasets/workflows/variantcomp/small.yaml --outdir data/variantcomp_small \
    --profile variantcomp --pipeline ../workflow-nxf-somatic-variantcomp --nf-profile conda
```
