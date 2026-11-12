# Testing a workflow: three tiers

The toolkit prepares the inputs of three kinds of workflow tests. The workflow repository
holds the Nextflow test profiles; the toolkit writes their samplesheets and params (with
the workflow's [profile](profiles.md)).

| Nextflow profile | Data | Internet | Purpose | Dataset config |
|---|---|---|---|---|
| `test` | a fraction of one sample: one pair, records of one contig (chr21) | yes | does the pipeline run through (setup check) | `access: remote`, `subset.contigs: [chr21]` |
| `test_full` | every record of one sample (IL_1) | yes | the same on a full sample | `access: remote` |
| `test_full_offline` | one to all SEQC2 samples (`--pairs` / `--centers`) or HG008 series (`select.series`), chosen by the user | no (after `somref fetch`) | validation on an offline HPC | `access: local` |

Remote datasets (`test`, `test_full`) are SEQC2 only; HG008 datasets are always built
locally (`somref fetch`).

```bash
# test / test_full: remote datasets; the files are committed into the pipeline (assets/)
somref samplesheet datasets/workflows/variantinterpretation/test.yaml \
    --profile variantinterpretation \
    --workdir ../variantinterpretation/assets/seqc2/test \
    --config-base '${projectDir}/assets/seqc2/test'

# test_full_offline: download, prepare, write the workflow files into the dataset folder
somref fetch datasets/workflows/variantinterpretation/offline.yaml --outdir data/vi \
    --pairs IL_1,FD_1,NV_1
somref samplesheet datasets/workflows/variantinterpretation/offline.yaml --outdir data/vi \
    --pairs IL_1,FD_1,NV_1 --profile variantinterpretation
nextflow run <pipeline> -profile test_full_offline,singularity \
    --pipelines_testdata_base_path data/vi/
```

For a new workflow: write its profile, put one config per tier into
`datasets/workflows/<workflow>/` (`test.yaml`, `test_full.yaml`, `offline.yaml`), and add
the Nextflow test profiles to the workflow that `includeConfig` the written
`params.config`.
