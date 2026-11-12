# Workflow profiles

The toolkit itself knows nothing about a specific workflow. A **profile** is a small YAML
file that says how one workflow takes a dataset: which samplesheet (one row per VCF, per
pair x technology or per alignment), which params, which extra files. `seqc2 samplesheet`
writes them into `<dataset>/<profile name>/` (or `--workdir`); `seqc2 run` also runs
Nextflow.

Two profiles are shipped, each written for one pipeline (`seqc2 profiles` lists them;
[README](../README.md#workflow-specific-profiles-and-configs)). For any other workflow
write your own profile and pass the file: `--profile path/to/myworkflow.yaml`. Copying a
shipped profile is a good start.

## Format

```yaml
name: myworkflow                           # default: the file name
description: what the workflow does        # shown by 'seqc2 profiles'
workflow:                                  # optional: the workflow the profile is for
  url: https://github.com/org/myworkflow   #   homepage
  configs: datasets/workflows/myworkflow/  #   its dataset configs
  notes: required pipeline revision, params, ...
pipeline: org/myworkflow                   # nextflow run target (or --pipeline)
revision: v1.0.0                           # -r (or --revision; left out for local paths)
nf_profile: conda                          # -profile (or --nf-profile)
nextflow_args: []                          # further nextflow arguments
accepts:                                   # the dataset must fit; each key optional
  sources: [seqc2, hg008]
  technologies: [WGS, WES]                 # SEQC2 technologies or HG008 series
  callers: [muTect2]
  aligners: [bwa]
requires: [vcf, reference]                 # dataset content needed (vcf, bam, resources)
paths: absolute                            # or relative (to the samplesheet folder)
samplesheet:
  file: samplesheet.csv
  rows: vcf                                # vcf: per VCF; pair: per pair x technology;
                                           # bam: per alignment
  missing: n/a                             # value of empty fields
  separator: ","
  columns:                                 # column: template
    sample: "{source_name}_{pair}_{tech}_{caller}"
    vcf: "{vcf}"
params:                                    # params.yaml / params.config;
  input: "{samplesheet}"                   #   a param with an empty variable is left out
  fasta: "{reference}"
params_remote:                             # replaces params for remote datasets;
  recompress_vcf: true                     #   null removes one
files:                                     # extra files: relative path -> template text
  known.csv: |
    name,variants
    truth,{truth_snv}
files_remote: {}                           # the same for remote datasets (null: not written)
```

Unknown keys and values stop with an error. `accepts` and `requires` are checked against
the dataset before anything is written.

## Template variables

Templates are Python format strings (`{variable}`).

| Variables | Content |
|---|---|
| per row | `pair`, `center`, `replicate`, `technology` (WGS, or the HG008 series), `tech` (lower case), `library_type` (wgs / capture), `aligner`, `caller`, `vcf`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `tumor_id`, `normal_id` (sample columns of the VCF), `target_bed`; rows `bam`: `role`, `bam`, `bai` |
| per dataset (also in rows) | `dataset`, `source` (seqc2 / hg008), `source_name` (SEQC2 / HG008), `dataset_dir`, `workdir`, `samplesheet`, `pairs`, `first_pair`, `all_bams` (true/false: every pair has tumor and normal alignments), `reference`, `reference_fai`, `target_bed_wes`, `truth_snv`, `truth_indel`, `truth_regions`, `gencode_gtf`, `giab_<name>`, `access` (local / remote), `contigs` (comma-separated: the records are, or for remote datasets should be, restricted to them) |
| HG008 only | per row: `series`, `platform`, `build`, `tumor_sample`, `normal_sample`, `tumor_coverage`, `normal_coverage` (`tumor_bai` etc. are `.crai` for CRAMs); per dataset: `reference_gzi`, `truth_sv`, `truth_sv_regions`, `truth_sv_regions_clonal` |

An empty variable gives `missing` (`n/a`) in a samplesheet column, leaves a param out, and
is an error in an extra file. Rendered params `true` / `false` become booleans.

## params.yaml, params.config and seqc2 run

Next to `params.yaml`, `params.config` holds the same params as a Nextflow config, for `-c`
or for `includeConfig` in a test profile of the workflow. Local files of the workflow folder
are written as `"<config base>/<path>"`; the base is the absolute folder path, or the
Nextflow expression given with `--config-base`, e.g. `'${projectDir}/assets/seqc2/test'` for
files committed into a pipeline, or `"${params.pipelines_testdata_base_path}myworkflow"` for
a dataset folder that is moved.

`seqc2 run` writes the workflow files, a params file with absolute paths into the run
folder and runs

```
nextflow run <pipeline> -r <revision> -profile <nf_profile> -params-file ... \
    --outdir <run>/results -resume [nextflow arguments after --]
```

`--dry-run` only prints the command.
