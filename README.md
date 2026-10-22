# seqc2-toolkit

Select, download and prepare the SEQC2 somatic reference data (breast cancer cell line
HCC1395 and its matched normal HCC1395BL) for testing bioinformatics workflows, and write
the samplesheets and parameters that the workflows need. The same works for the GIAB cancer
reference HG008 (pancreatic cancer cell line HG008-T, normal tissues HG008-N-D/-N-P):
Illumina, Element, Onso, PacBio HiFi, ONT and Ultima alignments, caller VCFs and the GIAB
somatic draft benchmarks (`source: hg008`).

- **Select** what a workflow needs from the SEQC2 data: technologies, sequencing
  centers and replicates, aligners, callers, VCFs and/or BAMs, contig subsets and BAM
  slices.
- **Download only what is missing**: one shared cache, md5-checked against the server's
  lists; derived files are rebuilt only when their recipe changes.
- **Or download nothing** (`access: remote`): the samplesheet points to the files on the
  SEQC2 server and the workflow reads them online (quick smoke tests with internet).
- **Write workflow inputs** with small profiles: samplesheet, params and extra files for
  any (nf-core style) workflow; optionally run it with Nextflow.

The toolkit knows nothing about a specific workflow: everything workflow-specific is in a
profile (`src/seqc2_toolkit/profiles/`, or your own YAML file).

## Install

```bash
conda env create -f envs/seqc2-toolkit.yaml   # python, pyyaml, samtools, bcftools, nextflow
conda activate seqc2-toolkit                   # installs this package (pip -e) as well
seqc2 --help
```

Without conda: `pip install -e .` (needs samtools, bcftools, tabix and bgzip on `PATH`
for building datasets).

## Quick start

```bash
seqc2 catalog                                     # what SEQC2 provides
seqc2 plan datasets/variantcomp_small.yaml        # what a config selects, sizes, cache
seqc2 fetch datasets/variantcomp_small.yaml --outdir data/variantcomp_small
seqc2 samplesheet datasets/variantcomp_small.yaml --outdir data/variantcomp_small \
    --profile variantcomp                         # -> data/variantcomp_small/variantcomp/
seqc2 run datasets/variantcomp_small.yaml --outdir data/variantcomp_small \
    --profile variantcomp --pipeline ../workflow-nxf-somatic-variantcomp --nf-profile conda
```

## Testing a workflow: three tiers

The toolkit prepares the inputs of three kinds of workflow tests. The workflow repository
holds the Nextflow test profiles; the toolkit writes their samplesheets and params.

| Profile | Data | Internet | Purpose | Dataset config (examples) |
|---|---|---|---|---|
| `test` | a fraction of one sample: one pair, records of one contig (chr21) | yes | does the pipeline run through (setup check) | `*_test.yaml` (`access: remote`) |
| `test_full` | every record of one sample (IL_1) | yes | the same on a full sample | `*_test_full.yaml` (`access: remote`) |
| `test_full_offline` | one to all SEQC2 samples, chosen by the user | no (after `seqc2 fetch`) | validation on an offline HPC | `*_offline.yaml` + `--pairs` / `--centers` |

```bash
# test / test_full: remote datasets; the files are committed into the pipeline (assets/)
seqc2 samplesheet datasets/variantinterpretation_test.yaml --profile variantinterpretation \
    --workdir ../variantinterpretation/assets/seqc2/test --config-base '${projectDir}/assets/seqc2/test'
# test_full_offline: download, prepare, write the workflow files into the dataset folder
seqc2 fetch datasets/variantinterpretation_offline.yaml --outdir data/vi --pairs IL_1,FD_1,NV_1
seqc2 samplesheet datasets/variantinterpretation_offline.yaml --outdir data/vi --pairs IL_1,FD_1,NV_1 \
    --profile variantinterpretation
nextflow run <pipeline> -profile test_full_offline,singularity --pipelines_testdata_base_path data/vi/
```

## What SEQC2 provides

The SEQC2 Somatic Mutation Working Group sequenced HCC1395 / HCC1395BL at several
centers, with several technologies, and called somatic variants with several pipelines
(Fang et al. 2021, Xiao et al. 2021). The toolkit covers the paired tumor/normal series:

| | WGS | WES |
|---|---|---|
| Pairs (center_replicate) | 21: EA_1, FD_1-3, IL_1-3, LL_1, NC_1, NS_1-9, NV_1-3 | 12: EA_1, FD_1-3, IL_1-3, LL_1, NC_1, NV_1-3 |
| VCFs | 3 aligners (bwa, bowtie, novo) x 3 callers (muTect2, strelka, somaticSniper) | same |
| BAMs | bwa, tumor + normal, every pair (10-160 GB each) | same |

Plus: the GRCh38.d1.vd1 reference, the exome target BED, the SEQC2 somatic truth set
(v1.2, v1.2.1) with its high-confidence regions; optional GENCODE and GIAB
stratification files. Not yet in the catalog: the FFPE (FFG, FFX), library-prep (LBP),
purity/depth (SPP) and AmpliSeq series, and FASTQs (SRA only).

`seqc2 catalog` shows the current state of the server.

## What GIAB HG008 provides (`source: hg008`)

HG008 is the Genome in a Bottle (GIAB) cancer reference: the pancreatic ductal
adenocarcinoma cell line HG008-T (batch 0823p23) and two normal tissues of the same donor,
HG008-N-D (duodenal) and HG008-N-P (pancreatic). A dataset config with `source: hg008`
selects from it; `seqc2 catalog --source hg008` shows the current state of the server.

| Series | Platform (center, aligner) | Samples |
|---|---|---|
| `illumina-bcm`, `illumina-nygc` | Illumina PCR-free 2x150 (BCM, NYGC; bwa-mem) | T, N-D |
| `element-2401`, `element-2406`, `element-2412` | Element AVITI standard insert (bwa-mem) | T, N-D (+ N-P: 2406, 2412) |
| `element-long-2406` | Element AVITI long insert (bwa-mem) | T, N-D, N-P |
| `onso` | PacBio Onso short reads (bwa-mem) | T, N-D |
| `hifi-pacbio`, `hifi-bcm` | PacBio HiFi Revio (PacBio, BCM; pbmm2) | T + N-P (PacBio), T + N-D (BCM) |
| `ont-ucsc`, `ont-ul-ucsc`, `ont-ul-ne` | ONT R10.4.1 standard / ultra-long (UCSC, Northeastern; minimap2) | T; normals from `ont-ne` |
| `ont-ne` | ONT R10.4.1 standard, normals only (Northeastern; minimap2) | N-D, N-P |
| `ultima` | Ultima UG100 (ua), CRAM | T, N-D |

- **Alignments** (BAM, Ultima: CRAM) on the GIAB mapping references GRCh38-GIABv3 and
  GRCh37 (hs37d5), md5-checked against each folder's list. 30-300x; one BAM is 50-420 GB,
  so read slices with `subset.bam_windows` (CRAMs are decoded with the remote reference
  and written as BAM).
- **Caller VCFs**: `muTect2`, `strelka2-snv`, `strelka2-indel` (NYGC pipeline,
  `illumina-nygc`, T vs N-D, GRCh38 + GRCh37) and `deepsomatic` (BCM, `hifi-bcm`, tumor
  column only, GRCh38). The sample columns differ per caller; profiles get them per VCF as
  `tumor_id` / `normal_id`.
- **Truth**: GIAB somatic small-variant draft benchmark v0.4 (GRCh38, GRCh37; one VCF,
  split by the toolkit into SNVs and indels; regions tier `all`, `nogermlineinterference`
  or `nogermlinewithin50bp`) and the somatic SV draft benchmark v0.6 (GRCh38; PASS VCF,
  regions `clonal_and_subclonal` and `clonal`). Both are drafts: read their READMEs on the
  server before drawing conclusions.
- Not in the catalog: CHM13 alignments (two different CHM13 FASTAs, not on the GIAB
  server), series on other references (Ultima ppmSeq, BCM DRAGEN 2025, Hartwig/Roche
  Axelios, Dovetail), Hi-C, Bionano, single cell, other or superseded passages, NIST
  clones, FASTQs, further callers (ClairS, Severus, GRIDSS2, DRAGEN, ...).

```yaml
source: hg008
select:
  series: [illumina-nygc]   # series (table above), or all
  normals: [N-D]            # N-D and/or N-P: pairs T-ND, T-NP
  build: GRCh38             # GRCh38 (GIABv3) or GRCh37 (hs37d5)
  callers: [muTect2]        # for the selected series that have them; [] for none
data:
  vcf: true
  bam: false                # false, true or a list of pairs (T-ND, T-NP)
  reference: true           # GIAB reference of the build (+ .fai, .gzi)
  truth: v0.4               # small-variant draft benchmark, or false
  truth_tier: all           # all, nogermlineinterference, nogermlinewithin50bp
  truth_sv: false           # v0.6 (GRCh38 only), or false
  gencode: false            # GRCh38 only
  giab: []                  # GIAB stratifications of the build
subset: {}                  # as for SEQC2 (no capture targets)
```

**Default (`source: hg008` only):** T vs N-D, NYGC Illumina, MuTect2 VCF, GRCh38 reference
and the small-variant benchmark, genome-wide. A combination that HG008 doesn't provide (e.g.
`illumina-nygc` with N-P, or `muTect2` without `illumina-nygc`) stops with an error.
Examples: `datasets/hg008_default.yaml`, `datasets/hg008_small.yaml` (chr22, 4 platforms,
BAM slices, both benchmarks; about 60 MB), `datasets/hg008_variantinterpretation_small.yaml`.

HG008 dataset folders use the same layout as SEQC2 with these differences:
`reference/GRCh38_GIABv3[_<contigs>].fa.gz` (+ `.fai`, `.gzi`),
`vcf/<pair>_<series>.<aligner>.<caller>.vcf.gz`, alignments named by sample
`bam/HG008-<T|N-D|N-P>_<series>.bam|cram` (a tumor alignment is shared by T-ND and T-NP),
`truth/..._tumorvariants_snps|indels.vcf.gz`, `truth_sv/`.

## Dataset config

A YAML file with three optional sections and a few top-level keys; everything left out
takes the default.

```yaml
access: local             # local: download and build (seqc2 fetch); remote: see below
select:
  technologies: [WGS]     # WGS and/or WES
  centers: [IL]           # sequencing centers, or all
  replicates: [1]         # replicate numbers, or all
  pairs: null             # explicit pairs, e.g. [IL_1, FD_2] (replaces centers/replicates)
  aligners: [bwa]         # bwa, bowtie, novo
  callers: [muTect2]      # muTect2, strelka, somaticSniper
data:
  vcf: true               # the somatic VCFs of the selection
  bam: false              # false, true (every selected pair) or a list of pairs
  reference: true         # GRCh38.d1.vd1 FASTA + .fai (3.2 GB, or only subset.contigs)
  reference_plain: false   # contig subsets as plain FASTA (default: bgzip + .gzi)
  target_bed: true        # exome target BED (when WES is selected)
  truth: v1.2.1           # SEQC2 truth set release, or false
  gencode: false          # GENCODE basic GTF release (e.g. 46), or false
  giab: []                # GIAB stratifications: homopolymer, tandem_repeat,
                          #   segmental_duplication, low_mappability, all_difficult
subset:
  contigs: null           # only these contigs in every file, e.g. [chr22]
  bam_windows: null       # BAMs: only these windows (read remotely, no full download)
  bam_targets: true       # BAMs: within the windows only the WES targets + padding
  bam_padding: 100
  bam_subsample: null     # fraction of read pairs (fixed seed 42)
vcf:                      # optional changes of the VCFs (local datasets only)
  sample_order: null      # sample columns in this order, e.g. [NORMAL, TUMOR]
  format_dp: false        # add FORMAT/DP = sum(FORMAT/AD) where a record has none
name: null                # default: the config file name
outdir: null              # default: ./<name>; relative to the config file
cache: null               # default: $SEQC2_CACHE or ~/.cache/seqc2
```

**Default (empty config):** WGS of one center and replicate, IL_1, aligned with bwa and
called with MuTect2; VCF, reference and truth set, genome-wide.

Unknown keys, unknown values (e.g. a misspelled caller) and combinations that SEQC2
doesn't provide (e.g. NS pairs with WES) stop with an error that names the allowed values.

`vcf` is for workflows that expect a sample order (e.g. filters that index the samples by
position) or a read depth per sample: the SEQC2 MuTect2 (GATK 3.8) records have no
FORMAT/DP. Records with FORMAT/DP are kept; the header documents the change
(`##seqc2_toolkit_vcf=...`).

### Choosing the samples on the command line

`plan`, `fetch`, `samplesheet`, `run` and `verify` take `--technologies`, `--pairs`,
`--centers` and `--replicates` (comma-separated, or `all`); they replace the values of the
config's `select` section. `--centers` / `--replicates` clear the config's `pairs`. Use the
same options for `fetch` and `samplesheet`: `samplesheet` stops if the dataset folder was
built with another selection.

```bash
seqc2 fetch datasets/variantcomp_offline.yaml --pairs IL_1,FD_1,NV_1   # three pairs
seqc2 fetch datasets/variantcomp_offline.yaml --centers all            # first replicate of every center
seqc2 fetch datasets/variantinterpretation_offline.yaml --centers all --replicates all --technologies WES,WGS
```

### Remote datasets (`access: remote`)

Nothing is downloaded (except the server's md5 lists for the catalog). `seqc2 fetch` (or
`seqc2 samplesheet` directly) writes only `dataset.json` and `README.md`; every file in the
samplesheet and params is a URL that the workflow reads itself (Nextflow stages it):

- VCFs: the original SEQC2 files (plain gzip, not bgzip, no index; sample columns TUMOR,
  NORMAL; all contigs). The workflow must accept plain gzip.
- Reference: without `subset.contigs` the whole GRCh38.d1.vd1 FASTA (3.2 GB). With
  `subset.contigs` a smaller online reference with exactly these contigs
  (`remote_references` in `sources.yaml`; currently `chr21` from nf-core/test-datasets,
  sequence identical to GRCh38.d1.vd1). The VCFs still hold every contig: the workflow has
  to drop the other records (template variable `{contigs}`, e.g. `bcftools view --targets`).
- Target BED, truth set, GENCODE GTF, GIAB files: the original files (URLs).
- Not possible: BAM slices (`bam_windows`, `bam_subsample`); whole BAMs are possible
  but 10-160 GB each.

### How much data does a workflow need?

| Question | Recommendation |
|---|---|
| Replicates | One per center (`replicates: [1]`). Replicates of one center are the same DNA, library kit and sequencer; they only matter for reproducibility studies. Several centers show real technical variability. |
| Technologies | WGS for genome-wide tests; WES for capture workflows (with `target_bed`); both for cross-technology comparisons. |
| Aligners / callers | One (bwa + muTect2) unless the workflow compares callers or aligners. MuTect2 writes FORMAT/AD and AF; Strelka only tier counts (AU/CU/GU/TU), SomaticSniper DP4/BCOUNT: check which fields your workflow reads. |
| BAMs | Only for workflows that read alignments (coverage, read evidence, calling). With `subset.bam_windows` the toolkit reads small slices remotely instead of downloading 10-160 GB per BAM. |
| Size | For CI / laptops: `contigs: [chr22]`, BAM windows of 1-2 Mb, `bam_subsample: 0.5`: about 150 MB for 3 pairs x 2 technologies. |

Examples in `datasets/`: `default.yaml`; for workflow-nxf-somatic-variantcomp
`variantcomp_test.yaml`, `variantcomp_test_full.yaml` (remote), `variantcomp_offline.yaml`,
`variantcomp_small.yaml` (chr22 + BAM slices), `variantcomp_full.yaml`; for
variantinterpretation `variantinterpretation_test.yaml`, `variantinterpretation_test_full.yaml`
(remote), `variantinterpretation_offline.yaml`, `variantinterpretation_small.yaml`.

## Downloads, cache and checksums

- Every remote file is downloaded once into the cache
  (`<cache>/<host>/<path>`, resumable) and checked against the md5 of the server's lists
  (VCFs, BAMs, BAM indices, GIAB files). `<cache>/registry.json` records size, time and
  md5 of every cached file, so it is not hashed again unless it changed (`--verify`
  hashes everything again).
- The dataset folder holds the derived files (bgzipped and indexed VCFs, subsets, BAM
  slices, ...). `seqc2_manifest.json` records the recipe (source, md5, settings) and md5
  of each; `seqc2 fetch` rebuilds only files whose recipe changed or that were modified,
  and removes files that the config no longer selects.
- `checksums.md5`, `README.md` (sources, settings) and `dataset.json` (the dataset for
  profiles) are written into the dataset folder. `seqc2 verify` checks every file.
- HTTPS: downloads are done in Python; the CA bundle is taken from `--ca-bundle`,
  `$CURL_CA_BUNDLE`, `$SSL_CERT_FILE`, the Python package certifi or the system paths, and
  exported for samtools/bcftools reading remote files. This avoids
  `curl: (77) error adding trust anchors` with conda builds of curl.

## Dataset folder

```
reference/GRCh38.d1.vd1[_<contigs>].fa[.gz]  (+ .fai, .gzi)
bed/S07604624_Covered_human_all_v6_plus_UTR_hg38.bed
vcf/<pair>_<technology>.<aligner>.<caller>.vcf.gz  (+ .tbi; sample columns TUMOR, NORMAL)
bam/<pair>_<technology>_<T|N>.bam  (+ .bam.bai)
truth/high-confidence_sSNV|sINDEL_in_HC_regions_<release>.vcf.gz, High-Confidence_Regions_*.bed
annotation/gencode.v<release>.basic.annotation.gtf.gz (+ .tbi)
annotation/giab/GRCh38/<GIAB release path>
dataset.json, checksums.md5, README.md, seqc2_manifest.json
<profile>/samplesheet.csv, params.yaml, ...   (seqc2 samplesheet / run)
```

## Workflow profiles

A profile describes how one workflow takes the data. Shipped: `variantcomp`
(core-unit-bioinformatics/workflow-nxf-somatic-variantcomp) and `variantinterpretation`
(cio-abcd/variantinterpretation); `seqc2 profiles` lists them. Your own profile is a
YAML file passed with `--profile path/to/profile.yaml`:

```yaml
name: variantinterpretation
pipeline: cio-abcd/variantinterpretation   # nextflow run target (or --pipeline)
revision: v1.2.0                           # -r (or --revision)
nf_profile: conda                          # -profile (or --nf-profile)
accepts: {callers: [muTect2]}              # the dataset must fit (technologies, callers, aligners)
requires: [vcf, reference]                 # dataset content needed
paths: absolute                            # or relative (to the samplesheet folder)
samplesheet:
  file: samplesheet.csv
  rows: vcf                                # vcf: per VCF; pair: per pair x technology; bam: per BAM
  columns:
    sample: "SEQC2_{pair}_{tech}_{aligner}_{caller}"
    vcf: "{vcf}"
params:                                    # params.yaml; a param with an empty variable is left out
  input: "{samplesheet}"
  fasta: "{reference}"
params_remote:                             # replaces params for remote datasets; null removes one
  recompress_vcf: true
files: {}                                  # extra files: relative path -> template text
files_remote: {}                           # the same for extra files (null: not written)
```

Templates are Python format strings with these variables:

| Variables | Content |
|---|---|
| per row | `pair`, `center`, `replicate`, `technology` (WGS), `tech` (wgs), `library_type` (wgs / capture), `aligner`, `caller`, `vcf`, `tumor_bam`, `tumor_bai`, `normal_bam`, `normal_bai`, `tumor_id` (TUMOR), `normal_id` (NORMAL), `target_bed`; rows `bam`: `role`, `bam`, `bai` |
| per dataset | `dataset`, `source` (seqc2 / hg008), `source_name` (SEQC2 / HG008), `dataset_dir`, `workdir`, `samplesheet`, `pairs`, `first_pair`, `all_bams` (true/false), `reference`, `reference_fai`, `target_bed_wes`, `truth_snv`, `truth_indel`, `truth_regions`, `gencode_gtf`, `giab_<name>`, `access` (local/remote), `contigs` (comma-separated: the records are, or for remote datasets should be, restricted to them) |
| HG008 only | per row: `series`, `platform`, `build`, `tumor_sample`, `normal_sample`, `tumor_coverage`, `normal_coverage` (`tumor_bai` etc. are `.crai` for CRAMs); per dataset: `reference_gzi`, `truth_sv`, `truth_sv_regions`, `truth_sv_regions_clonal` |

`accepts` can also restrict `sources: [seqc2]` or `[hg008]`. The shipped `variantcomp`
profile accepts only SEQC2 WGS/WES; `variantinterpretation` takes HG008 MuTect2 VCFs too
(sample names `HG008_<pair>_...`).

An empty variable gives `n/a` in a samplesheet column (configurable: `samplesheet.missing`),
leaves a param out, and is an error in an extra file. Rendered params `true` / `false`
become booleans.

Next to `params.yaml`, `params.config` holds the same params as a Nextflow config, for `-c`
or for `includeConfig` in a test profile of the workflow. Local files of the workflow folder
are written as `"<config base>/<path>"`; the base is the absolute folder path, or the
Nextflow expression given with `--config-base`, e.g. `'${projectDir}/assets/seqc2/test'` for
files committed into a pipeline, or `"${params.pipelines_testdata_base_path}variantinterpretation"`
for a dataset folder that is moved. `seqc2 run` writes the workflow files, a params file with absolute paths
into the run folder and runs `nextflow run <pipeline> -r <revision> -profile <profile>
-params-file ... --outdir <run>/results -resume`; further nextflow arguments follow `--`.

## Development

```bash
pytest                        # unit tests (no network)
black src tests && ruff check src tests
```

Git: CUBI process (`main` releases, `dev`, `feat-*` branches merged fast-forward).

## Citation

SEQC2 data: Fang LT, Zhu B, Zhao Y, et al. Establishing community reference samples, data
and call sets for benchmarking cancer mutation detection using whole-genome sequencing.
*Nat Biotechnol* 2021;39:1151-1160 (doi:10.1038/s41587-021-00993-6). Xiao W, Ren L, Chen
Z, et al. Toward best practice in cancer mutation detection with whole-genome and
whole-exome sequencing. *Nat Biotechnol* 2021;39:1141-1150 (doi:10.1038/s41587-021-00994-5).
HG008 data: "Development and extensive sequencing of a broadly-consented Genome in a Bottle
matched tumor-normal pair for somatic benchmarks", *Sci Data* (doi:10.1038/s41597-025-05438-2);
the benchmarks: see their READMEs on the GIAB server.
GIAB stratifications: Dwarshuis N, et al. *Nat Commun* 2024;15:9029. GENCODE: Frankish A,
et al. *Nucleic Acids Res* 2023;51:D942.
