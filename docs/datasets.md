# Dataset configs and dataset folders

A dataset config is a YAML file that says which reference data a test needs. `somref plan`
shows what it selects, `somref fetch` builds the dataset folder, `somref samplesheet` /
`somref run` turn the folder into the inputs of a workflow ([profiles.md](profiles.md)).

Every key is optional; everything left out takes the default. Unknown keys, unknown values
(e.g. a misspelled caller) and combinations that the source doesn't provide (e.g. SEQC2 NS
pairs with WES, HG008 `illumina-nygc` with N-P) stop with an error that names the allowed
values.

Examples: [`datasets/`](../datasets) holds general examples (`default.yaml`,
`hg008_default.yaml`, `hg008_small.yaml`); the configs of the shipped workflow profiles
are in `datasets/workflows/<workflow>/` (see the [README](../README.md#workflow-specific-profiles-and-configs)).

## SEQC2 config (`source: seqc2`, the default)

```yaml
source: seqc2             # seqc2 (default) or hg008 (next section)
access: local             # local: download and build (somref fetch); remote: see below
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
  reference_plain: false  # contig subsets as plain FASTA (default: bgzip + .gzi)
  target_bed: true        # exome target BED (when WES is selected); remote: local / false
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
cache: null               # default: $SOMREF_CACHE or ~/.cache/somref
```

**Default (empty config):** WGS of one center and replicate, IL_1, aligned with bwa and
called with MuTect2; VCF, reference and truth set, genome-wide.

`vcf` is for workflows that expect a sample order (e.g. filters that index the samples by
position) or a read depth per sample: the SEQC2 MuTect2 (GATK 3.8) records have no
FORMAT/DP. Records with FORMAT/DP are kept; the header documents the change
(`##somatic_refdata_vcf=...`).

## HG008 config (`source: hg008`)

```yaml
source: hg008
select:
  series: [illumina-nygc]   # sequencing series (docs/sources.md), or all
  normals: [N-D]            # N-D and/or N-P: pairs T-ND, T-NP
  build: GRCh38             # GRCh38 (GIABv3) or GRCh37 (hs37d5)
  callers: [muTect2]        # for the selected series that have them; [] for none
data:
  vcf: true
  bam: false                # false, true or a list of pairs (T-ND, T-NP)
  reference: true           # GIAB reference of the build (+ .fai, .gzi)
  truth: v0.4               # small-variant draft benchmark, or false
  truth_tier: all           # all, nogermlineinterference, nogermlinewithin50bp
  truth_sv: false           # SV draft benchmark v0.6 (GRCh38 only), or false
  gencode: false            # GENCODE basic GTF release (GRCh38 only), or false
  giab: []                  # GIAB stratifications of the build
subset: {}                  # as for SEQC2 (contigs, bam_windows, bam_subsample; no targets)
```

**Default (`source: hg008` only):** T vs N-D, NYGC Illumina, MuTect2 VCF, GRCh38 reference
and the small-variant benchmark, genome-wide. A caller that no selected series has (e.g.
`muTect2` without `illumina-nygc`) is an error: set `callers` or `data.vcf: false`. The ONT
tumor series take their normals from the series `ont-ne`.

`datasets/hg008_small.yaml` is a multi-technology example: chr22, Illumina, HiFi, ONT
ultra-long and Ultima slices of 100 kb, MuTect2 and DeepSomatic VCFs, both benchmarks and
the GIAB difficult regions, about 60 MB.

## Choosing the samples on the command line (SEQC2)

`plan`, `fetch`, `samplesheet`, `run` and `verify` take `--technologies`, `--pairs`,
`--centers` and `--replicates` (comma-separated, or `all`); they replace the values of the
config's `select` section. `--centers` / `--replicates` clear the config's `pairs`. Use the
same options for `fetch` and `samplesheet`: `samplesheet` stops if the dataset folder was
built with another selection.

```bash
somref fetch datasets/workflows/variantcomp/offline.yaml --pairs IL_1,FD_1,NV_1   # three pairs
somref fetch datasets/workflows/variantcomp/offline.yaml --centers all            # first replicate of every center
somref fetch datasets/workflows/variantinterpretation/offline.yaml --centers all --replicates all --technologies WES,WGS
```

## Remote datasets (`access: remote`, SEQC2)

Nothing is downloaded (except the server's md5 lists for the catalog). `somref fetch` (or
`somref samplesheet` directly) writes only `dataset.json` and `README.md`; every file in the
samplesheet and params is a URL that the workflow reads itself (Nextflow stages it):

- VCFs: the original SEQC2 files (plain gzip, not bgzip, no index; sample columns TUMOR,
  NORMAL; all contigs). The workflow must accept plain gzip.
- Reference: without `subset.contigs` the whole GRCh38.d1.vd1 FASTA (3.2 GB). With
  `subset.contigs` a smaller online reference with exactly these contigs
  (`remote_references` in `sources_seqc2.yaml`; currently `chr21` from nf-core/test-datasets,
  sequence identical to GRCh38.d1.vd1). The VCFs still hold every contig: the workflow has
  to drop the other records (template variable `{contigs}`, e.g. `bcftools view --targets`).
- Target BED, truth set, GENCODE GTF, GIAB files: the original files (URLs).
  `data.target_bed: local` instead writes the WES target BED into the dataset folder (BED3
  of the reference contigs, sorted, merged, gzip: `bed/<kit>_hg38.bed.gz`; chr21: 24 KB,
  genome: 2.2 MB), for workflows that need a `.bed`/`.bed.gz` (the remote file is a BED6
  named `*.bed6.gz`).
- Not possible: BAM slices (`bam_windows`, `bam_subsample`); whole BAMs are possible
  but 10-160 GB each.

## How much data does a workflow need?

| Question | Recommendation |
|---|---|
| Replicates (SEQC2) | One per center (`replicates: [1]`). Replicates of one center are the same DNA, library kit and sequencer; they only matter for reproducibility studies. Several centers show real technical variability. |
| Technologies | SEQC2: WGS for genome-wide tests; WES for capture workflows (with `target_bed`); both for cross-technology comparisons. HG008: one short-read series (e.g. `illumina-nygc`) plus long reads (`hifi-bcm`, `ont-ul-ucsc`) for multi-platform tests. |
| Aligners / callers | One (bwa + muTect2) unless the workflow compares callers or aligners. Check which FORMAT fields your workflow reads ([sources.md](sources.md)). |
| BAMs | Only for workflows that read alignments (coverage, read evidence, calling). With `subset.bam_windows` the toolkit reads small slices remotely instead of downloading 10-420 GB per file. |
| Size | For CI / laptops: `contigs: [chr22]`, BAM windows of 0.1-2 Mb, `bam_subsample: 0.5`: about 150 MB for 3 SEQC2 pairs x 2 technologies, about 60 MB for `hg008_small.yaml`. |

## Downloads, cache and checksums

- Every remote file is downloaded once into the cache
  (`<cache>/<host>/<path>`, resumable) and checked against the md5 of the server's lists
  (SEQC2 VCFs, BAMs, BAM indices; HG008 alignments, VCFs and benchmark files; GIAB
  stratifications). `<cache>/registry.json` records size, time and md5 of every cached
  file, so it is not hashed again unless it changed (`--verify` hashes everything again).
- The dataset folder holds the derived files (bgzipped and indexed VCFs, subsets, BAM
  slices, ...). `somref_manifest.json` records the recipe (source, md5, settings) and md5
  of each; `somref fetch` rebuilds only files whose recipe changed or that were modified,
  and removes files that the config no longer selects.
- `checksums.md5`, `README.md` (sources, settings, citation) and `dataset.json` (the
  dataset for profiles) are written into the dataset folder. `somref verify` checks every
  file.
- HTTPS: downloads are done in Python; the CA bundle is taken from `--ca-bundle`,
  `$CURL_CA_BUNDLE`, `$SSL_CERT_FILE`, the Python package certifi or the system paths, and
  exported for samtools/bcftools reading remote files. This avoids
  `curl: (77) error adding trust anchors` with conda builds of curl.

## Dataset folder

SEQC2:

```
reference/GRCh38.d1.vd1[_<contigs>].fa[.gz]  (+ .fai, .gzi)
bed/S07604624_Covered_human_all_v6_plus_UTR_hg38.bed
vcf/<pair>_<technology>.<aligner>.<caller>.vcf.gz  (+ .tbi; sample columns TUMOR, NORMAL)
bam/<pair>_<technology>_<T|N>.bam  (+ .bam.bai)
truth/high-confidence_sSNV|sINDEL_in_HC_regions_<release>.vcf.gz, High-Confidence_Regions_*.bed
annotation/gencode.v<release>.basic.annotation.gtf.gz (+ .tbi)
annotation/giab/GRCh38/<GIAB release path>
dataset.json, checksums.md5, README.md, somref_manifest.json
<profile>/samplesheet.csv, params.yaml, params.config, ...   (somref samplesheet / run)
```

HG008, where it differs:

```
reference/GRCh38_GIABv3[_<contigs>].fa.gz (+ .fai, .gzi)   or hs37d5 for GRCh37
vcf/<pair>_<series>.<aligner>.<caller>.vcf.gz              pair: T-ND, T-NP
bam/HG008-<T|N-D|N-P>_<series>.bam|cram (+ .bai / .crai)   one file per sample and series
truth/HG008-T_somatic_smvar_benchmark_<release>_<build>_tumorvariants_snps|indels.vcf.gz
truth/HG008-T_somatic_smvar_benchmark_<release>_<build>_<tier>.bed
truth_sv/<build>_HG008-T-V0.6_somatic-stvar_PASS.draftbenchmark.vcf.gz, *.bed
annotation/giab/<build>/<GIAB release path>
```
