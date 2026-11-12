# Reference data sources

The toolkit knows two public somatic reference datasets. A dataset config chooses one with
the top-level key `source` (default `seqc2`); `seqc2 catalog [--source hg008]` shows what the
server provides right now (built from the server's md5 lists).

| | SEQC2 (`source: seqc2`) | GIAB HG008 (`source: hg008`) |
|---|---|---|
| Tumor / normal | breast cancer cell line HCC1395 / B-lymphocyte line HCC1395BL | pancreatic cancer cell line HG008-T / normal duodenal (N-D) and pancreatic (N-P) tissue |
| Technologies | WGS and WES (Illumina) of 7 centers, up to 9 replicates | Illumina, Element AVITI, PacBio Onso, PacBio HiFi, ONT, Ultima (WGS) |
| Alignments | BAM (bwa), GRCh38.d1.vd1 | BAM / CRAM, GRCh38-GIABv3 and GRCh37 (hs37d5) |
| Somatic VCFs | 3 aligners x 3 callers per pair | MuTect2, Strelka2 (Illumina), DeepSomatic (HiFi) |
| Truth set | SEQC2 high-confidence somatic SNVs/indels v1.2, v1.2.1 | GIAB draft benchmarks: small variants v0.4, SVs v0.6 |
| Server | NCBI FTP `ReferenceSamples/seqc/Somatic_Mutation_WG` | NCBI FTP `ReferenceSamples/giab/data_somatic/HG008` |
| Description | `src/seqc2_toolkit/sources.yaml` | `src/seqc2_toolkit/sources_hg008.yaml` |

Both are public reference data; check the data use and citation notes of the providers
(see [Citation](../README.md#citation)).

## SEQC2 (HCC1395 / HCC1395BL)

The SEQC2 Somatic Mutation Working Group sequenced HCC1395 / HCC1395BL at several
centers, with several technologies, and called somatic variants with several pipelines
(Fang et al. 2021, Xiao et al. 2021). The toolkit covers the paired tumor/normal series:

| | WGS | WES |
|---|---|---|
| Pairs (center_replicate) | 21: EA_1, FD_1-3, IL_1-3, LL_1, NC_1, NS_1-9, NV_1-3 | 12: EA_1, FD_1-3, IL_1-3, LL_1, NC_1, NV_1-3 |
| VCFs | 3 aligners (bwa, bowtie, novo) x 3 callers (muTect2, strelka, somaticSniper) | same |
| BAMs | bwa, tumor + normal, every pair (10-160 GB each) | same |

Plus: the GRCh38.d1.vd1 reference, the exome target BED (Agilent SureSelect v6 + UTR), the
SEQC2 somatic truth set (v1.2, v1.2.1) with its high-confidence regions; optional GENCODE
and GIAB stratification files.

- Every SEQC2 VCF has the sample columns `TUMOR`, `NORMAL`. MuTect2 writes FORMAT/AD and
  AF (but no FORMAT/DP, see `vcf.format_dp` in [datasets.md](datasets.md)); Strelka only
  tier counts (AU/CU/GU/TU), SomaticSniper DP4/BCOUNT.
- Not in the catalog: the FFPE (FFG, FFX), library-prep (LBP), purity/depth (SPP) and
  AmpliSeq series, and FASTQs (SRA only).

## GIAB HG008 (HG008-T / HG008-N-D / HG008-N-P)

HG008 is the Genome in a Bottle (GIAB) cancer reference: the pancreatic ductal
adenocarcinoma cell line HG008-T (batch 0823p23) and two normal tissues of the same donor,
HG008-N-D (duodenal) and HG008-N-P (pancreatic). A dataset unit is a pair, `T-ND` or
`T-NP`, sequenced in one *series* (platform + center + run):

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

- **Alignments** (BAM; Ultima: CRAM) on the GIAB mapping references GRCh38-GIABv3 and
  GRCh37 (hs37d5), md5-checked against each folder's list. 30-300x; one file is
  50-420 GB, so read slices with `subset.bam_windows` (CRAMs are decoded with the remote
  reference and written as BAM).
- **Caller VCFs**: `muTect2`, `strelka2-snv`, `strelka2-indel` (NYGC pipeline,
  `illumina-nygc`, T vs N-D, GRCh38 + GRCh37) and `deepsomatic` (BCM, `hifi-bcm`, tumor
  column only, GRCh38). The sample columns differ per caller (`HG008-T`/`HG008-N`,
  `HG008.tumor`); profiles get them per VCF as `tumor_id` / `normal_id`.
- **Truth**: GIAB somatic small-variant draft benchmark v0.4 (GRCh38, GRCh37; one VCF,
  split by the toolkit into SNVs and indels; region tier `all` (recommended),
  `nogermlineinterference` or `nogermlinewithin50bp`) and the somatic SV draft benchmark
  v0.6 (GRCh38; PASS VCF, regions `clonal_and_subclonal` (for batch 0823p23 data) and
  `clonal`). **Both are drafts**: read their READMEs on the GIAB server before drawing
  conclusions.
- **Not in the catalog**: CHM13 alignments (two different CHM13 FASTAs, not on the GIAB
  server), series on other references (Ultima ppmSeq, BCM DRAGEN 2025, Hartwig/Roche
  Axelios, Dovetail), Hi-C, Bionano, single cell, other or superseded passages, NIST
  clones, FASTQs, further callers (ClairS, Severus, GRIDSS2, DRAGEN, ...).
- The GIAB server throttles parallel requests; the toolkit downloads one file at a time.
