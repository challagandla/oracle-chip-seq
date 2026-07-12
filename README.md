# oracle-chip-seq

Snakemake workflow for ChIP-seq differential binding. Handles histone modifications
and transcription factors, single-end and paired-end libraries, and selects narrow
or broad peak calling per target rather than per run.

## The peak-mode rule

Whether a target is punctate or domain-forming changes far more than one flag on
the peak caller, so the classification lives in one file,
`config/mark_registry.yaml`, and every downstream step reads it from there.

| | narrow | broad |
|---|---|---|
| examples | H3K4me3, H3K27ac, H3K9ac, H3K4me2, H2A.Z, CTCF and other TFs | H3K4me1, H3K27me3, H3K9me3, H3K36me3, H3K79me2, H4K20me1 |
| peak caller | `macs3 callpeak` default | `macs3 callpeak --broad --broad-cutoff` |
| reproducibility | IDR across replicates | naive overlap against pooled-replicate peaks |
| width filter | drop peaks above a few kb (artefacts) | allow domains up to 1 Mb |
| profile geometry | reference-point on the summit, or on the TSS for promoter marks | scale-regions across the domain body |
| counting window | fixed window on the called summit | the full called interval |
| size factors | median-of-ratios on peak counts | genome-wide background bins |
| motif enrichment | yes | no |
| FRiP and depth floors | per mark | per mark |

Classification follows the [ENCODE histone ChIP-seq spec](https://www.encodeproject.org/chip-seq/histone/).
Two points it gets right that are commonly got wrong:

- **Acetylation and H3K4me3 are narrow.** Calling H3K27ac with `--broad` fuses
  neighbouring enhancers into multi-kb blocks and destroys the summits that motif
  analysis and summit-centred counting depend on.
- **H3K4me1 is broad**, despite marking enhancers. Its signal is a wide, low
  shoulder around the nucleosome-depleted region, not a sharp peak.

Unlisted targets fall back to narrow, with a warning. Narrow is the safer default:
a punctate caller on a broad mark loses breadth, whereas a broad caller on a
punctate mark silently merges distinct regulatory elements.

Why each downstream step branches:

- **IDR only for narrow.** IDR models rank consistency within a ranked peak list.
  Broad peaks are wide with compressed scores and the model does not hold; ENCODE
  uses overlap for broad marks for exactly this reason.
- **Background-bin size factors for broad.** Median-of-ratios assumes most peaks
  are unchanged. A repressive domain mark can move globally — the textbook case is
  EZH2 inhibition — and when it does, normalising on reads-in-peaks absorbs the
  biology into the size factors and reports nothing. Estimating depth from
  background bins keeps the global component in the result. Where a truly global
  shift is expected an exogenous spike-in is the only unbiased normaliser;
  background bins are a proxy.
- **No motifs for broad.** Scanning a 50 kb Polycomb domain for 8-mers recovers its
  base composition, and the software will report that with confident p-values.
- **Per-mark QC thresholds.** A FRiP of 2% is a failure for H3K4me3 and entirely
  normal for H3K9me3. One global cutoff cannot serve both.

Inspect the resolved settings for any target:

```bash
python3 scripts/marks.py             # the whole table
python3 scripts/marks.py H3K27me3    # one target, fully resolved
```

## Steps

FASTQ (SRA or local) → FastQC → Trim Galore → bowtie2 → filter (MAPQ ≥ 30, primary
chromosomes, proper pairs, deduplicate) → fragment-length estimate → MACS3 (mode
from the registry) → blacklist and width filter → IDR or naive overlap →
per-target consensus → featureCounts → DESeq2 → ChIPseeker and GO →
monaLisa/JASPAR motifs → figures, browser tracks, MultiQC, QC gate.

Also computed: FRiP, NRF/PBC1/PBC2 library complexity, deepTools fingerprint,
log2(ChIP/Input) bigwigs, Spearman correlation and PCA on binned signal.

## Running it

```bash
snakemake --use-conda --cores 32
```

Reference paths and the contrast are set in `config.yaml`; the cohort is
`samples.tsv`:

```
sample_id  srr  target  assay  condition  replicate  control_id  layout  cell_line
```

`assay` is `chip` or `input`. `layout` is `single` or `paired` **per sample** —
libraries of both kinds can coexist in one run, and both are counted as fragments.
Leave `srr` blank and drop FASTQs into `data/raw/` to skip the download.

The design is validated before anything runs: missing controls, fewer than two
replicates in a condition, and unknown layouts all fail at parse time rather than
three hours into an alignment.

## Demo dataset

`samples.tsv` points at [GSE277460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277460)
(PRJNA1162470): Jurkat T cells, resting vs anti-CD3/CD28 stimulated, 28 libraries
over six targets — H3K4me3, H3K27ac (punctate histone), H3K4me1, H3K27me3 (domain
histone), CTCF (TF), POLR2A (polymerase) — with matched inputs, two conditions,
two biological replicates each.

It exercises every branch of the registry. It is also, as deposited, mostly
unusable, and that is the more instructive half:

- **SRA declares every replicate-1 run SINGLE. They are 2×150 paired** (92.7%
  concordant alignment). `samples.tsv` records the true layout, and
  `scripts/fetch_sra.py` fails loudly when the data disagree with the declaration
  rather than silently discarding R2.
- **Most ChIP libraries are contaminated with high-GC non-human DNA.** Alignment
  rate tracks GC content almost perfectly, and the human genome is ~41% GC:

  | library | GC | aligned |
  |---|---|---|
  | Input_resting_rep2 | 40.5% | 96.6% |
  | H3K4me3_resting_rep1 | 51.8% | 53.2% |
  | CTCF_resting_rep1 | 54.9% | 19.5% |
  | POLR2A_resting_rep2 | 56.3% | 12.7% |

  Library complexity is normal throughout (NRF 0.76–0.90), so this is not a PCR
  bottleneck — the reads are simply not human. The QC gate fails 20 of 28
  libraries on alignment rate and FRiP.

- **The CTCF and H3K27me3 ChIPs did not enrich at all** (FRiP ≈ 0.0000, 0–24
  peaks). Not a calling artefact: called against genome background alone, with no
  control and `--nolambda`, H3K27me3 yields 18 "domains" and CTCF 52 peaks, and
  the top hits sit in pericentromeric satellite.

- **Consequently the contrast is underpowered and returns almost nothing** — one
  significant region across the cohort, H3K4me3 at *MYC* (log2FC 2.27, padj 7e-4),
  which is at least the right gene for T-cell activation. Every target has at
  least one contaminated replicate. This is reported as a failed experiment, not
  dressed up as an absence of regulation.

The point of shipping it this way: a workflow is only worth as much as its
behaviour on bad data, and this dataset produces bad data in four different ways.

## Outputs

```
results/
  qc/            qc_gate.tsv (per-mark PASS/WARN/FAIL), multiqc, fingerprint, correlation, PCA
  bam/           filtered alignments
  bigwig/        *.cpm.bw and *.log2ratio.bw (ChIP/Input)
  peaks/         per-sample; reproducible/ (IDR or overlap); consensus/ per target
  profiles/      heatmaps and meta-profiles, geometry per peak mode
  differential/  DESeq2 results, up/down BEDs, size factors
  annotation/    ChIPseeker annotation, differential genes, GO
  motifs/        JASPAR enrichment in up/down peaks
  figures/       publication figures and genome-browser panels
  summary/       analysis_summary.md
```

Nothing under `results/` or `data/` is tracked by git.

## Licence

Pipeline code is MIT (`LICENSE`). No third-party code or reference data is bundled;
tools are conda-installed and invoked as separate processes. Every tool used is
open-source and free for commercial use with citation. HOMER is deliberately **not**
used — it is academic-use-only and not redistributable — so motif enrichment uses
monaLisa + JASPAR. See [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).
