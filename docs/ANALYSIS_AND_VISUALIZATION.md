# Analysis and visualization

This guide describes the analysis contract implemented by `Snakefile`. It explains why the workflow separates peak geometry, peak-caller genome size, signal-track normalization, and factor-specific visualization. Paths below use the example factor slug `H3K27ac`; substitute the slug generated from your configured `factor` value.

## Analysis overview

For every ChIP and input library, the workflow performs raw-read QC, contamination screening, paired-end trimming, alignment, filtering, duplicate removal, BAM indexing, and alignment summaries. It then branches as follows:

1. Each ChIP sample is peak-called against its matched input with the configured narrow or broad geometry.
2. Sample peaks are blacklist-filtered, reduced to replicate-supported peaks within each factor/condition, and then unioned across the two conditions for each factor.
3. Every ChIP and input library receives an RPGC coverage track; every ChIP sample also receives a matched log2 ChIP/input track.
4. The log2 tracks and factor consensus BED drive factor-specific matrices, heatmaps, and meta-profiles.
5. BAM fingerprints and peak-level signal summaries drive factor-specific fingerprint, Spearman-correlation, and PCA diagnostics.
6. Individual sample peaks and BAMs drive a separate DiffBind run for each factor.
7. Narrow-factor consensus peaks drive HOMER motif enrichment. Motif enrichment is not a default broad-factor analysis.

## Configuration that changes the analysis

### Peak geometry is explicit

Every record under `chip_samples` must contain one of:

```yaml
peak_mode: narrow
```

```yaml
peak_mode: broad
```

All samples assigned to one factor must use the same mode. The mode controls both the MACS output type and the default deepTools matrix geometry:

| `peak_mode` | MACS behavior | Peak file | Default matrix |
|---|---|---|---|
| `narrow` | q-value threshold, no `--broad` | `narrowPeak` | centered reference point |
| `broad` | `--broad --broad-cutoff` | `broadPeak` | scaled region body plus flanks |

H3K27ac is explicitly `narrow` in the supplied configuration because ENCODE classifies it as a narrow histone mark. A domain-like mark such as H3K27me3 should normally be configured as `broad`. The factor name alone never silently chooses a mode.

The relevant thresholds are:

```yaml
peak_calling:
  narrow_qvalue: 0.01
  broad_cutoff: 0.1
  consensus_min_replicates: 2
```

The peak-calling rules use duplicate-removed paired-end BAMs, the matched input BAM, `--format BAMPE`, and `--keep-dup all`; duplicates have already been removed upstream by `samtools markdup -r`.

### MACS genome size and deepTools effective genome size are separate

Each selected reference must define both values:

```yaml
references:
  human:
    gsize: "2.7e9"
    # Example for uniquely mappable 150-bp GRCh38 reads.
    effective_genome_size: 2862010428
```

- `gsize` is passed to MACS as `--gsize` when estimating the genome-wide background.
- `effective_genome_size` is passed to `bamCoverage --effectiveGenomeSize` for RPGC normalization. It must be a positive YAML integer, not text or scientific notation. Because this workflow applies a MAPQ filter, use a read-length-specific uniquely mappable value; non-N genome lengths are intended for analyses that retain multimappers.

They are separate configuration fields because they are consumed by different tools and have different input contracts. Select an assembly- and mapping-policy-appropriate value rather than copying a number between assemblies.

### Alignment filtering

The default alignment policy is:

```yaml
alignment:
  min_mapq: 30
  max_insert_size: 2000
```

Bowtie2 runs with `--no-mixed --no-discordant` and an explicit maximum insert size of 2,000 bp. The latter avoids Bowtie2's smaller implicit default silently excluding longer valid fragments and can be adjusted for a library with a justified fragment-length distribution. `samtools view` applies the MAPQ threshold, excludes flags represented by `-F 1804`, and always requires proper-pair flag `-f 2`. This is an invariant rather than a misleading toggle because the downstream fragment-level track and peak operations require concordant pairs. The workflow then repairs mate tags, coordinate-sorts, removes marked duplicates, indexes the final BAM, and writes `results/qc/samtools/<sample>.flagstat.txt`. `samtools markdup -s` separately records duplicate-removal counts in `results/qc/samtools/<sample>.markdup.txt`; post-removal `flagstat` alone cannot reconstruct that fraction.

This is a paired-end-only contract: every ChIP and control record must provide exactly two distinct gzipped FASTQ paths, and one FASTQ cannot be assigned to multiple sample records. Their original basenames can be arbitrary. Before trimming, the rule creates temporary symlinks named `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`; this gives Trim Galore deterministic output names without renaming or copying the source data. The temporary staging directory is removed when the rule exits.

## Peaks and reproducibility

### Per-sample peaks

MACS outputs are written to:

```text
results/peaks/raw/<sample>_peaks.narrowPeak
results/peaks/raw/<sample>_peaks.broadPeak
```

Only the path matching the sample's `peak_mode` is built. Any region overlapping the assembly-matched blacklist is removed with `bedtools intersect -v`, producing:

```text
results/peaks/<sample>_peaks.narrowPeak
results/peaks/<sample>_peaks.broadPeak
```

Peak and consensus files are sorted against the configured two-column
`chrom_sizes` file. This makes reference order deterministic and causes unknown
chromosome names to fail early instead of being carried into matrices or tracks.

### Condition-aware consensus peaks

For each factor/condition group, `bedtools multiinter` finds intervals shared across that group's blacklist-filtered replicate peak files. Intervals present in at least `peak_calling.consensus_min_replicates` replicate peak sets are retained, sorted, and merged:

```text
results/peaks/consensus/<factor>/<condition>.bed
```

The workflow validates that the threshold does not exceed the configured replicate count in either condition. If no interval survives for a condition, the rule fails rather than emitting an apparently valid empty reference set.

The two condition-level BEDs are then concatenated, sorted, and merged into the factor-level union:

```text
results/peaks/consensus/<factor>.bed
```

For example, with two control and two treated H3K27ac replicates and a threshold of two, `H3K27ac/control.bed` requires support from both control replicates and `H3K27ac/treated.bed` requires support from both treated replicates. Their factor union retains loci reproducible within either condition as well as loci shared by both. A singleton peak from one control replicate cannot become "reproducible" merely by overlapping a singleton from one treated replicate.

The factor BEDs are also merged into:

```text
results/peaks/consensus_peaks.bed
```

That cross-factor BED is a compatibility/reporting artifact. deepTools and narrow-factor motifs use the factor-level condition union; DiffBind uses the factor's individual sample peaks and BAMs.

The per-condition consensus procedure is an at-least-N replicate-overlap rule, followed by a union across conditions. It is not IDR, does not preserve MACS scores or summits, and should not be reported as IDR-based reproducibility.

## Signal tracks

### RPGC coverage tracks

`bamCoverage` creates one track for every ChIP and input library:

```text
results/bigwig/<sample>.rpgc.bw
```

The command uses RPGC normalization, the configured effective genome size and bin size, read extension, first mates only (`--samFlagInclude 64`) to avoid double-counting paired fragments, the assembly blacklist, and the configured deepTools thread count. These tracks are useful for genome-browser inspection and factor-level correlation/PCA.

### Matched log2 ChIP/input tracks

`bamCompare` creates one track for every ChIP sample:

```text
results/bigwig/<chip-sample>.log2ratio.bw
```

Each ChIP BAM is compared with its configured matched control using `--operation log2`, the configured pseudocount, bin size, read extension, first-mate selection, and blacklist. With the supplied setting:

```yaml
deeptools:
  log2_pseudocount: 1
  log2_scale_method: "None"
```

the workflow disables a pairwise scaling-factor method and independently CPM-normalizes the ChIP and input before calculating the log2 ratio. The accepted alternative scaling methods are `readCount` and `SES`. The factor heatmaps and profiles intentionally use these control-adjusted log2 tracks, not the RPGC tracks.

In a log2 ChIP/input plot, zero denotes equal normalized signal, positive values denote relative ChIP enrichment, and negative values denote relatively greater input signal. This remains a relative library-normalized quantity.

## Factor-specific heatmaps and profiles

For each factor, `computeMatrix` combines its ChIP log2-ratio BigWigs over its consensus BED.

- Narrow factors use `reference-point --referencePoint center`, with the configured upstream and downstream windows. The default is 3 kb on each side.
- Broad factors use `scale-regions`, with the configured region-body length and flanks. The default is a 5 kb scaled body plus 3 kb on each side.
- `deeptools.factor_modes` can optionally override display geometry for a named factor with `reference_point` or `scale_regions` (`narrow` and `broad` are accepted aliases). Keys are original factor names, not filesystem slugs. For example:

  ```yaml
  deeptools:
    factor_modes:
      H3K27ac: reference_point
      H3K27me3: scale_regions
  ```

  This changes visualization geometry, not peak calling.

The reusable and inspectable matrix outputs are:

```text
results/deeptools/<factor>/matrix.gz
results/deeptools/<factor>/matrix.tsv
results/deeptools/<factor>/regions.bed
```

`matrix.gz` is the deepTools binary input to downstream plots, `matrix.tsv` contains numeric matrix values for auditing or custom analysis, and `regions.bed` records the plotted regions in the same descending-mean order as the matrix. Sorting happens once in `computeMatrix`; both plotting commands preserve that order. `computeMatrix` assigns sample labels, while the factor region label is supplied later to `plotHeatmap` and `plotProfile`, where deepTools supports `--regionsLabel`.

`plotHeatmap` produces:

```text
results/deeptools/<factor>/heatmap.png
results/deeptools/<factor>/heatmap.pdf
```

Rows retain the descending-mean order established by `computeMatrix`. The configured diverging color map and z-axis bounds make depletion and enrichment around zero visually comparable; the configuration requires the bounds to straddle zero.

`plotProfile` produces:

```text
results/deeptools/<factor>/profile.png
results/deeptools/<factor>/profile.pdf
results/deeptools/<factor>/profile.tsv
```

The profile is a mean summary and can hide distinct peak subclasses. Inspect the heatmap, matrix TSV, and region BED before interpreting a smooth average as a uniform response. PNGs use the configured DPI; PDFs are preferable for scalable figures.

## Factor-specific visualization QC

All QC paths live under `results/deeptools/<factor>/qc/`.

### Fingerprint

`plotFingerprint` receives every ChIP BAM for the factor plus its unique matched controls. It reapplies the configured minimum MAPQ, excludes blacklisted intervals, extends paired fragments, counts first mates, and writes:

```text
fingerprint.png
fingerprint.tsv
fingerprint_metrics.tsv
```

Use the curve and quality table to compare enrichment and library behavior between ChIP and input samples. A fingerprint is a diagnostic, not a substitute for replicate reproducibility, FRiP, antibody validation, or biological review.

### Peak-signal correlation and PCA

`multiBigwigSummary BED-file` summarizes the factor's ChIP RPGC tracks over its consensus peaks:

```text
peak_summary.npz
peak_signal.tsv
```

The same summary drives:

```text
spearman_heatmap.png
spearman.tsv
pca.png
pca.tsv
```

`plotCorrelation` uses Spearman correlation without discarding zero-only regions. `plotPCA` runs with `--transpose` on the RPGC peak-signal matrix. The workflow deliberately omits deepTools 3.5.x `--log2` in transposed mode because that release transposes a pre-transform copy and silently ignores the requested transform. Replicates should be assessed for unexpected disagreement, condition structure, and isolated samples, while remembering that both plots are restricted to a peak set derived from the same experiment. Review raw QC, BAM metrics, browser tracks, and experimental metadata before labeling a sample an outlier.

## Factor-specific DiffBind and motifs

### DiffBind

Contrast direction and narrow-peak recentering are explicit:

```yaml
differential_binding:
  numerator_condition: treated
  reference_condition: control
  narrow_summits: 200
```

Every factor must use those two condition labels. Reported fold changes therefore have a stable treated-versus-control direction rather than depending on sample-sheet row order. Change both labels when your experiment uses different condition names.

`results/diffbind/sample_sheet.csv` contains all configured ChIP samples, their matched controls, individual sample peaks, factors, conditions, and replicate numbers. A separate R process subsets this sheet to each factor and writes:

```text
results/diffbind/<factor>/diffbind_summary.csv
results/diffbind/<factor>/diffbind_plots.pdf
results/diffbind/<factor>/diffbind.rds
results/diffbind/<factor>/contrast.tsv
```

The design contract requires exactly two conditions per factor, at least two ChIP replicates per condition, unique replicate numbers within each factor/condition group, one explicit tissue label per factor, and a control carrying the same condition label as each ChIP sample. This yields one unambiguous simple condition contrast. The script runs `dba.count(minOverlap = 2)`, creates the explicit numerator-versus-reference contrast with at least two members, analyzes it with DiffBind defaults, and reports sites at a 0.05 threshold. `contrast.tsv` and columns in the summary CSV record that direction. DiffBind uses the individual blacklist-filtered sample peak files and BAMs rather than the visualization consensus BED.

For broad factors, `dba.count` uses `summits = FALSE` so broad domains retain their called widths. Narrow factors explicitly use the configured `narrow_summits` value (200 gives a 401 bp recentered interval in DiffBind). This prevents DiffBind's narrow-oriented default from silently collapsing every broad domain to a short summit window.

The PDF contains overview, PCA, heatmap, and MA plots. Inspect the contrast definition and sample sheet before interpreting the CSV, especially for designs with batches, pairing, or donor effects; those covariates are not modeled by the supplied simple condition contrast. Designs with more or fewer than two conditions are rejected and require an explicitly extended analysis model.

### HOMER motifs

Motif enrichment is available only for factors configured with `peak_mode: narrow`; the HOMER rule's wildcard is constrained to those factor slugs, so a broad-factor motif target is rejected. HOMER analyzes the narrow factor's consensus BED with a 200 bp window and motif lengths 8, 10, and 12, then writes:

```text
results/motifs/<factor>/homer/knownResults.txt
results/motifs/<factor>/motif_summary.csv
results/motifs/<factor>/motif_summary.pdf
```

The workflow deliberately does not run default motif enrichment across broad domains, where a single fixed local sequence window is usually difficult to interpret. Motif enrichment is associative and should be validated against factor biology, background construction, and independent evidence.

The summary ranks motifs using HOMER's finite natural-log P-value column,
converted to `-log10(P-value)`. This preserves extremely significant motifs
whose printed P-values underflow to zero in ordinary double precision. HOMER
headers contain `#` characters in target/background count fields, so the TSV is
read without treating `#` as a comment marker.

## Useful targets

Run these from the repository root. Replace example sample and factor names with those in `config.yaml`.

```bash
# Validate configuration and the planned DAG without computing
bash run.sh --dry-run

# Full rule-all workflow
bash run.sh --cores 12

# One blacklist-filtered peak set, one condition consensus, and the factor union
bash run.sh --cores 8 results/peaks/H3K27ac_control_rep1_peaks.narrowPeak
bash run.sh --cores 8 results/peaks/consensus/H3K27ac/control.bed
bash run.sh --cores 8 results/peaks/consensus/H3K27ac.bed

# Coverage and matched-control tracks
bash run.sh --cores 8 results/bigwig/H3K27ac_control_rep1.rpgc.bw
bash run.sh --cores 8 results/bigwig/H3K27ac_control_rep1.log2ratio.bw

# Matrix, primary plots, and exported numeric data
bash run.sh --cores 8 results/deeptools/H3K27ac/matrix.tsv
bash run.sh --cores 8 results/deeptools/H3K27ac/heatmap.pdf
bash run.sh --cores 8 results/deeptools/H3K27ac/profile.tsv

# Factor QC
bash run.sh --cores 8 results/deeptools/H3K27ac/qc/fingerprint.png
bash run.sh --cores 8 results/deeptools/H3K27ac/qc/spearman_heatmap.png
bash run.sh --cores 8 results/deeptools/H3K27ac/qc/pca.png

# Statistical and motif analyses
bash run.sh --cores 8 results/diffbind/H3K27ac/diffbind_summary.csv
bash run.sh --cores 8 results/motifs/H3K27ac/motif_summary.pdf

# Consolidated and execution reports
bash run.sh --cores 12 results/multiqc/multiqc_report.html
bash run.sh --cores 4 results/report/snakemake_report.html
```

Rule logs are written under `results/logs/`. The Snakemake HTML report is optional and is not part of the default `rule all` target.

## Interpretation limits

- The workflow accepts paired-end data only; it is not a single-end ChIP-seq workflow.
- Condition consensus means overlap in at least N replicate peak sets within that factor/condition group. The factor BED is the union of its two condition BEDs. Neither step is IDR or establishes ENCODE pipeline compliance.
- The factor union deliberately contains condition-specific regions supported by sufficient within-condition replicates. It is not an intersection of peaks required to be present in both conditions.
- CPM, RPGC, and the supplied DiffBind defaults are relative normalization approaches. They cannot establish a global gain or loss of chromatin signal when most of the genome changes in one direction. Use an appropriate experimental spike-in and a validated spike-in-aware normalization strategy when global-shift inference is required.
- A matching genome assembly, chromosome naming convention, blacklist, Bowtie2 index, MACS genome size, and deepTools effective genome size are required. Mixed references invalidate the analysis.
- Scaling broad regions to one body length improves comparability but distorts absolute domain widths. Centering narrow regions depends on the accuracy and reproducibility of the called intervals.
- The current QC set does not implement every ENCODE metric or threshold. Fingerprints, correlations, and PCA are complementary diagnostics rather than a certification step.
- The supplied DiffBind model is an exactly-two-condition, one-tissue-per-factor analysis with an explicit contrast direction. Extend and validate the design before using it for paired, blocked, batch-confounded, multi-tissue, longitudinal, multicondition, or multifactor experiments.

## Primary method references

The rule environment pins deepTools 3.5.6 and NumPy 2.x. The upstream 3.5.6 documentation below was used for the command contracts and method guidance. The same environment pins MACS3 3.0.3, whose `callpeak` interface is documented by the linked upstream MACS3 reference.

- [deepTools 3.5.6 tool index](https://deeptools.readthedocs.io/en/3.5.6/content/list_of_tools.html), including direct documentation for [bamCoverage](https://deeptools.readthedocs.io/en/3.5.6/content/tools/bamCoverage.html), [bamCompare](https://deeptools.readthedocs.io/en/3.5.6/content/tools/bamCompare.html), [computeMatrix](https://deeptools.readthedocs.io/en/3.5.6/content/tools/computeMatrix.html), [plotHeatmap](https://deeptools.readthedocs.io/en/3.5.6/content/tools/plotHeatmap.html), and [plotProfile](https://deeptools.readthedocs.io/en/3.5.6/content/tools/plotProfile.html).
- [deepTools effective-genome-size guidance](https://deeptools.readthedocs.io/en/3.5.6/content/feature/effectiveGenomeSize.html), including read-length-specific values for MAPQ-filtered alignments.
- deepTools 3.5.6 documentation for [multiBigwigSummary](https://deeptools.readthedocs.io/en/3.5.6/content/tools/multiBigwigSummary.html), [plotCorrelation](https://deeptools.readthedocs.io/en/3.5.6/content/tools/plotCorrelation.html), [plotPCA](https://deeptools.readthedocs.io/en/3.5.6/content/tools/plotPCA.html), and [plotFingerprint QC metrics](https://deeptools.readthedocs.io/en/3.5.6/content/feature/plotFingerprint_QC_metrics.html).
- [ENCODE4 histone ChIP-seq standards and processing guidance](https://www.encodeproject.org/chip-seq/histone-encode4/), including biological replicates, matched controls, narrow/broad mark classes, read depth, and IDR-based reproducibility guidance.
- [ENCODE4 transcription-factor ChIP-seq standards](https://www.encodeproject.org/chip-seq/transcription-factor-encode4/), including replicate, control, antibody, sequencing-depth, and reproducibility expectations for TF experiments.
- [MACS3 `callpeak` documentation](https://macs3-project.github.io/MACS/docs/callpeak.html), the current upstream description of BAMPE, genome-size, q-value, broad-call, and output semantics used by the MACS call-peak interface.
- [DiffBind reference manual](https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf), including contrast construction and `summits` recentering semantics.
- [HOMER motif output documentation](https://homer.ucsd.edu/homer/motif/motifFinding.html), including the log P-value stored with motif enrichment results.
