# Oracle ChIP-seq Analysis and Differential Binding Pipeline

A version-constrained `Snakemake` workflow for paired-end transcription-factor and histone ChIP-seq. It performs read QC and contamination screening, alignment and duplicate removal, explicit narrow- or broad-peak calling, factor-specific reproducibility filtering, normalized signal-track generation, `deepTools` visualization and sample-level QC, `DiffBind` differential binding, HOMER motif enrichment, and consolidated reporting for human, mouse, or rat data.

## Repository structure
- `Snakefile` - main Snakemake workflow
- `workflow/` - stage-specific Snakemake rules and environment smoke checks
- `config.yaml` - species references, ChIP-seq sample definitions, contamination settings, and design metadata
- `config.sample.yaml` - example config template with a valid two-condition, two-replicate DiffBind design
- `config/fastq_screen.conf.example` - template FastQ Screen contamination database config
- `sample_manifest.tsv` - example manifest for config generation
- `envs/` - Conda environment definitions for ChIP-seq processing and R analysis
- `scripts/` - helpers for building DiffBind sample sheets, reference downloads, and manifest-driven configs
- `analysis/` - R analysis scripts for DiffBind and motif summaries
- `setup.sh`, `run.sh`, `make_config.sh`, and `prepare_references.sh` - activation-free entry points
- `TUTORIAL.md` - beginner-oriented installation, configuration, execution, and interpretation guide
- `docs/ANALYSIS_AND_VISUALIZATION.md` - analysis rationale, configuration, output map, interpretation, and focused target commands
- `VALIDATION.md` - reproducible repository checks and known validation limits
- `tests/` - design, helper, setup, and synthetic-DAG tests
- `.gitignore` - files and folders excluded from Git tracking
- `LICENSE` - MIT license

## What this pipeline does
1. Raw and trimmed FASTQ QC with `FastQC`, plus raw-read contamination screening with `FastQ Screen`
2. Paired-end adapter trimming with `Trim Galore`
3. Paired-end `Bowtie2` alignment with configurable MAPQ and proper-pair filtering
4. Name/coordinate sorting, mate-tag repair, duplicate removal, indexing, `flagstat`, and duplicate-removal metrics with `samtools`
5. Per-sample `MACS3` narrow- or broad-peak calling from the explicit `peak_mode`
6. Blacklist removal, replicate-supported peaks within each factor/condition, and a factor-level union of those condition BEDs
7. RPGC coverage BigWigs for ChIP and input samples, plus matched log2 ChIP/input BigWigs
8. Factor-specific `computeMatrix`, `plotHeatmap`, and `plotProfile` outputs in PNG, PDF, and tabular forms
9. Factor-specific fingerprint, peak-signal Spearman-correlation, and PCA diagnostics
10. Factor-specific `DiffBind` differential binding and narrow-factor HOMER motif enrichment
11. Project-wide `MultiQC` and optional Snakemake HTML reports

## Scope
This repository is intentionally ChIP-seq only. It does not run RNA-seq quantification, differential expression, or transcriptomics integration. The statistical comparison step is differential binding analysis from ChIP-seq peaks and aligned reads.

For the reasoning behind the defaults, exact output paths, and plot interpretation, see [Analysis and visualization](docs/ANALYSIS_AND_VISUALIZATION.md).

## Quick start

From the repository root, use the same three commands for installation,
validation, and execution:

```bash
bash setup.sh
bash run.sh --dry-run
bash run.sh --cores 8
```

The installer is data-independent: it installs the small
`oracle-chip-runner`, both rule environments, and runs package smoke
checks without reading FASTQs or references. Configure the files in steps 1–3
below before the dry-run. The complete beginner end-to-end guide is in
[TUTORIAL.md](TUTORIAL.md).

The environment specifications resolve for Linux x86_64 and Linux aarch64.
Installation runs package smoke checks on the current host; native aarch64
execution is not covered by this repository's CI. To install Miniforge outside
`$HOME/miniforge3`, export the same custom prefix for setup and later
entrypoints, for example
`export MINIFORGE_HOME=/opt/miniforge3`.

## Configure the workflow
### 1. Fill in `config.yaml`
- Set `species` to `human`, `mouse`, or `rat`.
- Provide absolute paths for:
  - `genome`
  - `chrom_sizes`
  - `black_list`
  - `bt2_index`
- Set `gsize` for MACS peak calling and a separate, literal positive integer `effective_genome_size` for deepTools RPGC normalization. Because alignments are MAPQ-filtered, choose a read-length- and mappability-specific value; the examples assume 150-bp reads. Do not enter the latter in scientific notation.
- Configure `contamination.fastq_screen_conf` to point to a FastQ Screen config file with Bowtie2 index prefixes for the genomes or contaminants you want to screen.
- Review `alignment.min_mapq`, `alignment.max_insert_size`, and the settings under `peak_calling`, `differential_binding`, and `deeptools`. Proper concordant pairs are an invariant of this paired-end workflow; the default maximum insert size is 2,000 bp.
- Update ChIP sample names, exactly two unique paired gzipped FASTQ paths, input controls, conditions, replicates, histone mark/factor, tissue labels, and an explicit `peak_mode` of `narrow` or `broad`. Give each input control its biological `condition`; a ChIP sample may only reference a control from the same condition. FASTQ filenames need not already follow the sample naming convention; trimming stages them under deterministic temporary names.
- Keep one peak mode per factor. The supplied H3K27ac example is `narrow`, consistent with its punctate enrichment geometry.
- For every factor, configure exactly two conditions with at least two biological replicates per condition. Replicate numbers must be unique within each factor/condition group.
- Set `differential_binding.numerator_condition` and `reference_condition` to those two labels. The supplied simple model also requires one tissue label per factor; use an extended design for blocked, paired, batch-aware, or multi-tissue studies.

### 2. Configure contamination screening
Copy and edit the FastQ Screen template:
```bash
cp config/fastq_screen.conf.example config/fastq_screen.conf
```
Replace each `/path/to/...` database entry with a real Bowtie2 index prefix, then set this in `config.yaml`:
```yaml
contamination:
  fastq_screen_conf: "config/fastq_screen.conf"
  subset: 100000
```

### 3. Place raw data
Put FASTQ pairs in `data/raw/` or update paths in `config.yaml`.

### 4. Install the runner and rule environments
```bash
bash setup.sh
```

### 5. Run the full workflow
```bash
bash run.sh --dry-run
bash run.sh --cores 12
```

### 6. Run selected downstream analyses
Prefer Snakemake targets so the correct rule environments are selected automatically. Factor names are converted to filesystem-safe slugs; the example factor remains `H3K27ac`:

```bash
bash run.sh --cores 4 results/diffbind/H3K27ac/diffbind_summary.csv
bash run.sh --cores 4 results/motifs/H3K27ac/motif_summary.pdf
```

## How to use the pipeline
- `bash run.sh --cores 12` builds all core outputs in `results/`.
- Append `results/contamination/fastq_screen/H3K27ac_control_rep1_R1_screen.txt` to run one contamination check.
- Append `results/diffbind/sample_sheet.csv` to generate the DiffBind sample sheet.
- Append `results/peaks/consensus/H3K27ac/control.bed` to build one replicate-supported factor/condition BED.
- Append `results/peaks/consensus/H3K27ac.bed` to union both condition BEDs for factor-level visualization.
- Append `results/deeptools/H3K27ac/heatmap.png` or `results/deeptools/H3K27ac/profile.png` to make one primary visualization.
- Append `results/deeptools/H3K27ac/qc/spearman_heatmap.png` or `results/deeptools/H3K27ac/qc/pca.png` to make one cross-sample diagnostic.
- Append `results/diffbind/H3K27ac/diffbind_summary.csv` to run factor-specific differential binding.
- Append `results/motifs/H3K27ac/motif_summary.pdf` to run motif enrichment for a narrow-peak factor. Motif analysis is not a default output for broad factors.
- Append `results/multiqc/multiqc_report.html` to generate the MultiQC report.
- Append `results/report/snakemake_report.html` to generate the Snakemake HTML report.

See [Analysis and visualization](docs/ANALYSIS_AND_VISUALIZATION.md) for the complete output contract and more focused targets.

## Sample manifest and config generation
1. Edit `sample_manifest.tsv` with your sample names and file paths.
2. Generate `config.yaml` from the manifest:
```bash
bash make_config.sh sample_manifest.tsv \
  --species human \
  --genome /path/to/hg38.fa \
  --chrom-sizes /path/to/hg38.chrom.sizes \
  --bt2-index /path/to/hg38 \
  --black-list /path/to/hg38-blacklist.v2.bed \
  --effective-genome-size 2862010428 \
  --fastq-screen-conf config/fastq_screen.conf \
  --numerator-condition treated \
  --reference-condition control \
  --output config.yaml
```

The effective size above is the deepTools example for uniquely mappable 150-bp
GRCh38 reads. Replace it for a different assembly, read length, or mapping
policy. Change the two condition options when the manifest uses labels other
than `treated` and `control`; they define the reported DiffBind fold-change
direction.

## Reference download helper
Use the species-specific helper to download a genome FASTA, build a Bowtie2 index, and generate chromosome sizes:
```bash
bash prepare_references.sh human --outdir references
```

Completed compressed downloads and decompressed FASTAs are validated before
reuse. Interrupted downloads remain as resumable `.part` files; completed data
are moved into place atomically, and a size-bound completion marker detects
legacy or subsequently truncated FASTAs before indexing.

Build separate Bowtie2 indexes for any extra contamination-screening databases you add to `config/fastq_screen.conf`, such as PhiX, E. coli, yeast, or alternative host genomes.

## Generate reports
```bash
bash run.sh --cores 12 results/multiqc/multiqc_report.html
bash run.sh --cores 12 results/report/snakemake_report.html
```

## Expected outputs
- `results/fastqc/` - raw and trimmed FASTQ QC reports
- `results/contamination/fastq_screen/` - raw FASTQ contamination-screen reports
- `results/trimmed/` - trimmed FASTQ files
- `results/bam/` and `results/qc/samtools/` - duplicate-removed BAMs, indices, `flagstat`, and `markdup` removal metrics
- `results/peaks/raw/` - raw MACS3 `narrowPeak` or `broadPeak` files
- `results/peaks/` - blacklist-filtered sample peaks, replicate-supported factor/condition BEDs, factor-level condition unions, and the cross-factor compatibility BED
- `results/bigwig/` - per-library RPGC tracks and per-ChIP log2 ChIP/input tracks
- `results/deeptools/<factor>/` - matrix, plotted-region BED, heatmap/profile PNG and PDF, numeric TSVs, fingerprints, Spearman correlations, and PCA
- `results/diffbind/<factor>/` - directional contrast metadata, differential-binding table, diagnostic PDF, and serialized DiffBind object
- `results/motifs/<factor>/` - HOMER and summary outputs for narrow-peak factors
- `results/provenance/` - resolved configuration used to build the workflow
- `results/multiqc/` - consolidated MultiQC report and parsed metrics
- `results/report/` - optional Snakemake HTML report

## Notes and best practices
- Use FASTQ and reference files that match the selected species assembly.
- Confirm `bt2_index` points to a built Bowtie2 index set.
- Confirm `black_list` points to the matching genome blacklist file.
- Confirm FastQ Screen database entries point to valid Bowtie2 index prefixes.
- Check the DiffBind sample sheet before running contrasts; each factor must have exactly two conditions, at least two biological replicates per condition, and no duplicate replicate numbers within a condition.
- Treat each factor/condition consensus BED as an at-least-N replicate-overlap rule and the factor BED as their union, not as IDR output.
- Relative library normalization cannot establish global biological signal shifts; use an appropriate spike-in design and a validated spike-in-aware analysis when that inference is required.
- This implementation accepts paired-end libraries only.
- Use `bash setup.sh --check` to verify the runner and both rule environments.
