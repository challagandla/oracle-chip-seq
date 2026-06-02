# Oracle Histone ChIP-seq + Transcriptomics Pipeline

A reproducible `Snakemake` workflow for histone ChIP-seq, RNA-seq transcription analysis, `deepTools` visualization, `DiffBind` differential binding, HOMER motif enrichment, and integrative peak-to-gene analysis across human, mouse, and rat.

## Repository structure
- `Snakefile` — main Snakemake workflow
- `config.yaml` — species references, sample definitions, design metadata
- `envs/` — Conda environment definitions for bioinformatics and R analysis
- `scripts/` — Python helper for building sample sheets
- `analysis/` — R analysis scripts for DiffBind, DESeq2, integrative analysis, and motif summaries
- `.gitignore` — files and folders excluded from Git tracking

## What this pipeline does
1. Raw FASTQ QC with `FastQC`
2. Adapter trimming with `Trim Galore`
3. Paired-end alignment with `Bowtie2`
4. BAM sorting, indexing, and duplicate removal with `samtools`
5. Broad peak calling for histone marks with `MACS2`
6. BigWig coverage track generation with `deepTools`
7. Heatmap/profile generation with `computeMatrix`, `plotHeatmap`, `plotProfile`
8. Differential binding analysis with `DiffBind`
9. Motif enrichment analysis with `HOMER`
10. RNA quantification with `Salmon`
11. Differential expression analysis with `DESeq2`
12. Gene/peak integrative analysis and scatterplots

## Quick start
### 1. Fill in `config.yaml`
- Set `species` to `human`, `mouse`, or `rat`
- Provide absolute paths for:
  - `genome`
  - `annotation`
  - `transcriptome`
  - `chrom_sizes`
  - `bt2_index`
  - `black_list`
- Update sample names, FASTQ paths, controls, conditions, and replicates.

### 2. Place raw data
Put FASTQ pairs in `data/raw/` or update paths in `config.yaml`.

### 3. Create Conda environments
```bash
conda env create -f envs/chipseq.yaml
conda env create -f envs/r_analysis.yaml
# Optional: dedicated RNA environment
conda env create -f envs/rna_seq.yaml
```

### 4. Run the full workflow
```bash
conda activate chipseq
snakemake --use-conda --cores 12
```

If you already have Snakemake installed in a separate environment, for example:
```bash
/home/epigenetics/miniforge3/envs/snakemake/bin/snakemake --use-conda --cores 12
```

### 5. Run only R analysis steps
```bash
conda activate r_analysis
Rscript analysis/diffbind_analysis.R results/diffbind/sample_sheet.csv results/diffbind
Rscript analysis/rnaseq_DE.R results/rnaseq/sample_metadata.tsv results/rnaseq results/rnaseq/deseq2_results.tsv results/rnaseq/deseq2_plots.pdf
Rscript analysis/integrative_analysis.R results/diffbind/diffbind_summary.csv results/rnaseq/deseq2_results.tsv /path/to/annotation.gtf results/integrative
```

## How to use the pipeline
- `snakemake --use-conda --cores 12` builds all outputs in `results/`
- `snakemake results/diffbind/sample_sheet.csv` generates sample metadata files
- `snakemake results/peaks/consensus_peaks.bed` creates a merged peak set for deepTools
- `snakemake results/integrative/integrative_scatter.pdf` runs the final integrative step

## Expected outputs
- `results/fastqc/` — QC reports
- `results/trimmed/` — trimmed FASTQ files
- `results/bam/` — aligned and sorted BAM files plus indices
- `results/peaks/` — MACS2 broad peak BED files and consensus peaks
- `results/bigwig/` — normalized signal tracks
- `results/deeptools/` — heatmap/profile plots
- `results/diffbind/` — differential binding summary and plots
- `results/motifs/` — HOMER motif results and summaries
- `results/rnaseq/` — Salmon quantification and DESeq2 results
- `results/integrative/` — integrated peak/gene analysis and scatter plots

## Notes and best practices
- This workflow is designed for paired-end ChIP and RNA datasets.
- Use FASTQ and reference files that match the selected species assembly.
- Confirm `bt2_index` points to a built Bowtie2 index set.
- For genome-specific peak metrics, adjust `MACS2` parameters in `Snakefile`.
- If you use `mamba`, replace `conda env create` with `mamba env create`.

## GitHub repository
This repository has been initialized locally and pushed to GitHub using the authenticated `gh` CLI account.

## Support
If you want, I can also add:
- a `report` rule for a combined HTML analysis report
- `config.sample.yaml` and a sample data manifest
- species-specific reference download helpers
