# Oracle Histone ChIP-seq Differential Binding Pipeline

A reproducible `Snakemake` workflow for histone ChIP-seq quality control, paired-end alignment, duplicate removal, broad peak calling, blacklist filtering, `deepTools` visualization, `DiffBind` differential binding analysis, and HOMER motif enrichment across human, mouse, and rat.

## Repository structure
- `Snakefile` - main Snakemake workflow
- `config.yaml` - species references, ChIP-seq sample definitions, and design metadata
- `config.sample.yaml` - example config template with a valid two-condition, two-replicate DiffBind design
- `sample_manifest.tsv` - example manifest for config generation
- `envs/` - Conda environment definitions for ChIP-seq processing and R analysis
- `scripts/` - helpers for building DiffBind sample sheets, reference downloads, and manifest-driven configs
- `analysis/` - R analysis scripts for DiffBind and motif summaries
- `AUDIT.md` - end-to-end audit checklist and validation notes
- `.gitignore` - files and folders excluded from Git tracking
- `LICENSE` - MIT license

## What this pipeline does
1. Raw FASTQ QC with `FastQC`
2. Adapter trimming with `Trim Galore`
3. Paired-end ChIP-seq alignment with `Bowtie2`
4. BAM sorting, indexing, and duplicate removal with `samtools`
5. Broad peak calling for histone marks with `MACS2`
6. Blacklist filtering with `bedtools`
7. BigWig coverage track generation with `deepTools`
8. Heatmap/profile generation with `computeMatrix`, `plotHeatmap`, and `plotProfile`
9. Differential binding analysis with `DiffBind`
10. Motif enrichment analysis with `HOMER`

## Scope
This repository is intentionally ChIP-seq only. It does not run RNA-seq quantification, differential expression, or transcriptomics integration. The statistical comparison step is differential binding analysis from ChIP-seq peaks and aligned reads.

## Quick start
### 1. Fill in `config.yaml`
- Set `species` to `human`, `mouse`, or `rat`.
- Provide absolute paths for:
  - `genome`
  - `chrom_sizes`
  - `black_list`
  - `bt2_index`
- Update ChIP sample names, FASTQ paths, input controls, conditions, replicates, histone mark/factor, and tissue labels.
- Include at least two biological replicates per condition for reliable differential binding analysis.

### 2. Place raw data
Put FASTQ pairs in `data/raw/` or update paths in `config.yaml`.

### 3. Create Conda environments
```bash
conda env create -f envs/chipseq.yaml
conda env create -f envs/r_analysis.yaml
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
python3 scripts/build_sample_sheets.py --config config.yaml --diffbind results/diffbind/sample_sheet.csv
Rscript analysis/diffbind_analysis.R results/diffbind/sample_sheet.csv results/diffbind
Rscript analysis/motif_summary.R results/motifs/homer/knownResults.txt results/motifs
```

## How to use the pipeline
- `snakemake --use-conda --cores 12` builds all core outputs in `results/`.
- `snakemake results/diffbind/sample_sheet.csv` generates the DiffBind sample sheet.
- `snakemake results/peaks/consensus_peaks.bed` creates a merged peak set for deepTools.
- `snakemake results/diffbind/diffbind_summary.csv` runs differential binding analysis.
- `snakemake results/report/snakemake_report.html` generates a Snakemake HTML report of workflow execution.

## Sample manifest and config generation
1. Edit `sample_manifest.tsv` with your sample names and file paths.
2. Generate `config.yaml` from the manifest:
```bash
python3 scripts/manifest_to_config.py sample_manifest.tsv \
  --species human \
  --genome /path/to/hg38.fa \
  --chrom-sizes /path/to/hg38.chrom.sizes \
  --bt2-index /path/to/hg38 \
  --black-list /path/to/hg38-blacklist.v2.bed \
  --output config.yaml
```

## Reference download helper
Use the species-specific helper to download a genome FASTA, build a Bowtie2 index, and generate chromosome sizes:
```bash
python3 scripts/download_references.py human --outdir references
```

## Generate a combined Snakemake report
```bash
snakemake --use-conda --cores 12 results/report/snakemake_report.html
```

## Expected outputs
- `results/fastqc/` - QC reports
- `results/trimmed/` - trimmed FASTQ files
- `results/bam/` - aligned and sorted BAM files plus indices
- `results/peaks/raw/` - raw MACS2 broadPeak files
- `results/peaks/` - blacklist-filtered broadPeak files and consensus peaks
- `results/bigwig/` - normalized signal tracks
- `results/deeptools/` - heatmap/profile plots
- `results/diffbind/` - differential binding summary, plots, and serialized DiffBind object
- `results/motifs/` - HOMER motif results and summaries
- `results/report/` - optional Snakemake HTML report

## Notes and best practices
- Use FASTQ and reference files that match the selected species assembly.
- Confirm `bt2_index` points to a built Bowtie2 index set.
- Confirm `black_list` points to the matching genome blacklist file.
- Check the DiffBind sample sheet before running contrasts; each condition should have biological replicates.
- For genome-specific peak metrics, adjust `MACS2` parameters in `Snakefile`.
- If you use `mamba`, replace `conda env create` with `mamba env create`.
