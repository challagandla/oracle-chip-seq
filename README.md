# Oracle Histone ChIP-seq Differential Binding Pipeline

A reproducible `Snakemake` workflow for histone ChIP-seq quality control, contamination screening with `FastQ Screen`, paired-end alignment, duplicate removal, broad peak calling, blacklist filtering, `deepTools` visualization, `DiffBind` differential binding analysis, open-source motif enrichment with `monaLisa` + `JASPAR`, and `MultiQC` reporting across human, mouse, and rat.

## Repository structure
- `Snakefile` - main Snakemake workflow
- `config.yaml` - species references, ChIP-seq sample definitions, contamination settings, and design metadata
- `config.sample.yaml` - example config template with a valid two-condition, two-replicate DiffBind design
- `config/fastq_screen.conf.example` - template FastQ Screen contamination database config
- `sample_manifest.tsv` - example manifest for config generation
- `envs/` - Conda environment definitions for ChIP-seq processing and R analysis
- `scripts/` - helpers for building DiffBind sample sheets, reference downloads, and manifest-driven configs
- `analysis/` - R analysis scripts for DiffBind and motif enrichment
- `AUDIT.md` - end-to-end audit checklist and validation notes
- `.gitignore` - files and folders excluded from Git tracking
- `LICENSE` - MIT license

## What this pipeline does
1. Raw FASTQ QC with `FastQC`
2. Raw FASTQ contamination screening with `FastQ Screen`
3. Adapter trimming with `Trim Galore`
4. Trimmed FASTQ QC with `FastQC`
5. Paired-end ChIP-seq alignment with `Bowtie2`
6. BAM sorting, indexing, and duplicate removal with `samtools`
7. Broad peak calling for histone marks with `MACS2`
8. Blacklist filtering with `bedtools`
9. BigWig coverage track generation with `deepTools`
10. Heatmap/profile generation with `computeMatrix`, `plotHeatmap`, and `plotProfile`
11. Differential binding analysis with `DiffBind`
12. Motif enrichment analysis with `monaLisa` against `JASPAR` vertebrate core motifs (open-source; no species-specific background required)
13. Project-wide reporting with `MultiQC`

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
- Configure `contamination.fastq_screen_conf` to point to a FastQ Screen config file with Bowtie2 index prefixes for the genomes or contaminants you want to screen.
- Update ChIP sample names, FASTQ paths, input controls, conditions, replicates, histone mark/factor, and tissue labels.
- Include at least two biological replicates per condition for reliable differential binding analysis.

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

### 4. Create Conda environments
```bash
conda env create -f envs/chipseq.yaml
conda env create -f envs/r_analysis.yaml
```

### 5. Run the full workflow
```bash
conda activate chipseq
snakemake --use-conda --cores 12
```

If you already have Snakemake installed in a separate environment, for example:
```bash
/home/epigenetics/miniforge3/envs/snakemake/bin/snakemake --use-conda --cores 12
```

### 6. Run only R analysis steps
```bash
conda activate r_analysis
python3 scripts/build_sample_sheets.py --config config.yaml --diffbind results/diffbind/sample_sheet.csv
Rscript analysis/diffbind_analysis.R results/diffbind/sample_sheet.csv results/diffbind
Rscript analysis/motif_enrichment.R results/peaks/consensus_peaks.bed /path/to/genome.fa results/motifs
```

## How to use the pipeline
- `snakemake --use-conda --cores 12` builds all core outputs in `results/`.
- `snakemake results/contamination/fastq_screen/H3K27ac_control_rep1_R1_screen.txt` runs one FastQ Screen contamination check.
- `snakemake results/diffbind/sample_sheet.csv` generates the DiffBind sample sheet.
- `snakemake results/peaks/consensus_peaks.bed` creates a merged peak set for deepTools.
- `snakemake results/diffbind/diffbind_summary.csv` runs differential binding analysis.
- `snakemake results/multiqc/multiqc_report.html` generates the MultiQC report.
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

Build separate Bowtie2 indexes for any extra contamination-screening databases you add to `config/fastq_screen.conf`, such as PhiX, E. coli, yeast, or alternative host genomes.

## Generate reports
```bash
snakemake --use-conda --cores 12 results/multiqc/multiqc_report.html
snakemake --use-conda --cores 12 results/report/snakemake_report.html
```

## Expected outputs
- `results/fastqc/` - raw and trimmed FASTQ QC reports
- `results/contamination/fastq_screen/` - raw FASTQ contamination-screen reports
- `results/trimmed/` - trimmed FASTQ files
- `results/bam/` - aligned and sorted BAM files plus indices
- `results/peaks/raw/` - raw MACS2 broadPeak files
- `results/peaks/` - blacklist-filtered broadPeak files and consensus peaks
- `results/bigwig/` - normalized signal tracks
- `results/deeptools/` - heatmap/profile plots
- `results/diffbind/` - differential binding summary, plots, and serialized DiffBind object
- `results/motifs/` - motif enrichment table (`motif_enrichment.tsv`) and top-motif summary (`motif_summary.csv`/`.pdf`)
- `results/multiqc/` - consolidated MultiQC report and parsed metrics
- `results/report/` - optional Snakemake HTML report

## Notes and best practices
- Use FASTQ and reference files that match the selected species assembly.
- Confirm `bt2_index` points to a built Bowtie2 index set.
- Confirm `black_list` points to the matching genome blacklist file.
- Confirm FastQ Screen database entries point to valid Bowtie2 index prefixes.
- Check the DiffBind sample sheet before running contrasts; each condition should have biological replicates.
- For genome-specific peak metrics, adjust `MACS2` parameters in `Snakefile`.
- If you use `mamba`, replace `conda env create` with `mamba env create`.

## License & usage

The pipeline's own code is **MIT** (see [LICENSE](LICENSE)). It bundles no third-party code or data;
tools are conda-installed and invoked, so the MIT license is unaffected by the (incl. GPL) licenses of
those tools. Every orchestrated tool is open-source and freely usable (including commercially) with
citation — there are no academic-only or non-redistributable dependencies. Full breakdown:
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).
