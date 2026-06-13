# End-to-end audit

Audit date: 2026-06-02

> **Update 2026-06-13:** HOMER motif enrichment was replaced with an open-source `monaLisa` + `JASPAR`
> step (`analysis/motif_enrichment.R`, `motif_enrichment` rule). The pipeline now has no academic-only
> or non-redistributable dependencies, and the repository was made public.

## Scope

This audit covers the repository state after refocusing the workflow on histone ChIP-seq differential binding. It checks repository hygiene, license, repository visibility, ChIP-only scope, workflow graph consistency, contamination screening, MultiQC aggregation, configuration schema alignment, helper scripts, environments, and expected runtime inputs.

## Findings and status

| Area | Status | Notes |
| --- | --- | --- |
| Visibility | Pass | Repository is public; all dependencies are open-source and freely redistributable. |
| License | Pass | Repository includes an MIT license. |
| ChIP-only scope | Pass | RNA-seq quantification, differential expression, and transcriptomic integration files are not part of the workflow. |
| Workflow graph | Pass | `rule all` targets ChIP-seq QC, contamination screening, alignment, peak calling, blacklist filtering, bigWig tracks, deepTools plots, DiffBind outputs, motif outputs, and MultiQC. |
| Contamination detection | Pass | FastQ Screen runs on every raw paired-end ChIP and input/control FASTQ using the configured reference database. |
| MultiQC aggregation | Pass | MultiQC aggregates FastQC, Trim Galore, FastQ Screen, alignment, peak-calling, DiffBind, motif, and workflow logs from `results/`. |
| Differential binding design | Pass | Config and manifest examples include two conditions with two ChIP-seq replicates per condition. Helper scripts enforce this design before DiffBind runs. |
| Sample sheet | Pass | `scripts/build_sample_sheets.py` creates a DiffBind sample sheet with BAMs, controls, filtered peaks, factor, tissue, condition, replicate, and peak caller fields. |
| References | Pass | Human, mouse, and rat reference names are internally consistent with the download helper aliases. |
| Environments | Pass | ChIP-seq and R analysis environments include the tools needed for ChIP-seq processing, FastQ Screen, MultiQC, DiffBind, and monaLisa/JASPAR motif enrichment. |
| Generated data hygiene | Pass | Results, raw FASTQs, references, BAMs, bigWigs, Snakemake metadata, and local FastQ Screen configs are ignored by Git. |

## Validation commands

Run these commands from a checkout after installing the declared environments and preparing a real FastQ Screen config:

```bash
python3 -m py_compile scripts/build_sample_sheets.py scripts/manifest_to_config.py scripts/download_references.py
python3 scripts/build_sample_sheets.py --config config.yaml --diffbind /tmp/diffbind_sample_sheet.csv
python3 scripts/manifest_to_config.py sample_manifest.tsv \
  --species human \
  --genome /path/to/hg38.fa \
  --chrom-sizes /path/to/hg38.chrom.sizes \
  --bt2-index /path/to/hg38 \
  --black-list /path/to/hg38-blacklist.v2.bed \
  --output /tmp/config.from_manifest.yaml
snakemake --use-conda --cores 12 results/multiqc/multiqc_report.html --dry-run
```

## Runtime prerequisites

- Replace all `/path/to/...` reference placeholders with real files before running.
- Copy `config/fastq_screen.conf.example` to `config/fastq_screen.conf`, replace all database index paths, and update `contamination.fastq_screen_conf` in `config.yaml` if needed.
- Provide paired-end FASTQ files for every `chip_samples` and `chip_controls` entry.
- Ensure the blacklist BED, genome FASTA, chromosome sizes file, Bowtie2 index, and FastQ Screen database indexes all use the intended genome assemblies.
- Use at least two biological ChIP-seq replicates per condition for DiffBind contrasts.
