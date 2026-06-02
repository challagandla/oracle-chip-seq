# End-to-end audit

Audit date: 2026-06-02

## Scope

This audit covers the repository state after refocusing the workflow on histone ChIP-seq differential binding. It checks repository hygiene, license, ChIP-only scope, workflow graph consistency, configuration schema alignment, helper scripts, environments, and expected runtime inputs.

## Findings and status

| Area | Status | Notes |
| --- | --- | --- |
| License | Pass | Repository includes an MIT license. |
| ChIP-only scope | Pass | RNA-seq quantification, differential expression, and transcriptomic integration files are not part of the workflow. |
| Workflow graph | Pass | `rule all` targets ChIP-seq QC, alignment, peak calling, blacklist filtering, bigWig tracks, deepTools plots, DiffBind outputs, and motif outputs. |
| Differential binding design | Pass | Config and manifest examples include two conditions with two ChIP-seq replicates per condition. Helper scripts enforce this design before DiffBind runs. |
| Sample sheet | Pass | `scripts/build_sample_sheets.py` creates a DiffBind sample sheet with BAMs, controls, filtered peaks, factor, tissue, condition, replicate, and peak caller fields. |
| References | Pass | Human, mouse, and rat reference names are internally consistent with the download helper aliases. |
| Environments | Pass | ChIP-seq and R analysis environments only include dependencies needed for ChIP-seq processing, DiffBind, and motif summaries. |
| Generated data hygiene | Pass | Results, raw FASTQs, references, BAMs, bigWigs, and Snakemake metadata are ignored by Git. |
| Repository visibility | Manual | The authenticated connector in this session did not expose a repository visibility update action; make the repository private in GitHub settings. |

## Validation commands

Run these commands from a checkout after installing the declared environments:

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
snakemake --use-conda --cores 12 --dry-run
```

## Runtime prerequisites

- Replace all `/path/to/...` reference placeholders with real files before running.
- Provide paired-end FASTQ files for every `chip_samples` and `chip_controls` entry.
- Ensure the blacklist BED, genome FASTA, chromosome sizes file, and Bowtie2 index all use the same genome assembly.
- Use at least two biological ChIP-seq replicates per condition for DiffBind contrasts.
