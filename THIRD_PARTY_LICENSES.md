# Third-party tools, data & licenses

This pipeline's own code is **MIT** (see `LICENSE`). It **bundles no third-party source code and no
reference data** — tools are installed via conda/bioconda and invoked as separate processes, and
reference genomes/blacklists/indexes are provided/downloaded by the user. Invoking a separate program
is not a derivative work, so the MIT license is unaffected by the licenses of the tools it orchestrates.

## No commercial-use restrictions

Every tool this pipeline orchestrates is open-source and freely usable, including for commercial work,
with citation. There are **no academic-only, non-commercial, or non-redistributable dependencies**.
(Motif enrichment uses `monaLisa` + `JASPAR` rather than HOMER for exactly this reason.)

## Tools (installed via conda; invoked, not bundled)

| Tool | License | Role |
|---|---|---|
| Bowtie2 | GPL-3.0 | alignment |
| samtools / bedtools | MIT / GPL-2 | bam/interval ops |
| MACS2 | BSD-3-Clause | peak calling |
| deepTools | GPL-3.0 | coverage/QC |
| FastQC | GPL-3.0 | QC |
| FastQ Screen | GPL-3.0 | contamination QC |
| Trim Galore | GPL-3.0 | trimming |
| MultiQC | GPL-3.0 | report aggregation |

**R/Bioconductor:** DiffBind (Artistic-2.0), monaLisa (GPL-3.0), TFBSTools (GPL-2), JASPAR2020 motif
data (CC0 / public domain), GenomicRanges/rtracklayer/Rsamtools/Biostrings/SummarizedExperiment
(Artistic-2.0), ggplot2 (MIT), dplyr/readr (MIT). Invoked within R; not redistributed.

**On the GPL tools:** called as independent executables (Snakemake rules) — mere aggregation, not
linking — so they impose no copyleft obligation on this MIT pipeline.

## Reference data (provided/downloaded by the user; not redistributed)
Genome / Bowtie2 index / annotation from GENCODE/Ensembl (open); ENCODE blacklist from the Boyle Lab
(open, cite). None are bundled in this repository.

## Bottom line
No code-incorporation conflict, no bundled code/data, and no academic-only or non-redistributable
dependencies. Everything is freely usable (including commercially) with citation.
