# Third-party tools, data & licenses

This pipeline's own code is **MIT** (see `LICENSE`). It **bundles no third-party source code and no
reference data** — tools are installed via conda/bioconda and invoked as separate processes, and
reference genomes/blacklists/indexes are provided/downloaded by the user. Invoking a separate program
is not a derivative work, so the MIT license is unaffected by the licenses of the tools it orchestrates.

## ⚠️ Restriction that affects commercial use

| Dependency | Used by | Terms |
|---|---|---|
| **HOMER** | motif enrichment (`homer_motif` rule in `Snakefile`; `envs/chipseq.yaml`) | HOMER is **freeware for academic / non-profit use, is not open-source, and may not be redistributed**; commercial use requires contacting the author (C. Benner, Salk/UCSD). HOMER genome packages (`configureHomer.pl`) are derived from **UCSC** data (academic/non-profit; commercial needs a UCSC license). |

HOMER is conda-installed and invoked (not bundled). The **motif step is academic/non-commercial**; for
commercial use, obtain HOMER/UCSC permissions or skip the HOMER motif rule (alignment, peak calling,
coverage, QC and DiffBind do not depend on it).

## Tools (installed via conda; invoked, not bundled)

| Tool | License | Role |
|---|---|---|
| Bowtie2 | GPL-3.0 | alignment |
| samtools / bedtools | MIT / GPL-2 | bam/interval ops |
| MACS2 | BSD-3-Clause | peak calling |
| deepTools | GPL-3.0 | coverage/QC |
| **HOMER** | **academic/non-profit; not redistributable** (see above) | motif enrichment |
| FastQC | GPL-3.0 | QC |
| FastQ Screen | GPL-3.0 | contamination QC |
| Trim Galore | GPL-3.0 | trimming |
| MultiQC | GPL-3.0 | report aggregation |

**R/Bioconductor:** DiffBind (Artistic-2.0), ggplot2 (MIT), dplyr/readr (MIT). Invoked within R; not
redistributed.

**On the GPL tools:** called as independent executables (Snakemake rules) — mere aggregation, not
linking — so they impose no copyleft obligation on this MIT pipeline.

## Reference data (provided/downloaded by the user; not redistributed)
Genome / Bowtie2 index / annotation from GENCODE/Ensembl (open); ENCODE blacklist from the Boyle Lab
(open, cite). None are bundled in this repository.

## Bottom line
No code-incorporation conflict, no bundled code/data. The one operative constraint is **HOMER**
(academic/non-profit; not redistributable; commercial use needs author/UCSC permission). Everything
else is freely usable with citation.
