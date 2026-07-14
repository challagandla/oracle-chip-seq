# Data-independent installer workflow for every unique ChIP-seq rule environment.
# setup.sh uses this file only to create and smoke-test packages; it never reads
# biological data or writes pipeline results.

rule all:
    input:
        ".snakemake/setup-env-checks/chipseq.ok",
        ".snakemake/setup-env-checks/r-analysis.ok",

rule check_chipseq_env:
    output:
        touch(".snakemake/setup-env-checks/chipseq.ok")
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        for tool in snakemake fastqc fastq_screen trim_galore bowtie2 samtools \
            bedtools macs3 bamCoverage bamCompare computeMatrix plotHeatmap plotProfile \
            plotFingerprint multiBigwigSummary plotCorrelation plotPCA \
            findMotifsGenome.pl multiqc wget gunzip gzip; do command -v "$tool" >/dev/null; done
        macs3 --version >/dev/null
        bamCoverage --version >/dev/null
        computeMatrix --version >/dev/null
        python -c 'import numpy, pandas, yaml; assert 2 <= int(numpy.__version__.split(".")[0]) < 3'
        touch {output}
        """

rule check_r_analysis_env:
    output:
        touch(".snakemake/setup-env-checks/r-analysis.ok")
    conda:
        "../envs/r_analysis.yaml"
    shell:
        """
        mkdir -p .snakemake/setup-env-checks
        Rscript --vanilla -e 'suppressPackageStartupMessages({{library(DiffBind); library(readr); library(dplyr); library(ggplot2)}})'
        touch {output}
        """
