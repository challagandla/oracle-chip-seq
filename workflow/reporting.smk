rule multiqc:
    input:
        RAW_FASTQC_ZIP,
        TRIMMED_FASTQC_ZIP,
        FASTQ_SCREEN_TEXT,
        FASTQ_SCREEN_HTML,
        expand("results/bam/{sample}.sorted.bam.bai", sample=ALL_SAMPLES),
        expand("results/qc/samtools/{sample}.flagstat.txt", sample=ALL_SAMPLES),
        expand("results/qc/samtools/{sample}.markdup.txt", sample=ALL_SAMPLES),
        PEAKS,
        RPGC_BIGWIGS,
        LOG2_BIGWIGS,
        FACTOR_CONSENSUS,
        "results/peaks/consensus_peaks.bed",
        DIFFBIND_OUTPUTS,
        MOTIF_OUTPUTS,
        DEEPTOOLS_OUTPUTS
    output:
        html=MULTIQC_REPORT
    log:
        "results/logs/multiqc.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/multiqc results/logs
        multiqc results --outdir results/multiqc --filename multiqc_report.html --force > {log} 2>&1
        """


rule snakemake_report:
    input:
        "results/peaks/consensus_peaks.bed",
        DIFFBIND_OUTPUTS,
        MOTIF_OUTPUTS,
        DEEPTOOLS_OUTPUTS,
        MULTIQC_REPORT,
        config="results/provenance/resolved_config.yaml"
    output:
        "results/report/snakemake_report.html"
    threads: 1
    log:
        "results/report/snakemake_report.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/report
        snakemake --snakefile Snakefile --configfile {input.config:q} --use-conda \
            --cores {threads} --report {output:q} --nolock > {log:q} 2>&1
        """
