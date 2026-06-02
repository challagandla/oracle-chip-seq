import os
from pathlib import Path

configfile: "config.yaml"

SPECIES = config["species"]
REF = config["references"][SPECIES]

CHIP_SAMPLES = [s["id"] for s in config["chip_samples"]]
CHIP_CONTROLS = {c["id"]: c for c in config["chip_controls"]}
RNA_SAMPLES = [s["id"] for s in config["rna_samples"]]
ALL_SAMPLES = CHIP_SAMPLES + list(CHIP_CONTROLS.keys()) + RNA_SAMPLES

FASTQ = {s["id"]: s["fastq"] for s in config["chip_samples"] + config["chip_controls"] + config["rna_samples"]}

CONTROL_MAP = {s["id"]: s["control"] for s in config["chip_samples"]}
PEAKS = expand("results/peaks/{sample}_peaks.broadPeak", sample=CHIP_SAMPLES)
BIGWIGS = expand("results/bigwig/{sample}.bw", sample=CHIP_SAMPLES)
SALMON_QUANTS = expand("results/rnaseq/salmon/{sample}/quant.sf", sample=RNA_SAMPLES)

rule all:
    input:
        expand("results/fastqc/raw/{sample}_R1_fastqc.html", sample=ALL_SAMPLES),
        expand("results/fastqc/trimmed/{sample}_R1_val_1_fastqc.html", sample=ALL_SAMPLES),
        expand("results/bam/{sample}.sorted.bam.bai", sample=ALL_SAMPLES),
        PEAKS,
        BIGWIGS,
        "results/peaks/consensus_peaks.bed",
        "results/deeptools/matrix.gz",
        "results/deeptools/heatmap.png",
        "results/deeptools/profile.png",
        "results/diffbind/diffbind_summary.csv",
        "results/motifs/homer/knownResults.txt",
        "results/motifs/motif_summary.pdf",
        "results/rnaseq/deseq2_results.tsv",
        "results/rnaseq/deseq2_plots.pdf",
        "results/integrative/integrative_summary.csv",
        "results/integrative/integrative_scatter.pdf"

rule fastqc_raw:
    input:
        R1=lambda wc: FASTQ[wc.sample][0],
        R2=lambda wc: FASTQ[wc.sample][1]
    output:
        "results/fastqc/raw/{sample}_R1_fastqc.html",
        "results/fastqc/raw/{sample}_R1_fastqc.zip",
        "results/fastqc/raw/{sample}_R2_fastqc.html",
        "results/fastqc/raw/{sample}_R2_fastqc.zip"
    threads: 2
    log:
        "results/logs/fastqc_raw_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        fastqc -o results/fastqc/raw -t {threads} {input.R1} {input.R2} > {log} 2>&1
        """

rule trim_galore:
    input:
        R1=lambda wc: FASTQ[wc.sample][0],
        R2=lambda wc: FASTQ[wc.sample][1]
    output:
        trimmed1="results/trimmed/{sample}_R1_val_1.fq.gz",
        trimmed2="results/trimmed/{sample}_R2_val_2.fq.gz"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.trimmed1))
    threads: 4
    log:
        "results/logs/trim_galore_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        trim_galore --paired --cores {threads} --output_dir {params.outdir} {input.R1} {input.R2} > {log} 2>&1
        """

rule align_bowtie2:
    input:
        trimmed1="results/trimmed/{sample}_R1_val_1.fq.gz",
        trimmed2="results/trimmed/{sample}_R2_val_2.fq.gz"
    output:
        temp("results/bam/{sample}.bam")
    threads: 8
    params:
        index=REF["bt2_index"]
    log:
        "results/logs/bowtie2_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        bowtie2 -x {params.index} -1 {input.trimmed1} -2 {input.trimmed2} -p {threads} 2> {log} | samtools view -bS - > {output}
        """

rule sort_markdup:
    input:
        bam="results/bam/{sample}.bam"
    output:
        bam="results/bam/{sample}.sorted.bam",
        bai="results/bam/{sample}.sorted.bam.bai"
    threads: 4
    log:
        "results/logs/sort_markdup_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        samtools sort -@ {threads} -o results/bam/{wildcards.sample}.sorted.raw.bam {input.bam} >> {log} 2>&1
        samtools fixmate -m results/bam/{wildcards.sample}.sorted.raw.bam results/bam/{wildcards.sample}.fixmate.bam >> {log} 2>&1
        samtools markdup -r results/bam/{wildcards.sample}.fixmate.bam {output.bam} >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        rm -f results/bam/{wildcards.sample}.sorted.raw.bam
        """

rule call_peaks:
    input:
        bam="results/bam/{sample}.sorted.bam",
        control=lambda wc: f"results/bam/{CONTROL_MAP[wc.sample]}.sorted.bam"
    output:
        broadpeak="results/peaks/{sample}_peaks.broadPeak"
    params:
        gsize=REF["gsize"]
    log:
        "results/logs/macs2_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        macs2 callpeak -t {input.bam} -c {input.control} --format BAMPE --name {wildcards.sample} --broad --broad-cutoff 0.1 --gsize {params.gsize} --outdir results/peaks > {log} 2>&1
        """

rule merge_peaks:
    input:
        PEAKS
    output:
        "results/peaks/consensus_peaks.bed"
    log:
        "results/logs/merge_peaks.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output} 2> {log}
        """

rule bamcoverage:
    input:
        bam="results/bam/{sample}.sorted.bam"
    output:
        bw="results/bigwig/{sample}.bw"
    threads: 4
    params:
        gsize=REF["gsize"]
    log:
        "results/logs/bamcoverage_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} --normalizeUsing RPGC --effectiveGenomeSize {params.gsize} --binSize 25 --extendReads 200 > {log} 2>&1
        """

rule deeptools_matrix:
    input:
        bigwigs=BIGWIGS,
        peaks="results/peaks/consensus_peaks.bed"
    output:
        "results/deeptools/matrix.gz"
    log:
        "results/logs/deeptools_matrix.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        computeMatrix scale-regions -S {input.bigwigs} -R {input.peaks} --regionBodyLength 10000 --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --skipZeros -o {output} > {log} 2>&1
        """

rule deeptools_heatmap:
    input:
        matrix="results/deeptools/matrix.gz"
    output:
        "results/deeptools/heatmap.png"
    log:
        "results/logs/deeptools_heatmap.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        plotHeatmap -m {input.matrix} -out {output} --plotTitle "ChIP-seq signal heatmap" --colorMap RdBu > {log} 2>&1
        """

rule deeptools_profile:
    input:
        matrix="results/deeptools/matrix.gz"
    output:
        "results/deeptools/profile.png"
    log:
        "results/logs/deeptools_profile.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        plotProfile -m {input.matrix} -out {output} --plotTitle "ChIP-seq meta-profile" > {log} 2>&1
        """

rule build_sample_sheets:
    input:
        config="config.yaml"
    output:
        diffbind="results/diffbind/sample_sheet.csv",
        rna="results/rnaseq/sample_metadata.tsv"
    log:
        "results/logs/build_sample_sheets.log"
    conda:
        "envs/chipseq.yaml"
    script:
        "scripts/build_sample_sheets.py"

rule homer_motif:
    input:
        peaks="results/peaks/consensus_peaks.bed"
    output:
        "results/motifs/homer/knownResults.txt"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output)),
        genome=REF["name"]
    log:
        "results/logs/homer_motif.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        findMotifsGenome.pl {input.peaks} {params.genome} {params.outdir} -size 200 -len 8,10,12 > {log} 2>&1
        """

rule run_diffbind:
    input:
        sheet="results/diffbind/sample_sheet.csv",
        bams=expand("results/bam/{sample}.sorted.bam", sample=CHIP_SAMPLES + list(CHIP_CONTROLS.keys())),
        peaks=PEAKS
    output:
        summary="results/diffbind/diffbind_summary.csv"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.summary))
    log:
        "results/logs/diffbind.log"
    conda:
        "envs/r_analysis.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript analysis/diffbind_analysis.R {input.sheet} {params.outdir} > {log} 2>&1
        """

rule motif_summary:
    input:
        homer="results/motifs/homer/knownResults.txt"
    output:
        "results/motifs/motif_summary.pdf"
    log:
        "results/logs/motif_summary.log"
    conda:
        "envs/r_analysis.yaml"
    shell:
        """
        Rscript analysis/motif_summary.R {input.homer} results/motifs > {log} 2>&1
        """

rule salmon_index:
    input:
        transcriptome=REF["transcriptome"]
    output:
        touch("results/rnaseq/salmon_index.done")
    threads: 4
    log:
        "results/logs/salmon_index.log"
    conda:
        "envs/rna_seq.yaml"
    shell:
        """
        salmon index -t {input.transcriptome} -i results/rnaseq/salmon_index --type quasi > {log} 2>&1
        touch {output}
        """

rule salmon_quant:
    input:
        index="results/rnaseq/salmon_index.done",
        fastq1=lambda wc: FASTQ[wc.sample][0],
        fastq2=lambda wc: FASTQ[wc.sample][1]
    output:
        "results/rnaseq/salmon/{sample}/quant.sf"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output)),
    threads: 8
    log:
        "results/logs/salmon_quant_{sample}.log"
    conda:
        "envs/rna_seq.yaml"
    shell:
        """
        salmon quant -i results/rnaseq/salmon_index -l A -1 {input.fastq1} -2 {input.fastq2} -p {threads} -o {params.outdir} > {log} 2>&1
        """

rule rnaseq_deseq2:
    input:
        metadata="results/rnaseq/sample_metadata.tsv",
        quant=SALMON_QUANTS
    output:
        results="results/rnaseq/deseq2_results.tsv",
        plot="results/rnaseq/deseq2_plots.pdf"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.results))
    log:
        "results/logs/rnaseq_deseq2.log"
    conda:
        "envs/r_analysis.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript analysis/rnaseq_DE.R {input.metadata} {params.outdir} {output.results} {output.plot} > {log} 2>&1
        """

rule integrative_analysis:
    input:
        diffbind="results/diffbind/diffbind_summary.csv",
        deseq="results/rnaseq/deseq2_results.tsv",
        annotation=REF["annotation"]
    output:
        summary="results/integrative/integrative_summary.csv",
        plot="results/integrative/integrative_scatter.pdf"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.summary))
    log:
        "results/logs/integrative_analysis.log"
    conda:
        "envs/r_analysis.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript analysis/integrative_analysis.R {input.diffbind} {input.deseq} {input.annotation} {params.outdir} > {log} 2>&1
        """

rule snakemake_report:
    input:
        "results/peaks/consensus_peaks.bed",
        "results/diffbind/diffbind_summary.csv",
        "results/rnaseq/deseq2_results.tsv"
    output:
        "results/report/snakemake_report.html"
    log:
        "results/report/snakemake_report.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/report
        snakemake --snakefile Snakefile --configfile config.yaml --report {output} --nolock > {log} 2>&1
        """
