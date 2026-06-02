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
        lambda wc: FASTQ[wc.sample]
    output:
        html=expand("results/fastqc/raw/{fq_basename}_fastqc.html", fq_basename=lambda wildcards: [Path(f).name.replace(".fastq.gz", "") for f in FASTQ[wildcards.sample]]),
        zip=expand("results/fastqc/raw/{fq_basename}_fastqc.zip", fq_basename=lambda wildcards: [Path(f).name.replace(".fastq.gz", "") for f in FASTQ[wildcards.sample]])
    threads: 2
    shell:
        """
        fastqc -o results/fastqc/raw -t {threads} {input}
        """

rule trim_galore:
    input:
        lambda wc: FASTQ[wc.sample]
    output:
        trimmed1="results/trimmed/{sample}_R1_val_1.fq.gz",
        trimmed2="results/trimmed/{sample}_R2_val_2.fq.gz"
    params:
        outdir="results/trimmed"
    threads: 4
    shell:
        """
        trim_galore --paired --cores {threads} --output_dir {params.outdir} {input[0]} {input[1]}
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
    shell:
        """
        bowtie2 -x {params.index} -1 {input.trimmed1} -2 {input.trimmed2} -p {threads} 2> results/bam/{wildcards.sample}.bowtie2.log | samtools view -bS - > {output}
        """

rule sort_markdup:
    input:
        bam="results/bam/{sample}.bam"
    output:
        bam="results/bam/{sample}.sorted.bam",
        bai="results/bam/{sample}.sorted.bam.bai"
    threads: 4
    shell:
        """
        samtools sort -@ {threads} -o results/bam/{wildcards.sample}.sorted.raw.bam {input.bam}
        samtools fixmate -m results/bam/{wildcards.sample}.sorted.raw.bam results/bam/{wildcards.sample}.fixmate.bam
        samtools markdup -r results/bam/{wildcards.sample}.fixmate.bam {output.bam}
        samtools index {output.bam}
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
    shell:
        """
        macs2 callpeak -t {input.bam} -c {input.control} --format BAMPE --name {wildcards.sample} --broad --broad-cutoff 0.1 --gsize {params.gsize} --outdir results/peaks
        """

rule merge_peaks:
    input:
        PEAKS
    output:
        "results/peaks/consensus_peaks.bed"
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bedtools merge > {output}
        """

rule bamcoverage:
    input:
        bam="results/bam/{sample}.sorted.bam"
    output:
        bw="results/bigwig/{sample}.bw"
    threads: 4
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} --normalizeUsing RPGC --effectiveGenomeSize {REF["gsize"]} --binSize 25 --extendReads 200
        """

rule deeptools_matrix:
    input:
        bigwigs=BIGWIGS,
        peaks="results/peaks/consensus_peaks.bed"
    output:
        "results/deeptools/matrix.gz"
    shell:
        """
        computeMatrix scale-regions -S {input.bigwigs} -R {input.peaks} --regionBodyLength 10000 --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --skipZeros -o {output}
        """

rule deeptools_heatmap:
    input:
        matrix="results/deeptools/matrix.gz"
    output:
        "results/deeptools/heatmap.png"
    shell:
        """
        plotHeatmap -m {input.matrix} -out {output} --plotTitle "ChIP-seq signal heatmap" --colorMap RdBu
        """

rule deeptools_profile:
    input:
        matrix="results/deeptools/matrix.gz"
    output:
        "results/deeptools/profile.png"
    shell:
        """
        plotProfile -m {input.matrix} -out {output} --plotTitle "ChIP-seq meta-profile"
        """

rule build_sample_sheets:
    input:
        config="config.yaml"
    output:
        diffbind="results/diffbind/sample_sheet.csv",
        rna="results/rnaseq/sample_metadata.tsv"
    script:
        "scripts/build_sample_sheets.py"

rule homer_motif:
    input:
        peaks="results/peaks/consensus_peaks.bed"
    output:
        "results/motifs/homer/knownResults.txt"
    params:
        outdir="results/motifs/homer",
        genome=REF["name"]
    shell:
        """
        mkdir -p {params.outdir}
        findMotifsGenome.pl {input.peaks} {params.genome} {params.outdir} -size 200 -len 8,10,12
        """

rule run_diffbind:
    input:
        sheet="results/diffbind/sample_sheet.csv",
        bams=expand("results/bam/{sample}.sorted.bam", sample=CHIP_SAMPLES + list(CHIP_CONTROLS.keys())),
        peaks=PEAKS
    output:
        summary="results/diffbind/diffbind_summary.csv"
    params:
        outdir="results/diffbind"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript analysis/diffbind_analysis.R {input.sheet} {params.outdir}
        """

rule motif_summary:
    input:
        homer="results/motifs/homer/knownResults.txt"
    output:
        "results/motifs/motif_summary.pdf"
    shell:
        """
        Rscript analysis/motif_summary.R {input.homer} results/motifs
        """

rule salmon_index:
    input:
        transcriptome=REF["transcriptome"]
    output:
        touch("results/rnaseq/salmon_index.done")
    threads: 4
    shell:
        """
        salmon index -t {input.transcriptome} -i results/rnaseq/salmon_index --type quasi
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
        outdir="results/rnaseq/salmon/{sample}"
    threads: 8
    shell:
        """
        salmon quant -i results/rnaseq/salmon_index -l A -1 {input.fastq1} -2 {input.fastq2} -p {threads} -o {params.outdir}
        """

rule rnaseq_deseq2:
    input:
        metadata="results/rnaseq/sample_metadata.tsv",
        quant=SALMON_QUANTS
    output:
        results="results/rnaseq/deseq2_results.tsv",
        plot="results/rnaseq/deseq2_plots.pdf"
    params:
        outdir="results/rnaseq"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript analysis/rnaseq_DE.R {input.metadata} {params.outdir} {output.results} {output.plot}
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
        outdir="results/integrative"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript analysis/integrative_analysis.R {input.diffbind} {input.deseq} {input.annotation} {params.outdir}
        """
