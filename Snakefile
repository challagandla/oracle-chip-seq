import os

configfile: "config.yaml"

SPECIES = config["species"]
REF = config["references"][SPECIES]
CONTAMINATION = config.get("contamination", {})
FASTQ_SCREEN_CONF = CONTAMINATION.get("fastq_screen_conf", "config/fastq_screen.conf.example")
FASTQ_SCREEN_SUBSET = CONTAMINATION.get("subset", 100000)

CHIP_SAMPLES = [s["id"] for s in config["chip_samples"]]
CHIP_CONTROLS = {c["id"]: c for c in config["chip_controls"]}
ALL_SAMPLES = CHIP_SAMPLES + list(CHIP_CONTROLS.keys())

FASTQ = {s["id"]: s["fastq"] for s in config["chip_samples"] + config["chip_controls"]}
CONTROL_MAP = {s["id"]: s["control"] for s in config["chip_samples"]}

RAW_FASTQC_HTML = expand("results/fastqc/raw/{sample}_{read}_fastqc.html", sample=ALL_SAMPLES, read=["R1", "R2"])
RAW_FASTQC_ZIP = expand("results/fastqc/raw/{sample}_{read}_fastqc.zip", sample=ALL_SAMPLES, read=["R1", "R2"])
TRIMMED_FASTQC_HTML = expand("results/fastqc/trimmed/{sample}_R1_val_1_fastqc.html", sample=ALL_SAMPLES) + expand("results/fastqc/trimmed/{sample}_R2_val_2_fastqc.html", sample=ALL_SAMPLES)
TRIMMED_FASTQC_ZIP = expand("results/fastqc/trimmed/{sample}_R1_val_1_fastqc.zip", sample=ALL_SAMPLES) + expand("results/fastqc/trimmed/{sample}_R2_val_2_fastqc.zip", sample=ALL_SAMPLES)
FASTQ_SCREEN_TEXT = expand("results/contamination/fastq_screen/{sample}_{read}_screen.txt", sample=ALL_SAMPLES, read=["R1", "R2"])
FASTQ_SCREEN_HTML = expand("results/contamination/fastq_screen/{sample}_{read}_screen.html", sample=ALL_SAMPLES, read=["R1", "R2"])
PEAKS = expand("results/peaks/{sample}_peaks.broadPeak", sample=CHIP_SAMPLES)
BIGWIGS = expand("results/bigwig/{sample}.bw", sample=CHIP_SAMPLES)
MULTIQC_REPORT = "results/multiqc/multiqc_report.html"


rule all:
    input:
        RAW_FASTQC_HTML,
        TRIMMED_FASTQC_HTML,
        FASTQ_SCREEN_TEXT,
        FASTQ_SCREEN_HTML,
        expand("results/bam/{sample}.sorted.bam.bai", sample=ALL_SAMPLES),
        PEAKS,
        BIGWIGS,
        "results/peaks/consensus_peaks.bed",
        "results/deeptools/matrix.gz",
        "results/deeptools/heatmap.png",
        "results/deeptools/profile.png",
        "results/diffbind/diffbind_summary.csv",
        "results/diffbind/diffbind_plots.pdf",
        "results/diffbind/diffbind.rds",
        "results/motifs/motif_enrichment.tsv",
        "results/motifs/motif_summary.pdf",
        MULTIQC_REPORT


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
        mkdir -p results/fastqc/raw results/logs
        fastqc -o results/fastqc/raw -t {threads} {input.R1} {input.R2} > {log} 2>&1
        """


rule fastq_screen_raw:
    input:
        R1=lambda wc: FASTQ[wc.sample][0],
        R2=lambda wc: FASTQ[wc.sample][1],
        conf=FASTQ_SCREEN_CONF
    output:
        R1txt="results/contamination/fastq_screen/{sample}_R1_screen.txt",
        R1html="results/contamination/fastq_screen/{sample}_R1_screen.html",
        R2txt="results/contamination/fastq_screen/{sample}_R2_screen.txt",
        R2html="results/contamination/fastq_screen/{sample}_R2_screen.html"
    params:
        outdir="results/contamination/fastq_screen",
        subset=FASTQ_SCREEN_SUBSET
    threads: 4
    log:
        "results/logs/fastq_screen_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        tmpdir=$(mktemp -d {params.outdir}/tmp.{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        ln -sf "$(readlink -f {input.R1})" "$tmpdir/{wildcards.sample}_R1.fastq.gz"
        ln -sf "$(readlink -f {input.R2})" "$tmpdir/{wildcards.sample}_R2.fastq.gz"
        fastq_screen --conf {input.conf} --aligner bowtie2 --threads {threads} --subset {params.subset} --outdir {params.outdir} "$tmpdir/{wildcards.sample}_R1.fastq.gz" "$tmpdir/{wildcards.sample}_R2.fastq.gz" > {log} 2>&1
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
        mkdir -p {params.outdir} results/logs
        trim_galore --paired --cores {threads} --output_dir {params.outdir} {input.R1} {input.R2} > {log} 2>&1
        """


rule fastqc_trimmed:
    input:
        R1="results/trimmed/{sample}_R1_val_1.fq.gz",
        R2="results/trimmed/{sample}_R2_val_2.fq.gz"
    output:
        "results/fastqc/trimmed/{sample}_R1_val_1_fastqc.html",
        "results/fastqc/trimmed/{sample}_R1_val_1_fastqc.zip",
        "results/fastqc/trimmed/{sample}_R2_val_2_fastqc.html",
        "results/fastqc/trimmed/{sample}_R2_val_2_fastqc.zip"
    threads: 2
    log:
        "results/logs/fastqc_trimmed_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/fastqc/trimmed results/logs
        fastqc -o results/fastqc/trimmed -t {threads} {input.R1} {input.R2} > {log} 2>&1
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
        mkdir -p results/bam results/logs
        (bowtie2 -x {params.index} -1 {input.trimmed1} -2 {input.trimmed2} -p {threads} 2> {log} | samtools view -bS - > {output}) 2>> {log}
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
        mkdir -p results/bam results/logs
        samtools sort -n -@ {threads} -o results/bam/{wildcards.sample}.name_sorted.bam {input.bam} > {log} 2>&1
        samtools fixmate -m results/bam/{wildcards.sample}.name_sorted.bam results/bam/{wildcards.sample}.fixmate.bam >> {log} 2>&1
        samtools sort -@ {threads} -o results/bam/{wildcards.sample}.coord_sorted.bam results/bam/{wildcards.sample}.fixmate.bam >> {log} 2>&1
        samtools markdup -r results/bam/{wildcards.sample}.coord_sorted.bam {output.bam} >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        rm -f results/bam/{wildcards.sample}.name_sorted.bam results/bam/{wildcards.sample}.fixmate.bam results/bam/{wildcards.sample}.coord_sorted.bam
        """


rule call_peaks:
    input:
        bam="results/bam/{sample}.sorted.bam",
        control=lambda wc: f"results/bam/{CONTROL_MAP[wc.sample]}.sorted.bam"
    output:
        broadpeak="results/peaks/raw/{sample}_peaks.broadPeak"
    params:
        gsize=REF["gsize"],
        outdir="results/peaks/raw"
    log:
        "results/logs/macs2_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        macs2 callpeak -t {input.bam} -c {input.control} --format BAMPE --name {wildcards.sample} --broad --broad-cutoff 0.1 --gsize {params.gsize} --outdir {params.outdir} > {log} 2>&1
        """


rule filter_blacklist:
    input:
        peaks="results/peaks/raw/{sample}_peaks.broadPeak",
        blacklist=lambda wc: REF["black_list"]
    output:
        "results/peaks/{sample}_peaks.broadPeak"
    log:
        "results/logs/filter_blacklist_{sample}.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/peaks results/logs
        bedtools intersect -v -a {input.peaks} -b {input.blacklist} > {output} 2> {log}
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
        mkdir -p results/peaks results/logs
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
        mkdir -p results/bigwig results/logs
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
        mkdir -p results/deeptools results/logs
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
        mkdir -p results/deeptools results/logs
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
        mkdir -p results/deeptools results/logs
        plotProfile -m {input.matrix} -out {output} --plotTitle "ChIP-seq meta-profile" > {log} 2>&1
        """


rule build_sample_sheets:
    input:
        config="config.yaml"
    output:
        diffbind="results/diffbind/sample_sheet.csv"
    log:
        "results/logs/build_sample_sheets.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/diffbind results/logs
        python3 scripts/build_sample_sheets.py --config {input.config} --diffbind {output.diffbind} > {log} 2>&1
        """


rule run_diffbind:
    input:
        sheet="results/diffbind/sample_sheet.csv",
        bams=expand("results/bam/{sample}.sorted.bam", sample=CHIP_SAMPLES + list(CHIP_CONTROLS.keys())),
        peaks=PEAKS
    output:
        summary="results/diffbind/diffbind_summary.csv",
        plots="results/diffbind/diffbind_plots.pdf",
        rds="results/diffbind/diffbind.rds"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.summary))
    log:
        "results/logs/diffbind.log"
    conda:
        "envs/r_analysis.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        Rscript analysis/diffbind_analysis.R {input.sheet} {params.outdir} > {log} 2>&1
        """


rule motif_enrichment:
    input:
        peaks="results/peaks/consensus_peaks.bed"
    output:
        table="results/motifs/motif_enrichment.tsv",
        summary="results/motifs/motif_summary.pdf"
    params:
        genome=REF["genome"],
        outdir=lambda wc, output: os.path.dirname(str(output.table))
    log:
        "results/logs/motif_enrichment.log"
    conda:
        "envs/r_analysis.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        Rscript analysis/motif_enrichment.R {input.peaks} {params.genome} {params.outdir} > {log} 2>&1
        """


rule multiqc:
    input:
        RAW_FASTQC_ZIP,
        TRIMMED_FASTQC_ZIP,
        FASTQ_SCREEN_TEXT,
        FASTQ_SCREEN_HTML,
        expand("results/bam/{sample}.sorted.bam.bai", sample=ALL_SAMPLES),
        PEAKS,
        BIGWIGS,
        "results/peaks/consensus_peaks.bed",
        "results/diffbind/diffbind_summary.csv",
        "results/diffbind/diffbind_plots.pdf",
        "results/motifs/motif_enrichment.tsv",
        "results/motifs/motif_summary.pdf"
    output:
        html=MULTIQC_REPORT
    log:
        "results/logs/multiqc.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/multiqc results/logs
        multiqc results --outdir results/multiqc --filename multiqc_report.html --force > {log} 2>&1
        """


rule snakemake_report:
    input:
        "results/peaks/consensus_peaks.bed",
        "results/diffbind/diffbind_summary.csv",
        "results/diffbind/diffbind_plots.pdf",
        "results/motifs/motif_summary.pdf",
        MULTIQC_REPORT
    output:
        "results/report/snakemake_report.html"
    threads: 1
    log:
        "results/report/snakemake_report.log"
    conda:
        "envs/chipseq.yaml"
    shell:
        """
        mkdir -p results/report
        snakemake --snakefile Snakefile --configfile config.yaml --use-conda --cores {threads} --report {output} --nolock > {log} 2>&1
        """
