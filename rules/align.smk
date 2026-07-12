# Trimming, alignment and read filtering.
#
# Library layout is resolved per sample, not per run: this dataset mixes
# single-end (replicate 1) and paired-end (replicate 2) libraries, which is
# common when a study is topped up after the fact. Every rule below branches on
# is_paired() rather than assuming one layout.


rule trim_pe:
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        report="results/trimmed/{sample}_trimming_report.txt",
    params:
        outdir="results/trimmed",
    threads: 4
    log:
        "results/logs/trim/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        mkdir -p {params.outdir} $(dirname {log})
        trim_galore --paired --cores {threads} --gzip \
            --output_dir {params.outdir} {input.r1} {input.r2} > {log} 2>&1
        mv {params.outdir}/{wildcards.sample}_R1_val_1.fq.gz {output.r1}
        mv {params.outdir}/{wildcards.sample}_R2_val_2.fq.gz {output.r2}
        cat {params.outdir}/{wildcards.sample}_R1.fastq.gz_trimming_report.txt \
            {params.outdir}/{wildcards.sample}_R2.fastq.gz_trimming_report.txt > {output.report}
        """


rule trim_se:
    input:
        r1="data/raw/{sample}.fastq.gz",
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        report="results/trimmed/{sample}_trimming_report.txt",
    params:
        outdir="results/trimmed",
    threads: 4
    log:
        "results/logs/trim/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        mkdir -p {params.outdir} $(dirname {log})
        trim_galore --cores {threads} --gzip \
            --output_dir {params.outdir} {input.r1} > {log} 2>&1
        mv {params.outdir}/{wildcards.sample}_trimmed.fq.gz {output.r1}
        mv {params.outdir}/{wildcards.sample}.fastq.gz_trimming_report.txt {output.report}
        """


ruleorder: trim_pe > trim_se


def trimmed_fastqs(wc):
    if is_paired(wc.sample):
        return {
            "r1": f"results/trimmed/{wc.sample}_R1.fastq.gz",
            "r2": f"results/trimmed/{wc.sample}_R2.fastq.gz",
        }
    return {"r1": f"results/trimmed/{wc.sample}_R1.fastq.gz"}


rule bowtie2:
    input:
        unpack(trimmed_fastqs),
        index=multiext(REF["bt2_index"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2"),
    output:
        bam=temp("results/bam/{sample}.raw.bam"),
    params:
        index=REF["bt2_index"],
        # -X 1000 admits the long fragments that broad-mark chromatin yields;
        # the default 500 would silently discard them as discordant.
        reads=lambda wc, input: (
            f"-1 {input.r1} -2 {input.r2} -X 1000 --no-mixed --no-discordant"
            if is_paired(wc.sample) else f"-U {input.r1}"
        ),
    threads: config["threads"]["align"]
    log:
        "results/logs/bowtie2/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        mkdir -p results/bam $(dirname {log})
        # Default end-to-end sensitivity, as in the ENCODE ChIP-seq pipeline.
        # --very-sensitive costs ~4x the runtime and buys nothing here: reads below
        # MAPQ 30 are discarded in the next rule anyway, so the extra seed effort is
        # spent on alignments that are about to be thrown away.
        bowtie2 -x {params.index} {params.reads} -p {threads} 2> {log} \
          | samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """


rule filter_bam:
    """Remove the reads that generate false peaks.

    In order: keep primary alignments only, drop MAPQ < 30 (multi-mappers land in
    repeats and produce reproducible artefact peaks), require properly-paired reads
    for PE libraries, restrict to primary chromosomes (chrM has no ChIP signal but
    huge coverage, and it distorts every library-size normalisation), then mark and
    remove PCR duplicates.

    Duplicates are marked before removal so the library-complexity metrics can be
    computed from the marked file — deleting them first would make NRF/PBC
    uncomputable.
    """
    input:
        bam="results/bam/{sample}.raw.bam",
    output:
        bam="results/bam/{sample}.filtered.bam",
        bai="results/bam/{sample}.filtered.bam.bai",
        markdup_stats="results/qc/markdup/{sample}.txt",
        flagstat="results/qc/flagstat/{sample}.txt",
        idxstats="results/qc/idxstats/{sample}.txt",
    params:
        mapq=config["alignment"]["min_mapq"],
        keep_re=config["alignment"]["keep_chroms_regex"],
        # 0x904 = unmapped + secondary + supplementary. For PE also require 0x2.
        flags=lambda wc: "-f 2 -F 0x904" if is_paired(wc.sample) else "-F 0x904",
        rmdup=lambda wc: "-r" if config["alignment"]["remove_duplicates"] else "",
    threads: config["threads"]["sort"]
    log:
        "results/logs/filter/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bam results/qc/markdup results/qc/flagstat results/qc/idxstats $(dirname {log})
        tmp=$(mktemp -d results/bam/tmp.{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT

        keep=$(samtools idxstats {input.bam} | cut -f1 | grep -E '{params.keep_re}' | tr '\n' ' ')
        if [ -z "$keep" ]; then
            echo "ERROR: regex '{params.keep_re}' matched no contigs in {input.bam}" >&2
            exit 1
        fi

        samtools view -b {params.flags} -q {params.mapq} -@ {threads} \
            {input.bam} $keep > "$tmp/clean.bam" 2> {log}

        samtools sort -n -@ {threads} -o "$tmp/name.bam" "$tmp/clean.bam" 2>> {log}
        samtools fixmate -m -@ {threads} "$tmp/name.bam" "$tmp/fix.bam" 2>> {log}
        samtools sort -@ {threads} -o "$tmp/coord.bam" "$tmp/fix.bam" 2>> {log}
        samtools markdup {params.rmdup} -@ {threads} -f {output.markdup_stats} \
            "$tmp/coord.bam" {output.bam} 2>> {log}

        samtools index -@ {threads} {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        samtools idxstats {output.bam} > {output.idxstats}
        """


rule library_complexity:
    """NRF / PBC1 / PBC2 (ENCODE). Computed on the pre-dedup, position-sorted BAM:
    once duplicates are removed these are all 1.0 by construction."""
    input:
        bam="results/bam/{sample}.raw.bam",
    output:
        tsv="results/qc/complexity/{sample}.tsv",
    params:
        mapq=config["alignment"]["min_mapq"],
        paired=lambda wc: "1" if is_paired(wc.sample) else "0",
    log:
        "results/logs/complexity/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/complexity $(dirname {log})
        python3 scripts/library_complexity.py \
            --bam {input.bam} --sample {wildcards.sample} \
            --mapq {params.mapq} --paired {params.paired} \
            --out {output.tsv} > {log} 2>&1
        """


rule fragment_length:
    """Estimate the fragment length actually sequenced.

    PE: measured directly from the insert-size distribution.
    SE: there is no observed fragment, so it is inferred from strand
        cross-correlation. This number is not cosmetic — it becomes MACS2's
        --extsize for single-end libraries, and a wrong value smears or splits
        every peak.
    """
    input:
        bam="results/bam/{sample}.filtered.bam",
        bai="results/bam/{sample}.filtered.bam.bai",
    output:
        tsv="results/qc/fragment/{sample}.tsv",
    params:
        paired=lambda wc: "1" if is_paired(wc.sample) else "0",
        gsize=REF["gsize"],
    threads: 4
    log:
        "results/logs/fragment/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/fragment $(dirname {log})
        python3 scripts/fragment_length.py \
            --bam {input.bam} --sample {wildcards.sample} \
            --paired {params.paired} --gsize {params.gsize} \
            --threads {threads} --out {output.tsv} > {log} 2>&1
        """
