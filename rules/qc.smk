# Read-level QC and the mark-aware acceptance gate.
#
# FastQC emits one report per FASTQ, so a single-end sample produces one zip and a
# paired-end sample produces two. Snakemake needs a static output list, hence the
# SE/PE rule pairs below rather than one rule with a computed output.


rule fastqc_raw_se:
    input:
        r1="data/raw/{sample}.fastq.gz",
    output:
        zip="results/qc/fastqc/raw/{sample}_R1_fastqc.zip",
        html="results/qc/fastqc/raw/{sample}_R1_fastqc.html",
    wildcard_constraints:
        sample=_SE_RE,
    threads: 2
    log:
        "results/logs/fastqc/{sample}.raw.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/fastqc/raw $(dirname {log})
        tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
        # FastQC names outputs after the input file, so stage under the canonical
        # {sample}_R1 name to get a predictable output path.
        ln -sf "$(readlink -f {input.r1})" "$tmp/{wildcards.sample}_R1.fastq.gz"
        fastqc -o results/qc/fastqc/raw -t {threads} "$tmp/{wildcards.sample}_R1.fastq.gz" > {log} 2>&1
        """


rule fastqc_raw_pe:
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
    output:
        zip1="results/qc/fastqc/raw/{sample}_R1_fastqc.zip",
        zip2="results/qc/fastqc/raw/{sample}_R2_fastqc.zip",
        html1="results/qc/fastqc/raw/{sample}_R1_fastqc.html",
        html2="results/qc/fastqc/raw/{sample}_R2_fastqc.html",
    wildcard_constraints:
        sample=_PE_RE,
    threads: 4
    log:
        "results/logs/fastqc/{sample}.raw.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/fastqc/raw $(dirname {log})
        fastqc -o results/qc/fastqc/raw -t {threads} {input.r1} {input.r2} > {log} 2>&1
        """


rule fastqc_trimmed_se:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
    output:
        zip="results/qc/fastqc/trimmed/{sample}_R1_fastqc.zip",
    wildcard_constraints:
        sample=_SE_RE,
    threads: 2
    log:
        "results/logs/fastqc/{sample}.trimmed.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/fastqc/trimmed $(dirname {log})
        fastqc -o results/qc/fastqc/trimmed -t {threads} {input.r1} > {log} 2>&1
        """


rule fastqc_trimmed_pe:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
    output:
        zip1="results/qc/fastqc/trimmed/{sample}_R1_fastqc.zip",
        zip2="results/qc/fastqc/trimmed/{sample}_R2_fastqc.zip",
    wildcard_constraints:
        sample=_PE_RE,
    threads: 4
    log:
        "results/logs/fastqc/{sample}.trimmed.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/fastqc/trimmed $(dirname {log})
        fastqc -o results/qc/fastqc/trimmed -t {threads} {input.r1} {input.r2} > {log} 2>&1
        """


rule qc_gate:
    """Collect every per-sample metric and judge it against the threshold for that
    mark. The thresholds live in the registry because one global cutoff is wrong by
    construction: FRiP 2% fails H3K4me3 and passes H3K9me3.

    The gate reports rather than aborts. A shallow library is still worth looking
    at — it just must not be quietly presented as though it were well powered.
    """
    input:
        frip=expand("results/qc/frip/{s}.tsv", s=CHIP_SAMPLES),
        complexity=expand("results/qc/complexity/{s}.tsv", s=ALL_SAMPLES),
        fragment=expand("results/qc/fragment/{s}.tsv", s=ALL_SAMPLES),
        flagstat=expand("results/qc/flagstat/{s}.txt", s=ALL_SAMPLES),
        peaks=all_chip_peaks(),
    output:
        tsv="results/qc/qc_gate.tsv",
        md="results/qc/qc_gate.md",
    params:
        registry=config["mark_registry"],
        samples=config["samples"],
    log:
        "results/logs/qc_gate.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc $(dirname {log})
        python3 scripts/qc_gate.py \
            --samples {params.samples} \
            --registry {params.registry} \
            --config config.yaml \
            --qcdir results/qc \
            --tsv {output.tsv} --md {output.md} > {log} 2>&1
        """


rule multiqc:
    input:
        fastqc_raw=[
            f"results/qc/fastqc/raw/{s}_{r}_fastqc.zip"
            for s in ALL_SAMPLES for r in reads_of(s)
        ],
        fastqc_trimmed=[
            f"results/qc/fastqc/trimmed/{s}_{r}_fastqc.zip"
            for s in ALL_SAMPLES for r in reads_of(s)
        ],
        flagstat=expand("results/qc/flagstat/{s}.txt", s=ALL_SAMPLES),
        markdup=expand("results/qc/markdup/{s}.txt", s=ALL_SAMPLES),
        bowtie2=expand("results/logs/bowtie2/{s}.log", s=ALL_SAMPLES),
        gate="results/qc/qc_gate.tsv",
    output:
        html="results/qc/multiqc/multiqc_report.html",
    log:
        "results/logs/multiqc.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/multiqc $(dirname {log})
        multiqc results/qc results/logs/bowtie2 results/trimmed \
            --outdir results/qc/multiqc --filename multiqc_report.html \
            --force > {log} 2>&1
        """
