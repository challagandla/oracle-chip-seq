# FASTQ retrieval from SRA.
#
# Declared as a rule so the FASTQs are workflow outputs rather than something the
# user must stage by hand: `snakemake` on a clean checkout reproduces the run from
# the accessions in samples.tsv alone.
#
# fetch_sra.py cross-checks the layout it actually receives against the layout the
# sample table declares and fails on a mismatch. That check is not defensive
# boilerplate — SRA's own metadata declares every replicate-1 run of this study
# SINGLE, and they are 2x150 paired (92.7% concordant alignment). Trusting the
# metadata would have silently discarded R2 from half the cohort.


rule fetch_sra_pe:
    output:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
    wildcard_constraints:
        sample=_PE_RE,
    params:
        samples=config["samples"],
    threads: 8
    log:
        "results/logs/sra/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p data/raw $(dirname {log})
        python3 scripts/fetch_sra.py --samples {params.samples} \
            --outdir data/raw --threads {threads} \
            --only {wildcards.sample} > {log} 2>&1
        """


rule fetch_sra_se:
    output:
        r1="data/raw/{sample}.fastq.gz",
    wildcard_constraints:
        sample=_SE_RE,
    params:
        samples=config["samples"],
    threads: 8
    log:
        "results/logs/sra/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p data/raw $(dirname {log})
        python3 scripts/fetch_sra.py --samples {params.samples} \
            --outdir data/raw --threads {threads} \
            --only {wildcards.sample} > {log} 2>&1
        """
