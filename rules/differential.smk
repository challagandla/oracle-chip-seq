# Differential binding.
#
# Counting and normalisation are both driven by the peak mode:
#
#   count regions   narrow -> fixed windows on the MACS2 summit (registry: summits)
#                   broad  -> the full called domain
#
#   size factors    narrow -> median-of-ratios on the peak counts themselves. Valid
#                             because most peaks are unchanged between conditions.
#                   broad  -> estimated from genome-wide background bins instead.
#                             A repressive mark can move globally (this is the whole
#                             point of, say, an EZH2 inhibitor), and when it does,
#                             normalising on reads-in-peaks divides the biology out
#                             and returns nothing. csaw makes the same argument.


def _target_summits(wc):
    """MACS2 summit files for a target. Narrow only — broad calling emits none."""
    if REG.is_broad(wc.target):
        return []
    return [
        f"results/peaks/raw/{s}_summits.bed"
        for s in CHIP[CHIP.target == wc.target].sample_id
    ]


rule count_regions:
    input:
        peaks="results/peaks/consensus/{target}.bed",
        summits=_target_summits,
        chrom_sizes=REF["chrom_sizes"],
    output:
        bed="results/differential/{target}/count_regions.bed",
        saf="results/differential/{target}/count_regions.saf",
    params:
        registry=config["mark_registry"],
    log:
        "results/logs/differential/{target}_regions.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/differential/{wildcards.target} $(dirname {log})
        python3 scripts/make_count_regions.py \
            --peaks {input.peaks} --target {wildcards.target} \
            --summits {input.summits} \
            --chrom-sizes {input.chrom_sizes} \
            --registry {params.registry} \
            --out-bed {output.bed} --out-saf {output.saf} > {log} 2>&1
        """


rule background_bins:
    """Genome-wide 10 kb bins, used to derive size factors for broad marks."""
    input:
        chrom_sizes=REF["chrom_sizes"],
        blacklist=REF["blacklist"],
    output:
        saf="results/differential/background_bins.saf",
    params:
        keep_re=config["alignment"]["keep_chroms_regex"],
    log:
        "results/logs/differential/background_bins.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/differential $(dirname {log})
        awk -v OFS='\t' '$1 ~ /{params.keep_re}/ {{print $1,$2}}' {input.chrom_sizes} \
          | sort -k1,1 > results/differential/.genome.txt
        bedtools makewindows -g results/differential/.genome.txt -w 10000 \
          | bedtools intersect -v -a - -b {input.blacklist} \
          | awk -v OFS='\t' 'BEGIN{{print "GeneID","Chr","Start","End","Strand"}}
                             {{print "bin_"NR,$1,$2+1,$3,"+"}}' > {output.saf}
        echo "$(( $(wc -l < {output.saf}) - 1 )) background bins" > {log}
        """


def _count_bams(wc):
    return expand("results/bam/{s}.filtered.bam", s=list(CHIP[CHIP.target == wc.target].sample_id))


rule count_peaks:
    """featureCounts over the peak regions.

    Single- and paired-end libraries are counted in separate invocations and then
    joined, because `-p --countReadPairs` is a run-level flag. Counting a PE library
    as single-end would double every count relative to the SE libraries and put a
    2x artefact straight into the contrast — the exact failure this dataset's mixed
    layout invites.
    """
    input:
        saf="results/differential/{target}/count_regions.saf",
        bams=_count_bams,
    output:
        counts="results/differential/{target}/counts.tsv",
    params:
        se=lambda wc: " ".join(
            f"results/bam/{s}.filtered.bam"
            for s in CHIP[CHIP.target == wc.target].sample_id if not is_paired(s)
        ),
        pe=lambda wc: " ".join(
            f"results/bam/{s}.filtered.bam"
            for s in CHIP[CHIP.target == wc.target].sample_id if is_paired(s)
        ),
    threads: 8
    log:
        "results/logs/differential/{target}_counts.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/differential/{wildcards.target} $(dirname {log})
        python3 scripts/count_features.py \
            --saf {input.saf} --out {output.counts} --threads {threads} \
            --se {params.se} --pe {params.pe} > {log} 2>&1
        """


rule count_background:
    input:
        saf="results/differential/background_bins.saf",
        bams=_count_bams,
    output:
        counts="results/differential/{target}/background_counts.tsv",
    params:
        se=lambda wc: " ".join(
            f"results/bam/{s}.filtered.bam"
            for s in CHIP[CHIP.target == wc.target].sample_id if not is_paired(s)
        ),
        pe=lambda wc: " ".join(
            f"results/bam/{s}.filtered.bam"
            for s in CHIP[CHIP.target == wc.target].sample_id if is_paired(s)
        ),
    threads: 8
    log:
        "results/logs/differential/{target}_background.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/differential/{wildcards.target} $(dirname {log})
        python3 scripts/count_features.py \
            --saf {input.saf} --out {output.counts} --threads {threads} \
            --se {params.se} --pe {params.pe} > {log} 2>&1
        """


rule sample_sheet:
    output:
        csv="results/differential/{target}/coldata.csv",
    run:
        sub = CHIP[CHIP.target == wildcards.target].copy()
        sub = sub.sort_values(["condition", "replicate"])
        Path(output.csv).parent.mkdir(parents=True, exist_ok=True)
        sub[["sample_id", "condition", "replicate", "layout", "target"]].to_csv(
            output.csv, index=False
        )


rule differential_binding:
    input:
        counts="results/differential/{target}/counts.tsv",
        background="results/differential/{target}/background_counts.tsv",
        coldata="results/differential/{target}/coldata.csv",
        regions="results/differential/{target}/count_regions.bed",
    output:
        results="results/differential/{target}/results.tsv",
        up="results/differential/{target}/up.bed",
        down="results/differential/{target}/down.bed",
        norm="results/differential/{target}/normalized_counts.tsv",
        rlog="results/differential/{target}/rlog.tsv",
        sizefactors="results/differential/{target}/size_factors.tsv",
    params:
        registry=config["mark_registry"],
        ref=REF_LEVEL,
        trt=TRT_LEVEL,
        fdr=config["differential"]["fdr"],
        lfc=config["differential"]["min_lfc"],
    log:
        "results/logs/differential/{target}_deseq2.log",
    conda:
        "../envs/r_analysis.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/differential/{wildcards.target} $(dirname {log})
        Rscript analysis/differential_binding.R \
            --counts {input.counts} \
            --background {input.background} \
            --coldata {input.coldata} \
            --regions {input.regions} \
            --target {wildcards.target} \
            --registry {params.registry} \
            --reference {params.ref} --treatment {params.trt} \
            --fdr {params.fdr} --lfc {params.lfc} \
            --outdir results/differential/{wildcards.target} > {log} 2>&1
        """
