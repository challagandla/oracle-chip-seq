# Peak calling and replicate reproducibility.
#
# Everything in this file branches on the peak mode from config/mark_registry.yaml:
#
#   caller         narrow -> MACS2 default;  broad -> MACS2 --broad --broad-cutoff
#   width filter   narrow -> drop peaks wider than a few kb (they are artefacts);
#                  broad  -> a 500 kb H3K27me3 domain is the signal, not an artefact
#   reproducibility
#                  narrow -> IDR. It models the rank consistency of a ranked peak
#                            list, which needs peaks with well-separated scores.
#                  broad  -> naive overlap against pooled-replicate peaks. IDR on
#                            broad peaks is not valid: the peaks are wide, their
#                            scores are compressed, and the rank-reproducibility
#                            model does not hold. ENCODE uses overlap for broad
#                            marks for exactly this reason.

NARROW_CHIP = [s for s in CHIP_SAMPLES if not REG.is_broad(target_of(s))]
BROAD_CHIP = [s for s in CHIP_SAMPLES if REG.is_broad(target_of(s))]

# target/condition groups that get pooled
GROUPS = [(t, c) for t in TARGETS for c in CONDITIONS]


def _frag(sample):
    return f"results/qc/fragment/{sample}.tsv"


def _macs2_cmd(sample, treat_bam, ctrl_bam, outdir, name, frag_tsv):
    """Build the MACS2 invocation for one sample from the registry.

    Kept in Python rather than shell so the peak mode, the q-value and the
    single-end fragment size are all resolved in one place.
    """
    return (
        f'frag=$(awk \'NR==2 {{print $3}}\' {frag_tsv}); '
        f"macs2 callpeak -t {treat_bam} -c {ctrl_bam} "
        f"--name {name} --outdir {outdir} "
        f"$(python3 scripts/macs2_args.py --target {target_of(sample)} "
        f"--gsize {REF['gsize']} --paired {int(is_paired(sample))} --extsize $frag)"
    )


rule macs2_narrow:
    input:
        bam="results/bam/{sample}.filtered.bam",
        bai="results/bam/{sample}.filtered.bam.bai",
        control=lambda wc: f"results/bam/{control_of(wc.sample)}.filtered.bam",
        frag=lambda wc: _frag(wc.sample),
    output:
        peaks="results/peaks/raw/{sample}_peaks.narrowPeak",
        summits="results/peaks/raw/{sample}_summits.bed",
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in NARROW_CHIP) if NARROW_CHIP else "$^",
    params:
        cmd=lambda wc, input: _macs2_cmd(
            wc.sample, input.bam, input.control, "results/peaks/raw", wc.sample, input.frag
        ),
    log:
        "results/logs/macs2/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        mkdir -p results/peaks/raw $(dirname {log})
        {params.cmd} > {log} 2>&1
        """


rule macs2_broad:
    input:
        bam="results/bam/{sample}.filtered.bam",
        bai="results/bam/{sample}.filtered.bam.bai",
        control=lambda wc: f"results/bam/{control_of(wc.sample)}.filtered.bam",
        frag=lambda wc: _frag(wc.sample),
    output:
        peaks="results/peaks/raw/{sample}_peaks.broadPeak",
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in BROAD_CHIP) if BROAD_CHIP else "$^",
    params:
        cmd=lambda wc, input: _macs2_cmd(
            wc.sample, input.bam, input.control, "results/peaks/raw", wc.sample, input.frag
        ),
    log:
        "results/logs/macs2/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        mkdir -p results/peaks/raw $(dirname {log})
        {params.cmd} > {log} 2>&1
        """


rule macs2_relaxed:
    """A permissive peak list used only as IDR input.

    IDR needs to see the noisy tail of the ranking to fit its two-component model;
    feeding it an already q-filtered list truncates the irreproducible component
    and the fit degenerates. ENCODE calls peaks at p < 0.01 for this purpose. Narrow
    targets only, since broad targets do not go through IDR.
    """
    input:
        bam="results/bam/{sample}.filtered.bam",
        bai="results/bam/{sample}.filtered.bam.bai",
        control=lambda wc: f"results/bam/{control_of(wc.sample)}.filtered.bam",
        frag=lambda wc: _frag(wc.sample),
    output:
        peaks="results/peaks/relaxed/{sample}_peaks.narrowPeak",
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in NARROW_CHIP) if NARROW_CHIP else "$^",
    params:
        fmt=lambda wc: "BAMPE" if is_paired(wc.sample) else "BAM",
        nomodel=lambda wc: "" if is_paired(wc.sample) else "--nomodel --extsize",
        gsize=REF["gsize"],
    log:
        "results/logs/macs2_relaxed/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/peaks/relaxed $(dirname {log})
        frag=$(awk 'NR==2 {{print $3}}' {input.frag})
        extra=""
        if [ -n "{params.nomodel}" ]; then extra="{params.nomodel} $frag"; fi
        macs2 callpeak -t {input.bam} -c {input.control} \
            --name {wildcards.sample} --outdir results/peaks/relaxed \
            --format {params.fmt} --gsize {params.gsize} \
            --pvalue 0.01 --keep-dup all $extra > {log} 2>&1
        """


rule filter_peaks:
    """Blacklist removal plus a peak-mode-aware width filter.

    The width filter is the point of the registry: a 40 kb "peak" in a CTCF track
    is a mapping artefact and must go, while a 40 kb H3K27me3 peak is a Polycomb
    domain and must stay. One global width cutoff cannot serve both.
    """
    input:
        peaks="results/peaks/raw/{sample}_peaks.{ext}",
        blacklist=REF["blacklist"],
    output:
        peaks="results/peaks/{sample}_peaks.{ext}",
    params:
        max_width=lambda wc: REG.get(target_of(wc.sample))["qc"]["max_peak_width"],
    log:
        "results/logs/filter_peaks/{sample}_{ext}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/peaks $(dirname {log})
        bedtools intersect -v -a {input.peaks} -b {input.blacklist} \
          | awk -v OFS='\t' -v maxw={params.max_width} '($3-$2) <= maxw' \
          | sort -k1,1 -k2,2n > {output.peaks} 2> {log}
        n=$(wc -l < {output.peaks})
        echo "{wildcards.sample}: $n peaks after blacklist and width<={params.max_width}" >> {log}
        """


# ------------------------------------------------------------ pooled replicates

rule pool_bam:
    """Merged replicate BAM per target x condition. Pooled depth is what makes the
    naive-overlap reproducible set for broad marks usable, and it also gives the
    browser tracks a signal track that is not replicate-specific."""
    input:
        bams=lambda wc: expand(
            "results/bam/{s}.filtered.bam", s=reps_for(wc.target, wc.condition)
        ),
    output:
        bam="results/bam/pooled/{target}_{condition}.bam",
        bai="results/bam/pooled/{target}_{condition}.bam.bai",
    threads: config["threads"]["sort"]
    log:
        "results/logs/pool/{target}_{condition}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bam/pooled $(dirname {log})
        samtools merge -f -@ {threads} {output.bam} {input.bams} > {log} 2>&1
        samtools index -@ {threads} {output.bam}
        """


rule pool_control:
    input:
        bams=lambda wc: expand(
            "results/bam/{s}.filtered.bam",
            s=sorted({control_of(s) for s in reps_for(wc.target, wc.condition)}),
        ),
    output:
        bam="results/bam/pooled/control_{target}_{condition}.bam",
        bai="results/bam/pooled/control_{target}_{condition}.bam.bai",
    threads: config["threads"]["sort"]
    log:
        "results/logs/pool/control_{target}_{condition}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bam/pooled $(dirname {log})
        samtools merge -f -@ {threads} {output.bam} {input.bams} > {log} 2>&1
        samtools index -@ {threads} {output.bam}
        """


rule macs2_pooled:
    input:
        bam="results/bam/pooled/{target}_{condition}.bam",
        bai="results/bam/pooled/{target}_{condition}.bam.bai",
        control="results/bam/pooled/control_{target}_{condition}.bam",
        frag=lambda wc: _frag(reps_for(wc.target, wc.condition)[0]),
    output:
        peaks="results/peaks/pooled/{target}_{condition}_peaks.{ext}",
    params:
        target=lambda wc: wc.target,
        # The pooled BAM mixes SE and PE replicates, so it cannot be read as BAMPE.
        # Treat the pool as single-end and extend by the fragment estimate: this is
        # the only layout-agnostic option, and it is why the per-replicate calls
        # above (which do use BAMPE where available) remain the primary peak set.
        gsize=REF["gsize"],
    log:
        "results/logs/macs2_pooled/{target}_{condition}_{ext}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/peaks/pooled $(dirname {log})
        frag=$(awk 'NR==2 {{print $3}}' {input.frag})
        args=$(python3 scripts/macs2_args.py --target {params.target} \
                 --gsize {params.gsize} --paired 0 --extsize $frag)
        macs2 callpeak -t {input.bam} -c {input.control} \
            --name {wildcards.target}_{wildcards.condition} \
            --outdir results/peaks/pooled $args > {log} 2>&1
        """


# ------------------------------------------------------- reproducibility filters

rule idr_narrow:
    """IDR across true replicates. Narrow targets only."""
    input:
        reps=lambda wc: expand(
            "results/peaks/relaxed/{s}_peaks.narrowPeak", s=reps_for(wc.target, wc.condition)
        ),
        blacklist=REF["blacklist"],
    output:
        idr="results/peaks/idr/{target}_{condition}.idr.txt",
        bed="results/peaks/reproducible/{target}_{condition}.bed",
    wildcard_constraints:
        target="|".join(re.escape(t) for t in IDR_TARGETS) if IDR_TARGETS else "$^",
    params:
        thr=config["peaks"]["idr_threshold"],
    log:
        "results/logs/idr/{target}_{condition}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/peaks/idr results/peaks/reproducible $(dirname {log})
        idr --samples {input.reps} \
            --input-file-type narrowPeak --rank p.value \
            --output-file {output.idr} \
            --idr-threshold {params.thr} \
            --plot > {log} 2>&1

        # IDR column 5 is int(-125*log2(IDR)); the global IDR itself is column 12
        # as -log10(IDR). Threshold on that directly rather than on the scaled score.
        thr_log10=$(python3 -c "import math;print(-math.log10({params.thr}))")
        awk -v OFS='\t' -v t="$thr_log10" '$12 >= t {{print $1,$2,$3,"peak_"NR,$5,"."}}' {output.idr} \
          | bedtools intersect -v -a - -b {input.blacklist} \
          | sort -k1,1 -k2,2n | bedtools merge -c 4,5 -o first,max \
          | awk -v OFS='\t' '{{print $1,$2,$3,$4,$5,"."}}' > {output.bed}
        echo "IDR<={params.thr}: $(wc -l < {output.bed}) reproducible peaks" >> {log}
        """


rule overlap_broad:
    """Naive overlap for broad targets: start from pooled-replicate peaks and keep
    only those supported by a peak in every individual replicate. This is ENCODE's
    reproducibility criterion for broad marks."""
    input:
        pooled=lambda wc: f"results/peaks/pooled/{wc.target}_{wc.condition}_peaks.{REG.peak_ext(wc.target)}",
        reps=lambda wc: [sample_peaks(s) for s in reps_for(wc.target, wc.condition)],
        blacklist=REF["blacklist"],
    output:
        bed="results/peaks/reproducible/{target}_{condition}.bed",
    wildcard_constraints:
        target="|".join(re.escape(t) for t in OVERLAP_TARGETS) if OVERLAP_TARGETS else "$^",
    params:
        frac=config["peaks"]["overlap_fraction"],
    log:
        "results/logs/overlap/{target}_{condition}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/peaks/reproducible $(dirname {log})
        tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

        bedtools intersect -v -a {input.pooled} -b {input.blacklist} \
          | sort -k1,1 -k2,2n | cut -f1-3 > "$tmp/cur.bed"

        for rep in {input.reps}; do
            sort -k1,1 -k2,2n "$rep" | cut -f1-3 > "$tmp/rep.bed"
            bedtools intersect -u -f {params.frac} -a "$tmp/cur.bed" -b "$tmp/rep.bed" > "$tmp/next.bed"
            mv "$tmp/next.bed" "$tmp/cur.bed"
        done

        bedtools merge -i "$tmp/cur.bed" \
          | awk -v OFS='\t' '{{print $1,$2,$3,"peak_"NR,"0","."}}' > {output.bed}
        echo "naive overlap (f={params.frac}): $(wc -l < {output.bed}) reproducible peaks" >> {log}
        """


rule consensus_peaks:
    """Union of the reproducible peaks from every condition. This is the interval
    set the differential test is run over — built per target, never across targets:
    merging CTCF and H3K27me3 intervals into one 'consensus' would be meaningless."""
    input:
        beds=lambda wc: [
            f"results/peaks/reproducible/{wc.target}_{c}.bed" for c in CONDITIONS
        ],
    output:
        bed="results/peaks/consensus/{target}.bed",
    log:
        "results/logs/consensus/{target}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/peaks/consensus $(dirname {log})
        cat {input.beds} | sort -k1,1 -k2,2n | bedtools merge \
          | awk -v OFS='\t' '{{print $1,$2,$3,"{wildcards.target}_"NR,"0","."}}' > {output.bed}
        echo "{wildcards.target}: $(wc -l < {output.bed}) consensus peaks" >> {log}
        """


rule frip:
    """Fraction of reads in peaks — the headline ChIP enrichment metric. Judged
    against a per-mark threshold, because 2% is a failure for H3K4me3 and entirely
    normal for H3K9me3."""
    input:
        bam="results/bam/{sample}.filtered.bam",
        bai="results/bam/{sample}.filtered.bam.bai",
        peaks=lambda wc: sample_peaks(wc.sample),
    output:
        tsv="results/qc/frip/{sample}.tsv",
    log:
        "results/logs/frip/{sample}.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/frip $(dirname {log})
        total=$(samtools view -c {input.bam})
        if [ "$total" -eq 0 ]; then echo "empty BAM" >&2; exit 1; fi
        inpeak=$(bedtools sort -i {input.peaks} | bedtools merge \
                 | bedtools intersect -u -a {input.bam} -b - -ubam | samtools view -c)
        frip=$(python3 -c "print(f'{{$inpeak/$total:.4f}}')")
        printf 'sample\ttotal_reads\treads_in_peaks\tFRiP\tn_peaks\n' > {output.tsv}
        printf '%s\t%s\t%s\t%s\t%s\n' "{wildcards.sample}" "$total" "$inpeak" "$frip" \
            "$(wc -l < {input.peaks})" >> {output.tsv}
        cat {output.tsv} > {log}
        """
