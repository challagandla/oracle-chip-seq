# Coverage tracks and signal-level QC.
#
# Two kinds of bigwig are produced deliberately:
#
#   *.cpm.bw       depth-normalised coverage of one library. Fine for looking at
#                  one track, misleading across antibodies — a weak antibody and a
#                  strong one give different dynamic ranges regardless of biology.
#   *.log2ratio.bw log2(ChIP / matched Input). Divides out chromatin accessibility,
#                  copy number and sonication bias. This is the track to compare
#                  across samples and the one the figures use.
#
# Single-end reads are extended to the fragment length estimated by
# cross-correlation; paired-end reads use the observed fragment. Skipping the
# extension for SE data halves the apparent peak width.


rule bigwig_cpm:
    input:
        bam="results/bam/{sample}.filtered.bam",
        bai="results/bam/{sample}.filtered.bam.bai",
        frag="results/qc/fragment/{sample}.tsv",
        blacklist=REF["blacklist"],
    output:
        bw="results/bigwig/{sample}.cpm.bw",
    params:
        binsize=config["signal"]["bin_size"],
        smooth=config["signal"]["smooth_length"],
        egs=REF["effective_genome_size"],
    threads: config["threads"]["deeptools"]
    log:
        "results/logs/bigwig/{sample}.cpm.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bigwig $(dirname {log})
        frag=$(awk 'NR==2 {{print $3}}' {input.frag})
        layout=$(awk 'NR==2 {{print $2}}' {input.frag})
        if [ "$layout" = "paired" ]; then ext="--extendReads"; else ext="--extendReads $frag"; fi
        bamCoverage -b {input.bam} -o {output.bw} \
            --normalizeUsing CPM \
            --effectiveGenomeSize {params.egs} \
            --binSize {params.binsize} --smoothLength {params.smooth} \
            --blackListFileName {input.blacklist} \
            $ext -p {threads} > {log} 2>&1
        """


rule bigwig_log2ratio:
    input:
        chip="results/bam/{sample}.filtered.bam",
        chip_bai="results/bam/{sample}.filtered.bam.bai",
        control=lambda wc: f"results/bam/{control_of(wc.sample)}.filtered.bam",
        control_bai=lambda wc: f"results/bam/{control_of(wc.sample)}.filtered.bam.bai",
        frag="results/qc/fragment/{sample}.tsv",
        blacklist=REF["blacklist"],
    output:
        bw="results/bigwig/{sample}.log2ratio.bw",
    params:
        binsize=config["signal"]["bin_size"],
        egs=REF["effective_genome_size"],
    threads: config["threads"]["deeptools"]
    log:
        "results/logs/bigwig/{sample}.log2.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/bigwig $(dirname {log})
        frag=$(awk 'NR==2 {{print $3}}' {input.frag})
        layout=$(awk 'NR==2 {{print $2}}' {input.frag})
        if [ "$layout" = "paired" ]; then ext="--extendReads"; else ext="--extendReads $frag"; fi
        bamCompare -b1 {input.chip} -b2 {input.control} -o {output.bw} \
            --operation log2 --scaleFactorsMethod SES \
            --pseudocount 1 \
            --binSize {params.binsize} \
            --effectiveGenomeSize {params.egs} \
            --blackListFileName {input.blacklist} \
            $ext -p {threads} > {log} 2>&1
        """


rule fingerprint:
    """plotFingerprint: cumulative read distribution across genome bins.

    A successful ChIP is strongly skewed (a small fraction of bins holds most
    reads) and its Input is nearly diagonal. If a ChIP curve tracks its Input,
    the immunoprecipitation did not enrich anything and no amount of peak calling
    will rescue it. Broad marks sit between the two by nature, which is why this
    plot is read per mark rather than against one universal shape.
    """
    input:
        bams=expand("results/bam/{s}.filtered.bam", s=ALL_SAMPLES),
        bais=expand("results/bam/{s}.filtered.bam.bai", s=ALL_SAMPLES),
        blacklist=REF["blacklist"],
    output:
        pdf="results/qc/fingerprint/fingerprint.pdf",
        metrics="results/qc/fingerprint/metrics.tsv",
        counts="results/qc/fingerprint/raw_counts.tsv",
    params:
        labels=" ".join(ALL_SAMPLES),
    threads: config["threads"]["deeptools"]
    log:
        "results/logs/fingerprint.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/fingerprint $(dirname {log})
        plotFingerprint -b {input.bams} --labels {params.labels} \
            --blackListFileName {input.blacklist} \
            --minMappingQuality 30 --skipZeros \
            --numberOfSamples 50000 \
            --outQualityMetrics {output.metrics} \
            --outRawCounts {output.counts} \
            --plotFile {output.pdf} \
            --plotTitle "ChIP enrichment fingerprint" \
            -p {threads} > {log} 2>&1
        """


rule multibigwig_summary:
    """Genome-wide binned signal matrix, used for the correlation heatmap and PCA.

    Built from the log2(ChIP/Input) tracks, not raw coverage: on raw coverage the
    dominant axis of variation is sequencing depth and antibody efficiency, and
    samples cluster by those rather than by biology.
    """
    input:
        bws=expand("results/bigwig/{s}.log2ratio.bw", s=CHIP_SAMPLES),
        blacklist=REF["blacklist"],
    output:
        npz="results/qc/correlation/signal_matrix.npz",
        tab="results/qc/correlation/signal_matrix.tsv",
    params:
        labels=" ".join(CHIP_SAMPLES),
    threads: config["threads"]["deeptools"]
    log:
        "results/logs/multibigwig_summary.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/correlation $(dirname {log})
        multiBigwigSummary bins -b {input.bws} --labels {params.labels} \
            --binSize 10000 \
            --blackListFileName {input.blacklist} \
            -o {output.npz} --outRawCounts {output.tab} \
            -p {threads} > {log} 2>&1
        """


rule correlation_plots:
    input:
        npz="results/qc/correlation/signal_matrix.npz",
    output:
        heatmap="results/qc/correlation/spearman_heatmap.pdf",
        pca="results/qc/correlation/pca.pdf",
        matrix="results/qc/correlation/spearman_matrix.tsv",
    log:
        "results/logs/correlation_plots.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/qc/correlation $(dirname {log})
        plotCorrelation -in {input.npz} \
            --corMethod spearman --whatToPlot heatmap --skipZeros \
            --plotNumbers --colorMap RdYlBu_r \
            --plotTitle "Spearman correlation, log2(ChIP/Input) 10 kb bins" \
            --outFileCorMatrix {output.matrix} \
            -o {output.heatmap} > {log} 2>&1
        plotPCA -in {input.npz} \
            --plotTitle "PCA, log2(ChIP/Input) 10 kb bins" \
            -o {output.pca} >> {log} 2>&1
        """


# ------------------------------------------------------------------- profiles
# Profile geometry is chosen by peak mode, and this is where the choice matters
# most visibly. Anchoring a 100 kb H3K27me3 domain on its "centre" produces a flat
# line: the domain is wider than any sensible window. Conversely, scaling a 300 bp
# CTCF peak to a fixed body length destroys the summit that is the entire point.
#
#   narrow / tss   -> reference-point, +/- 3 kb around peak centre (or TSS)
#   broad          -> scale-regions: body scaled to 10 kb with 5 kb flanks


def _profile_bws(wc):
    return expand(
        "results/bigwig/{s}.log2ratio.bw", s=list(CHIP[CHIP.target == wc.target].sample_id)
    )


def _profile_labels(wc):
    return " ".join(list(CHIP[CHIP.target == wc.target].sample_id))


rule compute_matrix:
    input:
        bws=_profile_bws,
        regions=lambda wc: (
            "results/annotation/tss.bed"
            if REG.get(wc.target)["profile"] == "tss"
            else f"results/peaks/consensus/{wc.target}.bed"
        ),
    output:
        matrix="results/profiles/{target}.matrix.gz",
    params:
        labels=_profile_labels,
        mode=lambda wc: REG.get(wc.target)["profile"],
    threads: config["threads"]["deeptools"]
    log:
        "results/logs/profiles/{target}_matrix.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/profiles $(dirname {log})
        case "{params.mode}" in
          scale_regions)
            computeMatrix scale-regions -S {input.bws} -R {input.regions} \
                --samplesLabel {params.labels} \
                --regionBodyLength 10000 \
                --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
                --binSize 100 --skipZeros --missingDataAsZero \
                -o {output.matrix} -p {threads} > {log} 2>&1
            ;;
          *)
            computeMatrix reference-point -S {input.bws} -R {input.regions} \
                --samplesLabel {params.labels} \
                --referencePoint center \
                --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
                --binSize 50 --skipZeros --missingDataAsZero \
                -o {output.matrix} -p {threads} > {log} 2>&1
            ;;
        esac
        """


rule plot_profile:
    input:
        matrix="results/profiles/{target}.matrix.gz",
    output:
        heatmap="results/profiles/{target}_heatmap.pdf",
        profile="results/profiles/{target}_profile.pdf",
    params:
        mode=lambda wc: REG.get(wc.target)["profile"],
        anchor=lambda wc: (
            "TSS" if REG.get(wc.target)["profile"] == "tss" else "peak centre"
        ),
    log:
        "results/logs/profiles/{target}_plot.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/profiles $(dirname {log})
        if [ "{params.mode}" = "scale_regions" ]; then
            startlab="start"; endlab="end"
        else
            startlab="{params.anchor}"; endlab=""
        fi
        plotHeatmap -m {input.matrix} -o {output.heatmap} \
            --colorMap RdBu_r --missingDataColor 1 \
            --startLabel "$startlab" --endLabel "$endlab" \
            --heatmapHeight 12 --heatmapWidth 4 \
            --yAxisLabel "log2 ChIP/Input" \
            --plotTitle "{wildcards.target} ({params.mode})" \
            --legendLocation none > {log} 2>&1
        plotProfile -m {input.matrix} -o {output.profile} \
            --perGroup --startLabel "$startlab" --endLabel "$endlab" \
            --yAxisLabel "log2 ChIP/Input" \
            --plotTitle "{wildcards.target} meta-profile" >> {log} 2>&1
        """
