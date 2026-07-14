rule deeptools_matrix:
    input:
        bigwigs=factor_log2_bigwigs,
        peaks="results/peaks/consensus/{factor}.bed"
    output:
        matrix="results/deeptools/{factor}/matrix.gz",
        values="results/deeptools/{factor}/matrix.tsv",
        regions="results/deeptools/{factor}/regions.bed"
    threads: DEEPTOOLS_THREADS
    params:
        args=matrix_args,
        labels=factor_labels,
        outdir=lambda wc, output: os.path.dirname(str(output.matrix))
    log:
        "results/logs/deeptools_matrix_{factor}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        computeMatrix {params.args} -S {input.bigwigs:q} -R {input.peaks:q} \
            --samplesLabel {params.labels} \
            --sortRegions descend --sortUsing mean --numberOfProcessors {threads} \
            --outFileNameMatrix {output.values:q} \
            --outFileSortedRegions {output.regions:q} \
            -o {output.matrix:q} > {log:q} 2>&1
        """


rule deeptools_heatmap:
    input:
        matrix="results/deeptools/{factor}/matrix.gz"
    output:
        png="results/deeptools/{factor}/heatmap.png",
        pdf="results/deeptools/{factor}/heatmap.pdf"
    params:
        axis=plot_axis_args,
        title=lambda wc: quote(f"{SLUG_TO_FACTOR[wc.factor]} log2 ChIP/Input"),
        factor=factor_label,
        outdir=lambda wc, output: os.path.dirname(str(output.png)),
        color_map=HEATMAP_COLOR_MAP,
        zmin=HEATMAP_ZMIN,
        zmax=HEATMAP_ZMAX,
        dpi=PLOT_DPI
    log:
        "results/logs/deeptools_heatmap_{factor}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        plotHeatmap -m {input.matrix:q} -out {output.png:q} \
            --plotTitle {params.title} --yAxisLabel 'log2 ChIP/Input' \
            --regionsLabel {params.factor} \
            --colorMap {params.color_map:q} --zMin {params.zmin} --zMax {params.zmax} \
            --sortRegions keep --missingDataColor 0.85 \
            {params.axis} --dpi {params.dpi} > {log:q} 2>&1
        plotHeatmap -m {input.matrix:q} -out {output.pdf:q} \
            --plotTitle {params.title} --yAxisLabel 'log2 ChIP/Input' \
            --regionsLabel {params.factor} \
            --colorMap {params.color_map:q} --zMin {params.zmin} --zMax {params.zmax} \
            --sortRegions keep --missingDataColor 0.85 \
            {params.axis} >> {log:q} 2>&1
        """


rule deeptools_profile:
    input:
        matrix="results/deeptools/{factor}/matrix.gz"
    output:
        png="results/deeptools/{factor}/profile.png",
        pdf="results/deeptools/{factor}/profile.pdf",
        data="results/deeptools/{factor}/profile.tsv"
    params:
        axis=plot_axis_args,
        title=lambda wc: quote(f"{SLUG_TO_FACTOR[wc.factor]} meta-profile"),
        factor=factor_label,
        outdir=lambda wc, output: os.path.dirname(str(output.png)),
        dpi=PLOT_DPI
    log:
        "results/logs/deeptools_profile_{factor}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        plotProfile -m {input.matrix:q} -out {output.png:q} \
            --plotTitle {params.title} --yAxisLabel 'log2 ChIP/Input' \
            --regionsLabel {params.factor} \
            --plotType lines --perGroup {params.axis} --dpi {params.dpi} \
            --outFileNameData {output.data:q} > {log:q} 2>&1
        plotProfile -m {input.matrix:q} -out {output.pdf:q} \
            --plotTitle {params.title} --yAxisLabel 'log2 ChIP/Input' \
            --regionsLabel {params.factor} \
            --plotType lines --perGroup {params.axis} >> {log:q} 2>&1
        """


rule deeptools_fingerprint:
    input:
        bams=factor_bams,
        bais=factor_bais,
        blacklist=lambda wc: REF["black_list"]
    output:
        plot="results/deeptools/{factor}/qc/fingerprint.png",
        counts="results/deeptools/{factor}/qc/fingerprint.tsv",
        metrics="results/deeptools/{factor}/qc/fingerprint_metrics.tsv"
    threads: DEEPTOOLS_THREADS
    params:
        labels=factor_bam_labels,
        title=lambda wc: quote(f"{SLUG_TO_FACTOR[wc.factor]} fingerprint"),
        outdir=lambda wc, output: os.path.dirname(str(output.plot)),
        min_mapq=MIN_MAPQ
    log:
        "results/logs/deeptools_fingerprint_{factor}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        plotFingerprint -b {input.bams:q} --labels {params.labels} \
            --plotFile {output.plot:q} --outRawCounts {output.counts:q} \
            --outQualityMetrics {output.metrics:q} --plotTitle {params.title} \
            --extendReads --samFlagInclude 64 --minMappingQuality {params.min_mapq} \
            --blackListFileName {input.blacklist:q} \
            --numberOfProcessors {threads} > {log:q} 2>&1
        """


rule deeptools_peak_summary:
    input:
        bigwigs=factor_rpgc_bigwigs,
        peaks="results/peaks/consensus/{factor}.bed"
    output:
        npz="results/deeptools/{factor}/qc/peak_summary.npz",
        counts="results/deeptools/{factor}/qc/peak_signal.tsv"
    threads: DEEPTOOLS_THREADS
    params:
        labels=factor_labels,
        outdir=lambda wc, output: os.path.dirname(str(output.npz))
    log:
        "results/logs/deeptools_peak_summary_{factor}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        multiBigwigSummary BED-file -b {input.bigwigs:q} --BED {input.peaks:q} \
            --labels {params.labels} -o {output.npz:q} --outRawCounts {output.counts:q} \
            --numberOfProcessors {threads} > {log:q} 2>&1
        """


rule deeptools_correlation:
    input:
        npz="results/deeptools/{factor}/qc/peak_summary.npz"
    output:
        plot="results/deeptools/{factor}/qc/spearman_heatmap.png",
        matrix="results/deeptools/{factor}/qc/spearman.tsv"
    params:
        title=lambda wc: quote(f"{SLUG_TO_FACTOR[wc.factor]} peak-signal Spearman correlation")
    log:
        "results/logs/deeptools_correlation_{factor}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        plotCorrelation -in {input.npz:q} --corMethod spearman --whatToPlot heatmap \
            --plotNumbers --colorMap viridis --plotTitle {params.title} \
            --plotFile {output.plot:q} --outFileCorMatrix {output.matrix:q} \
            > {log:q} 2>&1
        """


rule deeptools_pca:
    input:
        npz="results/deeptools/{factor}/qc/peak_summary.npz"
    output:
        plot="results/deeptools/{factor}/qc/pca.png",
        data="results/deeptools/{factor}/qc/pca.tsv"
    params:
        title=lambda wc: quote(f"{SLUG_TO_FACTOR[wc.factor]} peak-signal PCA")
    log:
        "results/logs/deeptools_pca_{factor}.log"
    conda:
        "../envs/chipseq.yaml"
    shell:
        """
        plotPCA -in {input.npz:q} --transpose --plotTitle {params.title} \
            --plotFile {output.plot:q} --outFileNameData {output.data:q} \
            > {log:q} 2>&1
        """
