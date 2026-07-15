rule build_sample_sheets:
    output:
        diffbind="results/diffbind/sample_sheet.csv",
        resolved_config="results/provenance/resolved_config.yaml"
    params:
        config_digest=CONFIG_DIGEST
    log:
        "results/logs/build_sample_sheets.log"
    conda:
        "../envs/chipseq.yaml"
    script:
        "../scripts/build_sample_sheets.py"


rule run_diffbind:
    input:
        sheet="results/diffbind/sample_sheet.csv",
        bams=factor_bams,
        bais=factor_bais,
        peaks=factor_peaks,
        script="analysis/diffbind_analysis.R"
    output:
        summary="results/diffbind/{factor}/diffbind_summary.csv",
        plots="results/diffbind/{factor}/diffbind_plots.pdf",
        rds="results/diffbind/{factor}/diffbind.rds",
        contrast="results/diffbind/{factor}/contrast.tsv"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.summary)),
        factor=factor_label,
        numerator=quote(DIFFBIND_NUMERATOR),
        reference=quote(DIFFBIND_REFERENCE),
        narrow_summits=DIFFBIND_NARROW_SUMMITS
    log:
        "results/logs/diffbind_{factor}.log"
    conda:
        "../envs/r_analysis.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        Rscript {input.script:q} {input.sheet:q} {params.outdir:q} {params.factor} \
            {params.numerator} {params.reference} {params.narrow_summits} > {log:q} 2>&1
        """


rule motif_enrichment:
    input:
        peaks="results/peaks/consensus/{factor}.bed",
        genome=lambda wc: REF["genome"],
        genome_index=lambda wc: f"{REF['genome']}.fai",
        blacklist=lambda wc: REF["black_list"],
        script="analysis/motif_enrichment.R"
    output:
        table="results/motifs/{factor}/motif_enrichment.tsv",
        csv="results/motifs/{factor}/motif_summary.csv",
        pdf="results/motifs/{factor}/motif_summary.pdf"
    wildcard_constraints:
        factor="|".join(re.escape(s) for s in NARROW_FACTOR_SLUGS) if NARROW_FACTOR_SLUGS else "$^"
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.pdf)),
        factor=lambda wc: SLUG_TO_FACTOR[wc.factor],
        window=MOTIF_WINDOW_BP,
        background_multiplier=MOTIF_BACKGROUND_MULTIPLIER,
        seed=MOTIF_SEED
    threads: 1
    log:
        "results/logs/motif_enrichment_{factor}.log"
    conda:
        "../envs/r_analysis.yaml"
    shell:
        """
        mkdir -p {params.outdir:q} results/logs
        Rscript {input.script:q} {input.peaks:q} {input.genome:q} \
            {input.blacklist:q} {params.outdir:q} {params.window} \
            {params.background_multiplier} {params.seed} {params.factor:q} \
            > {log:q} 2>&1
        """
