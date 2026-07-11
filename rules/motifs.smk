# Motif enrichment (monaLisa + JASPAR).
#
# Runs only for targets whose registry entry sets motifs.enabled. That flag is not
# a convenience switch — it states whether the question is coherent for the mark:
#
#   TF (CTCF)         The motif IS the validation. If the top motif in CTCF peaks
#                     is not the CTCF motif, the experiment failed, whatever the
#                     FRiP says.
#   H3K27ac/H3K4me3   Punctate enough that the summit window contains the binding
#                     sites of the factors that established the element. Read as
#                     "which factors drive the responsive enhancers", not as a
#                     validation of the antibody.
#   H3K27me3, broad   Disabled. Scanning a multi-kb Polycomb domain for 8-mers
#                     returns its base composition, and the software will report
#                     that with confident p-values — which is exactly why the flag
#                     exists rather than just letting it run.
#
# Up and down peaks are tested separately. Pooling them asks "what motifs are in
# peaks that changed", which mixes gain with loss and usually recovers neither.
#
# HOMER is deliberately not used: it is academic-use-only and not redistributable.


rule motif_enrichment:
    input:
        peaks="results/differential/{target}/{direction}.bed",
        background="results/peaks/consensus/{target}.bed",
        genome=REF["genome"],
    output:
        tsv="results/motifs/{target}/{direction}/motif_enrichment.tsv",
    wildcard_constraints:
        target="|".join(re.escape(t) for t in MOTIF_TARGETS) if MOTIF_TARGETS else "$^",
    params:
        outdir=lambda wc: f"results/motifs/{wc.target}/{wc.direction}",
        registry=config["mark_registry"],
    threads: 4
    log:
        "results/logs/motifs/{target}_{direction}.log",
    conda:
        "../envs/r_analysis.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        Rscript analysis/motif_enrichment.R \
            --peaks {input.peaks} \
            --background {input.background} \
            --genome {input.genome} \
            --target {wildcards.target} \
            --direction {wildcards.direction} \
            --registry {params.registry} \
            --outdir {params.outdir} > {log} 2>&1
        """
