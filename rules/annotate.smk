# Peak annotation and functional enrichment (tertiary analysis).


rule tss_bed:
    """TSS of protein-coding genes. Used as the anchor for promoter-mark profiles
    (H3K4me3, H3K9ac) — for those marks the biology is defined relative to the TSS,
    not to wherever the peak caller happened to place a summit."""
    input:
        gtf=REF["gtf"],
    output:
        bed="results/annotation/tss.bed",
    params:
        keep_re=config["alignment"]["keep_chroms_regex"],
    log:
        "results/logs/annotation/tss.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/annotation $(dirname {log})
        zcat -f {input.gtf} \
          | awk -v OFS='\t' '$3=="gene" && /gene_type "protein_coding"/ {{
                match($0, /gene_name "[^"]+"/);
                name = substr($0, RSTART+11, RLENGTH-12);
                if ($7=="+") {{ tss=$4-1 }} else {{ tss=$5-1 }}
                if (tss < 0) tss = 0
                print $1, tss, tss+1, name, 0, $7
            }}' \
          | awk -v OFS='\t' '$1 ~ /{params.keep_re}/' \
          | sort -k1,1 -k2,2n -u > {output.bed} 2> {log}
        echo "$(wc -l < {output.bed}) protein-coding TSS" >> {log}
        """


rule annotate_peaks:
    """ChIPseeker annotation of the consensus peaks plus GO enrichment on the genes
    near differential peaks.

    The peak-mode rule shows up here too. For a punctate mark the nearest-gene
    assignment is meaningful — a peak is a discrete element with a plausible target.
    For a broad domain spanning hundreds of kb, "nearest gene" is close to
    meaningless, so for broad marks the genes reported are every gene the domain
    overlaps rather than a single nearest neighbour.
    """
    input:
        consensus="results/peaks/consensus/{target}.bed",
        up="results/differential/{target}/up.bed",
        down="results/differential/{target}/down.bed",
        results="results/differential/{target}/results.tsv",
        tss="results/annotation/tss.bed",
    output:
        annotation="results/annotation/{target}/peak_annotation.tsv",
        genes="results/annotation/{target}/differential_genes.tsv",
        go="results/annotation/{target}/go_enrichment.tsv",
        distribution="results/annotation/{target}/feature_distribution.pdf",
    params:
        registry=config["mark_registry"],
        txdb=REF["txdb"],
        orgdb=REF["orgdb"],
    log:
        "results/logs/annotation/{target}.log",
    conda:
        "../envs/r_analysis.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/annotation/{wildcards.target} $(dirname {log})
        # Isolate the conda R from the host R installation. ~/.Rprofile is sourced on
        # every R startup and a .libPaths() call in it will prepend a library built
        # against a different R version, which then fails at dyn.load with an
        # undefined-symbol error. Clearing R_LIBS_USER alone does not help, because
        # the profile runs after the environment is read.
        export R_PROFILE_USER=/dev/null
        export R_ENVIRON_USER=/dev/null
        export R_LIBS_USER=""
        export R_LIBS_SITE=""
        Rscript analysis/annotate_peaks.R \
            --consensus {input.consensus} \
            --up {input.up} --down {input.down} \
            --target {wildcards.target} \
            --registry {params.registry} --tss {input.tss} \
            --txdb {params.txdb} --orgdb {params.orgdb} \
            --outdir results/annotation/{wildcards.target} > {log} 2>&1
        """
