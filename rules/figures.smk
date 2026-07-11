# Figures and the written summary.


rule figure_panels:
    input:
        gate="results/qc/qc_gate.tsv",
        cor="results/qc/correlation/spearman_matrix.tsv",
        peaks=all_chip_peaks(),
        diff=expand("results/differential/{t}/results.tsv", t=TARGETS),
        anno=expand("results/annotation/{t}/peak_annotation.tsv", t=TARGETS),
        go=expand("results/annotation/{t}/go_enrichment.tsv", t=TARGETS),
        motifs=expand("results/motifs/{t}/{d}/motif_enrichment.tsv",
                      t=MOTIF_TARGETS, d=["up", "down"]),
    output:
        done=touch("results/figures/figures.done"),
        fig1="results/figures/fig1_qc.pdf",
        fig2="results/figures/fig2_peak_landscape.pdf",
        fig4="results/figures/fig4_differential.pdf",
    params:
        registry=config["mark_registry"],
        samples=config["samples"],
    log:
        "results/logs/figures.log",
    conda:
        "../envs/r_analysis.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/figures $(dirname {log})
        # Isolate the conda R from the host R installation. ~/.Rprofile is sourced on
        # every R startup and a .libPaths() call in it will prepend a library built
        # against a different R version, which then fails at dyn.load with an
        # undefined-symbol error. Clearing R_LIBS_USER alone does not help, because
        # the profile runs after the environment is read.
        export R_PROFILE_USER=/dev/null
        export R_ENVIRON_USER=/dev/null
        export R_LIBS_USER=""
        export R_LIBS_SITE=""
        Rscript analysis/figures.R \
            --results results \
            --registry {params.registry} \
            --samples {params.samples} \
            --config config.yaml \
            --outdir results/figures > {log} 2>&1
        """


rule browser_tracks:
    """Genome-browser snapshots at the loci that decide whether the experiment
    worked. For TCR stimulation these are the immediate-early and activation genes:
    if H3K27ac does not increase over IL2, CD69, EGR2 and NR4A1, either the
    stimulation or the ChIP failed, and nothing else in the analysis is worth
    reading."""
    input:
        bws=expand("results/bigwig/{s}.log2ratio.bw", s=CHIP_SAMPLES),
        peaks=expand("results/peaks/consensus/{t}.bed", t=TARGETS),
        gtf=REF["gtf"],
    output:
        done=touch("results/figures/browser/tracks.done"),
    params:
        loci=lambda wc: " ".join(
            f"{name}:{v['chrom']}:{v['start']}-{v['end']}"
            for name, v in config["browser_loci"].items()
        ),
        samples=config["samples"],
        registry=config["mark_registry"],
    log:
        "results/logs/browser_tracks.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/figures/browser $(dirname {log})
        python3 scripts/browser_tracks.py \
            --samples {params.samples} \
            --registry {params.registry} \
            --gtf {input.gtf} \
            --outdir results/figures/browser \
            --loci {params.loci} > {log} 2>&1
        """


rule analysis_summary:
    input:
        gate="results/qc/qc_gate.tsv",
        gate_md="results/qc/qc_gate.md",
        diff=expand("results/differential/{t}/results.tsv", t=TARGETS),
        sf=expand("results/differential/{t}/size_factors.tsv", t=TARGETS),
        consensus=expand("results/peaks/consensus/{t}.bed", t=TARGETS),
        go=expand("results/annotation/{t}/go_enrichment.tsv", t=TARGETS),
        motifs=expand("results/motifs/{t}/{d}/motif_enrichment.tsv",
                      t=MOTIF_TARGETS, d=["up", "down"]),
        figures="results/figures/figures.done",
        browser="results/figures/browser/tracks.done",
    output:
        md="results/summary/analysis_summary.md",
    params:
        registry=config["mark_registry"],
        samples=config["samples"],
    log:
        "results/logs/summary.log",
    conda:
        "../envs/chipseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/summary $(dirname {log})
        python3 scripts/analysis_summary.py \
            --results results --config config.yaml \
            --registry {params.registry} --samples {params.samples} \
            --out {output.md} > {log} 2>&1
        """
