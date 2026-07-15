import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from shlex import quote

from snakemake.exceptions import WorkflowError

configfile: "config.yaml"

include: "workflow/common.smk"


rule all:
    input:
        RAW_FASTQC_HTML,
        TRIMMED_FASTQC_HTML,
        FASTQ_SCREEN_TEXT,
        FASTQ_SCREEN_HTML,
        expand("results/bam/{sample}.sorted.bam.bai", sample=ALL_SAMPLES),
        expand("results/qc/samtools/{sample}.flagstat.txt", sample=ALL_SAMPLES),
        expand("results/qc/samtools/{sample}.markdup.txt", sample=ALL_SAMPLES),
        PEAKS,
        RPGC_BIGWIGS,
        LOG2_BIGWIGS,
        FACTOR_CONSENSUS,
        "results/peaks/consensus_peaks.bed",
        DEEPTOOLS_OUTPUTS,
        DIFFBIND_OUTPUTS,
        MOTIF_OUTPUTS,
        "results/provenance/resolved_config.yaml",
        MULTIQC_REPORT


include: "workflow/qc.smk"
include: "workflow/alignment.smk"
include: "workflow/peaks.smk"
include: "workflow/signal.smk"
include: "workflow/deeptools.smk"
include: "workflow/analysis.smk"
include: "workflow/reporting.smk"
