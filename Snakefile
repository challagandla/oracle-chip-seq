import re
import sys
from pathlib import Path

import pandas as pd
from snakemake.exceptions import WorkflowError

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from marks import MarkRegistry  # noqa: E402

configfile: "config.yaml"

REG = MarkRegistry(config["mark_registry"])
REF = config["references"][config["species"]]

# ---------------------------------------------------------------- sample table

samples = pd.read_csv(config["samples"], sep="\t", dtype=str).fillna("")
samples["replicate"] = samples["replicate"].astype(int)
samples = samples.set_index("sample_id", drop=False)

CHIP = samples[samples.assay == "chip"]
INPUTS = samples[samples.assay == "input"]
ALL_SAMPLES = list(samples.sample_id)
CHIP_SAMPLES = list(CHIP.sample_id)
TARGETS = sorted(CHIP.target.unique())

# Layout is a per-sample property here, not a global one: this study sequenced
# replicate 1 single-end and replicate 2 paired-end.
SE_SAMPLES = list(samples[samples.layout == "single"].sample_id)
PE_SAMPLES = list(samples[samples.layout == "paired"].sample_id)

CONTRAST = config["contrast"]
REF_LEVEL, TRT_LEVEL = CONTRAST["reference"], CONTRAST["treatment"]
CONDITIONS = [REF_LEVEL, TRT_LEVEL]

NARROW_TARGETS = [t for t in TARGETS if not REG.is_broad(t)]
BROAD_TARGETS = [t for t in TARGETS if REG.is_broad(t)]
IDR_TARGETS = [t for t in TARGETS if REG.uses_idr(t)]
OVERLAP_TARGETS = [t for t in TARGETS if not REG.uses_idr(t)]
MOTIF_TARGETS = [t for t in TARGETS if REG.motifs_enabled(t)]


def _validate_design():
    """Fail at parse time rather than three hours into an alignment."""
    problems = []
    unknown = sorted({t for t in TARGETS if not REG.is_known(t)})
    if unknown:
        # Not fatal — defaults apply. But an unrecognised histone mark silently
        # treated as a narrow TF is exactly the mistake the registry exists to
        # prevent, so it has to be visible.
        print(
            "[mark-registry] WARNING: not in registry, defaulting to "
            f"{REG.defaults['peak_mode']}: {', '.join(unknown)}",
            file=sys.stderr,
        )
    for sid, row in CHIP.iterrows():
        if not row.control_id:
            problems.append(f"{sid}: no control_id")
        elif row.control_id not in samples.index:
            problems.append(f"{sid}: control '{row.control_id}' is not in the sample table")
    for t in TARGETS:
        for cond in CONDITIONS:
            n = len(CHIP[(CHIP.target == t) & (CHIP.condition == cond)])
            if n < 2:
                problems.append(
                    f"{t}/{cond}: {n} replicate(s); differential binding and "
                    "reproducibility filtering both need at least 2"
                )
    for sid, row in samples.iterrows():
        if row.layout not in ("single", "paired"):
            problems.append(f"{sid}: layout must be 'single' or 'paired', got '{row.layout}'")
    if problems:
        raise WorkflowError("Invalid design:\n  - " + "\n  - ".join(problems))


_validate_design()

# ------------------------------------------------------------------- accessors

def is_paired(sample):
    return samples.loc[sample, "layout"] == "paired"


def target_of(sample):
    return samples.loc[sample, "target"]


def control_of(sample):
    return samples.loc[sample, "control_id"]


def raw_fastqs(sample):
    if is_paired(sample):
        return [f"data/raw/{sample}_R1.fastq.gz", f"data/raw/{sample}_R2.fastq.gz"]
    return [f"data/raw/{sample}.fastq.gz"]


def reads_of(sample):
    return ["R1", "R2"] if is_paired(sample) else ["R1"]


def peak_ext(sample):
    """The peak-mode rule applied at the filename level, so a broad target cannot
    silently be consumed as a narrowPeak anywhere downstream."""
    return REG.peak_ext(target_of(sample))


def sample_peaks(sample):
    return f"results/peaks/{sample}_peaks.{peak_ext(sample)}"


def reps_for(target, condition):
    sub = CHIP[(CHIP.target == target) & (CHIP.condition == condition)]
    return list(sub.sort_values("replicate").sample_id)


def all_chip_peaks():
    return [sample_peaks(s) for s in CHIP_SAMPLES]


_SE_RE = "|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
_PE_RE = "|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"


wildcard_constraints:
    sample="|".join(re.escape(s) for s in ALL_SAMPLES),
    target="|".join(re.escape(t) for t in TARGETS) if TARGETS else "$^",
    condition="|".join(re.escape(c) for c in CONDITIONS),
    read="R1|R2",
    ext="narrowPeak|broadPeak",
    direction="up|down",


include: "rules/sra.smk"
include: "rules/align.smk"
include: "rules/qc.smk"
include: "rules/peaks.smk"
include: "rules/signal.smk"
include: "rules/differential.smk"
include: "rules/annotate.smk"
include: "rules/motifs.smk"
include: "rules/figures.smk"


rule all:
    input:
        # QC
        "results/qc/qc_gate.tsv",
        "results/qc/multiqc/multiqc_report.html",
        # Signal
        expand("results/bigwig/{sample}.cpm.bw", sample=ALL_SAMPLES),
        expand("results/bigwig/{sample}.log2ratio.bw", sample=CHIP_SAMPLES),
        "results/qc/correlation/spearman_heatmap.pdf",
        "results/qc/correlation/pca.pdf",
        "results/qc/fingerprint/fingerprint.pdf",
        # Peaks + reproducibility
        all_chip_peaks(),
        expand("results/peaks/consensus/{target}.bed", target=TARGETS),
        # Profiles — geometry chosen per peak mode
        expand("results/profiles/{target}_heatmap.pdf", target=TARGETS),
        # Differential binding
        expand("results/differential/{target}/results.tsv", target=TARGETS),
        # Annotation + enrichment
        expand("results/annotation/{target}/peak_annotation.tsv", target=TARGETS),
        # Motifs — narrow/TF targets only
        expand("results/motifs/{target}/{direction}/motif_enrichment.tsv",
               target=MOTIF_TARGETS, direction=["up", "down"]),
        # Figures
        "results/figures/figures.done",
        "results/summary/analysis_summary.md",
