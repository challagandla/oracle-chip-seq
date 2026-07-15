#!/usr/bin/env python3
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import re

import pandas as pd
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def infer_factor(sample):
    return str(sample.get("factor") or sample.get("mark") or "").strip()


def peak_ext(sample):
    mode = str(sample.get("peak_mode", "")).strip().lower()
    if mode not in {"narrow", "broad"}:
        raise ValueError(
            f"Sample {sample.get('id', '<unnamed>')} must set peak_mode to 'narrow' or 'broad'"
        )
    return "broadPeak" if mode == "broad" else "narrowPeak"


def validate_design(cfg):
    if not isinstance(cfg, dict):
        raise ValueError("Configuration must be a YAML mapping")
    controls = cfg.get("chip_controls", [])
    samples = cfg.get("chip_samples", [])
    if not isinstance(controls, list) or not all(isinstance(item, dict) for item in controls):
        raise ValueError("chip_controls must be a list of sample mappings")
    if not isinstance(samples, list) or not all(isinstance(item, dict) for item in samples):
        raise ValueError("chip_samples must be a list of sample mappings")
    if not controls:
        raise ValueError("No chip_controls were found in the config")
    if not samples:
        raise ValueError("No chip_samples were found in the config")

    differential = cfg.get("differential_binding", {})
    if not isinstance(differential, dict):
        raise ValueError("differential_binding must be a mapping")
    numerator = str(differential.get("numerator_condition", "treated")).strip()
    reference = str(differential.get("reference_condition", "control")).strip()
    narrow_summits = differential.get("narrow_summits", 200)
    if not numerator or not reference or numerator == reference:
        raise ValueError(
            "DiffBind numerator and reference conditions must be distinct and non-empty"
        )
    if (
        isinstance(narrow_summits, bool)
        or not isinstance(narrow_summits, int)
        or narrow_summits <= 0
    ):
        raise ValueError("differential_binding.narrow_summits must be a positive integer")

    safe_id = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
    fastq_owners = defaultdict(list)
    for item in controls + samples:
        raw_id = item.get("id", "")
        sample_id = str(raw_id).strip()
        if not isinstance(raw_id, str):
            raise ValueError(f"Sample ID {raw_id!r} must be text")
        if raw_id != sample_id:
            raise ValueError(f"Sample ID {raw_id!r} has leading or trailing whitespace")
        if not safe_id.fullmatch(sample_id):
            raise ValueError(
                f"Sample ID {sample_id!r} is unsafe; use letters, numbers, '.', '_' and '-'"
            )
        fastq = item.get("fastq")
        if (
            not isinstance(fastq, list)
            or len(fastq) != 2
            or not all(isinstance(path, str) and path.strip() for path in fastq)
        ):
            raise ValueError(f"Sample {sample_id} requires exactly two FASTQ paths")
        if not all(path.lower().endswith((".fastq.gz", ".fq.gz")) for path in fastq):
            raise ValueError(f"Sample {sample_id} FASTQs must end in .fastq.gz or .fq.gz")
        normalized = [str(Path(path).resolve(strict=False)) for path in fastq]
        if normalized[0] == normalized[1]:
            raise ValueError(f"Sample {sample_id} uses the same file for R1 and R2")
        for read, path in zip(("R1", "R2"), normalized):
            fastq_owners[path].append(f"{sample_id} {read}")

    reused_fastqs = {path: owners for path, owners in fastq_owners.items() if len(owners) > 1}
    if reused_fastqs:
        details = "; ".join(
            f"{path}: {', '.join(owners)}" for path, owners in sorted(reused_fastqs.items())
        )
        raise ValueError(f"FASTQ files are assigned more than once: {details}")

    ids = [item["id"] for item in controls + samples]
    duplicate_ids = sorted([sample_id for sample_id, count in Counter(ids).items() if count > 1])
    if duplicate_ids:
        raise ValueError(f"Duplicate sample IDs found: {', '.join(duplicate_ids)}")

    control_ids = {control["id"] for control in controls}
    control_conditions = {}
    for control in controls:
        raw_condition = control.get("condition")
        condition = raw_condition.strip() if isinstance(raw_condition, str) else ""
        if not condition:
            raise ValueError(f"Control {control['id']} is missing condition")
        control_conditions[control["id"]] = condition
    missing_controls = sorted(
        {
            sample.get("control", "")
            for sample in samples
            if not sample.get("control") or sample.get("control") not in control_ids
        }
    )
    if missing_controls:
        raise ValueError(f"Unknown or missing ChIP control IDs: {', '.join(missing_controls)}")

    factor_condition_counts = Counter()
    replicates_by_factor_condition = defaultdict(list)
    modes_by_factor = defaultdict(set)
    tissues_by_factor = defaultdict(set)
    for sample in samples:
        raw_factor = sample.get("factor") or sample.get("mark")
        factor = infer_factor(sample)
        if not isinstance(raw_factor, str) or not factor:
            raise ValueError(f"Sample {sample['id']} is missing factor")
        raw_condition = sample.get("condition", "")
        raw_tissue = sample.get("tissue")
        condition = str(raw_condition).strip()
        tissue = raw_tissue.strip() if isinstance(raw_tissue, str) else ""
        if not isinstance(raw_condition, str) or not condition:
            raise ValueError(f"Sample {sample['id']} is missing condition")
        control_id = sample.get("control")
        if control_id in control_conditions and control_conditions[control_id] != condition:
            raise ValueError(
                f"Sample {sample['id']} condition '{condition}' does not match "
                f"control {control_id} condition '{control_conditions[control_id]}'"
            )
        if not tissue:
            raise ValueError(f"Sample {sample['id']} is missing tissue")
        replicate = sample.get("replicate", 0)
        if isinstance(replicate, bool) or not isinstance(replicate, int) or replicate < 1:
            raise ValueError(f"Sample {sample['id']} replicate must be a positive integer")
        mode = str(sample.get("peak_mode", "")).strip().lower()
        if mode not in {"narrow", "broad"}:
            raise ValueError(f"Sample {sample['id']} must set peak_mode to 'narrow' or 'broad'")
        factor_condition_counts[(factor, condition)] += 1
        replicates_by_factor_condition[(factor, condition)].append(replicate)
        modes_by_factor[factor].add(mode)
        tissues_by_factor[factor].add(tissue)

    for factor, modes in modes_by_factor.items():
        if len(modes) != 1:
            raise ValueError(f"Factor {factor} mixes peak modes: {', '.join(sorted(modes))}")
        if len(tissues_by_factor[factor]) != 1:
            raise ValueError(
                f"Factor {factor} requires one tissue for the supplied simple DiffBind model"
            )
        conditions = sorted(condition for f, condition in factor_condition_counts if f == factor)
        if len(conditions) != 2:
            raise ValueError(
                f"Factor {factor} requires exactly two conditions for one DiffBind contrast"
            )
        if set(conditions) != {reference, numerator}:
            raise ValueError(
                f"Factor {factor} conditions must match DiffBind reference '{reference}' "
                f"and numerator '{numerator}'"
            )
        low_rep_conditions = [
            condition
            for condition in conditions
            if factor_condition_counts[(factor, condition)] < 2
        ]
        if low_rep_conditions:
            raise ValueError(
                f"Factor {factor} requires at least two replicates per condition; "
                f"insufficient replicates for: {', '.join(low_rep_conditions)}"
            )
        duplicate_rep_conditions = [
            condition
            for condition in conditions
            if len(set(replicates_by_factor_condition[(factor, condition)]))
            != len(replicates_by_factor_condition[(factor, condition)])
        ]
        if duplicate_rep_conditions:
            raise ValueError(
                f"Factor {factor} has duplicate replicate numbers in: "
                f"{', '.join(duplicate_rep_conditions)}"
            )


def build_diffbind_sheet(cfg):
    validate_design(cfg)
    rows = []
    for sample in cfg.get("chip_samples", []):
        control_id = sample["control"]
        rows.append(
            {
                "SampleID": sample["id"],
                "Tissue": str(sample["tissue"]).strip(),
                "Factor": infer_factor(sample),
                "Condition": str(sample["condition"]).strip(),
                "Replicate": sample["replicate"],
                "bamReads": f"results/bam/{sample['id']}.sorted.bam",
                "ControlID": control_id,
                "bamControl": f"results/bam/{control_id}.sorted.bam",
                "Peaks": f"results/peaks/{sample['id']}_peaks.{peak_ext(sample)}",
                "PeakCaller": "bed",
                "PeakMode": str(sample["peak_mode"]).strip().lower(),
            }
        )
    return pd.DataFrame(rows)


def write_sheet_from_config(cfg, diffbind_path, resolved_config_path=None):
    diffbind = build_diffbind_sheet(cfg)
    diffbind_path = Path(diffbind_path)
    diffbind_path.parent.mkdir(parents=True, exist_ok=True)
    diffbind.to_csv(diffbind_path, index=False)
    if resolved_config_path is not None:
        resolved_config_path = Path(resolved_config_path)
        resolved_config_path.parent.mkdir(parents=True, exist_ok=True)
        with open(resolved_config_path, "w") as fh:
            yaml.safe_dump(dict(cfg), fh, sort_keys=False)


def write_sheet(cfg_path, diffbind_path, resolved_config_path=None):
    write_sheet_from_config(load_config(cfg_path), diffbind_path, resolved_config_path)


def main():
    parser = argparse.ArgumentParser(description="Build a DiffBind sample sheet from config.yaml")
    parser.add_argument("--config", default="config.yaml", help="YAML config file")
    parser.add_argument(
        "--diffbind",
        default="results/diffbind/sample_sheet.csv",
        help="Output DiffBind sample sheet",
    )
    parser.add_argument(
        "--resolved-config",
        help="Optional path for a normalized copy of the configuration used",
    )
    args = parser.parse_args()
    write_sheet(args.config, args.diffbind, args.resolved_config)


if "snakemake" in globals():
    snakemake_context = globals()["snakemake"]
    write_sheet_from_config(
        snakemake_context.config,
        snakemake_context.output.diffbind,
        snakemake_context.output.resolved_config,
    )
elif __name__ == "__main__":
    main()
