#!/usr/bin/env python3
import argparse
from collections import Counter
from pathlib import Path

import pandas as pd
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def infer_factor(sample):
    return sample.get("factor") or sample.get("mark") or sample["id"].split("_")[0]


def peak_ext(sample):
    # histone -> MACS2 broad peaks, tf -> narrow peaks (must match the Snakefile).
    return "broadPeak" if sample.get("mark_type", "histone") == "histone" else "narrowPeak"


def validate_design(cfg):
    controls = cfg.get("chip_controls", [])
    samples = cfg.get("chip_samples", [])
    if not controls:
        raise ValueError("No chip_controls were found in the config")
    if not samples:
        raise ValueError("No chip_samples were found in the config")

    ids = [item["id"] for item in controls + samples]
    duplicate_ids = sorted([sample_id for sample_id, count in Counter(ids).items() if count > 1])
    if duplicate_ids:
        raise ValueError(f"Duplicate sample IDs found: {', '.join(duplicate_ids)}")

    control_ids = {control["id"] for control in controls}
    missing_controls = sorted(
        {
            sample.get("control", "")
            for sample in samples
            if not sample.get("control") or sample.get("control") not in control_ids
        }
    )
    if missing_controls:
        raise ValueError(f"Unknown or missing ChIP control IDs: {', '.join(missing_controls)}")

    condition_counts = Counter(sample["condition"] for sample in samples)
    if len(condition_counts) < 2:
        raise ValueError("Differential binding requires at least two ChIP-seq conditions")
    low_rep_conditions = sorted(condition for condition, count in condition_counts.items() if count < 2)
    if low_rep_conditions:
        raise ValueError(
            "Differential binding requires at least two ChIP-seq replicates per condition; "
            f"insufficient replicates for: {', '.join(low_rep_conditions)}"
        )


def build_diffbind_sheet(cfg):
    validate_design(cfg)
    rows = []
    for sample in cfg.get("chip_samples", []):
        control_id = sample["control"]
        rows.append(
            {
                "SampleID": sample["id"],
                "Tissue": sample.get("tissue", "example"),
                "Factor": infer_factor(sample),
                "Condition": sample["condition"],
                "Replicate": sample["replicate"],
                "bamReads": f"results/bam/{sample['id']}.sorted.bam",
                "ControlID": control_id,
                "bamControl": f"results/bam/{control_id}.sorted.bam",
                "Peaks": f"results/peaks/{sample['id']}_peaks.{peak_ext(sample)}",
                "PeakCaller": sample.get("peak_caller", "bed"),
            }
        )
    return pd.DataFrame(rows)


def write_sheet(cfg_path, diffbind_path):
    cfg = load_config(cfg_path)
    diffbind = build_diffbind_sheet(cfg)
    diffbind_path = Path(diffbind_path)
    diffbind_path.parent.mkdir(parents=True, exist_ok=True)
    diffbind.to_csv(diffbind_path, index=False)


def main():
    parser = argparse.ArgumentParser(description="Build a DiffBind sample sheet from config.yaml")
    parser.add_argument("--config", default="config.yaml", help="YAML config file")
    parser.add_argument("--diffbind", default="results/diffbind/sample_sheet.csv", help="Output DiffBind sample sheet")
    args = parser.parse_args()
    write_sheet(args.config, args.diffbind)


if __name__ == "__main__":
    main()
