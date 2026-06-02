#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def infer_factor(sample):
    return sample.get("factor") or sample.get("mark") or sample["id"].split("_")[0]


def build_diffbind_sheet(cfg):
    control_ids = {control["id"] for control in cfg.get("chip_controls", [])}
    rows = []
    for sample in cfg.get("chip_samples", []):
        control_id = sample.get("control")
        if not control_id:
            raise ValueError(f"ChIP sample {sample['id']} is missing a control")
        if control_id not in control_ids:
            raise ValueError(f"ChIP sample {sample['id']} references unknown control {control_id}")
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
                "Peaks": f"results/peaks/{sample['id']}_peaks.broadPeak",
                "PeakCaller": sample.get("peak_caller", "bed"),
            }
        )
    if not rows:
        raise ValueError("No chip_samples were found in the config")
    return pd.DataFrame(rows)


def write_sheet(cfg_path, diffbind_path):
    cfg = load_config(cfg_path)
    diffbind = build_diffbind_sheet(cfg)
    diffbind_path = Path(diffbind_path)
    diffbind_path.parent.mkdir(parents=True, exist_ok=True)
    diffbind.to_csv(diffbind_path, index=False)


def main():
    if "snakemake" in globals():
        write_sheet(snakemake.input.config, snakemake.output.diffbind)
        return

    parser = argparse.ArgumentParser(description="Build a DiffBind sample sheet from config.yaml")
    parser.add_argument("--config", default="config.yaml", help="YAML config file")
    parser.add_argument("--diffbind", default="results/diffbind/sample_sheet.csv", help="Output DiffBind sample sheet")
    args = parser.parse_args()
    write_sheet(args.config, args.diffbind)


if __name__ == "__main__":
    main()
