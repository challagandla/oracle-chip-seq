import argparse
import yaml
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(description="Build sample sheets for DiffBind and RNA-seq")
parser.add_argument("--config", default="config.yaml", help="YAML config file")
parser.add_argument("--diffbind", default="results/diffbind/sample_sheet.csv", help="Output DiffBind sample sheet")
parser.add_argument("--rna", default="results/rnaseq/sample_metadata.tsv", help="Output RNA sample metadata")
args = parser.parse_args()

with open(args.config) as fh:
    cfg = yaml.safe_load(fh)

rows = []
control_map = {c["id"]: c for c in cfg["chip_controls"]}
for sample in cfg["chip_samples"]:
    control_id = sample["control"]
    rows.append({
        "SampleID": sample["id"],
        "Condition": sample["condition"],
        "Replicate": sample["replicate"],
        "bamReads": f"results/bam/{sample['id']}.sorted.bam",
        "bamControl": f"results/bam/{control_id}.sorted.bam",
        "Peaks": f"results/peaks/{sample['id']}_peaks.broadPeak",
        "ControlID": control_id,
    })

pd.DataFrame(rows).to_csv(args.diffbind, index=False)

rna_rows = []
for sample in cfg["rna_samples"]:
    rna_rows.append({
        "SampleID": sample["id"],
        "Condition": sample["condition"],
        "Replicate": sample["replicate"],
        "QuantDir": f"results/rnaseq/salmon/{sample['id']}",
    })

pd.DataFrame(rna_rows).to_csv(args.rna, sep="\t", index=False)
