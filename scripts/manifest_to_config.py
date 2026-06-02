#!/usr/bin/env python3
import argparse
import csv
import yaml
from pathlib import Path

REQUIRED_COLUMNS = ["sample_id", "type", "condition", "replicate", "fastq_r1", "fastq_r2"]


def parse_manifest(path):
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"Missing manifest columns: {', '.join(missing)}")
        for row in reader:
            rows.append({k: row[k].strip() for k in reader.fieldnames})
    return rows


def build_config(rows, species, genome, annotation, transcriptome, chrom_sizes, bt2_index, black_list):
    chip_controls = []
    chip_samples = []
    rna_samples = []

    for row in rows:
        sample = {
            "id": row["sample_id"],
            "fastq": [row["fastq_r1"], row["fastq_r2"]],
            "condition": row["condition"],
            "replicate": int(row["replicate"]),
        }
        sample_type = row["type"].lower()
        if sample_type == "control":
            chip_controls.append(sample)
        elif sample_type == "chip":
            sample["control"] = row.get("control", "") or ""
            chip_samples.append(sample)
        elif sample_type == "rna":
            rna_samples.append(sample)
        else:
            raise ValueError(f"Unknown type '{row['type']}' for sample {row['sample_id']}")

    return {
        "species": species,
        "references": {
            species: {
                "genome": genome,
                "annotation": annotation,
                "transcriptome": transcriptome,
                "chrom_sizes": chrom_sizes,
                "black_list": black_list,
                "bt2_index": bt2_index,
                "name": species,
                "gsize": "2.7e9" if species == "human" else "1.87e9" if species == "mouse" else "2.53e9",
            }
        },
        "chip_controls": chip_controls,
        "chip_samples": chip_samples,
        "rna_samples": rna_samples,
    }


def main():
    parser = argparse.ArgumentParser(description="Convert a sample manifest TSV into Snakemake config.yaml")
    parser.add_argument("manifest", help="Sample manifest TSV")
    parser.add_argument("--species", choices=["human", "mouse", "rat"], required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--transcriptome", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--bt2-index", required=True)
    parser.add_argument("--black-list", required=True)
    parser.add_argument("--output", default="config.yaml")
    args = parser.parse_args()

    rows = parse_manifest(args.manifest)
    config = build_config(
        rows,
        args.species,
        args.genome,
        args.annotation,
        args.transcriptome,
        args.chrom_sizes,
        args.bt2_index,
        args.black_list,
    )

    with open(args.output, "w") as fh:
        yaml.safe_dump(config, fh, sort_keys=False)

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
