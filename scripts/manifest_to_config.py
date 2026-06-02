#!/usr/bin/env python3
import argparse
import csv
import yaml

REQUIRED_COLUMNS = ["sample_id", "type", "condition", "replicate", "fastq_r1", "fastq_r2"]
OPTIONAL_COLUMNS = ["control", "mark", "factor", "tissue", "peak_caller"]
SPECIES_META = {
    "human": {"name": "hg38", "gsize": "2.7e9"},
    "mouse": {"name": "mm10", "gsize": "1.87e9"},
    "rat": {"name": "rn6", "gsize": "2.53e9"},
}


def clean(value):
    return "" if value is None else str(value).strip()


def parse_manifest(path):
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"Missing manifest columns: {', '.join(missing)}")
        for row in reader:
            rows.append({k: clean(row.get(k)) for k in reader.fieldnames})
    return rows


def sample_from_row(row):
    sample = {
        "id": row["sample_id"],
        "fastq": [row["fastq_r1"], row["fastq_r2"]],
        "condition": row["condition"],
        "replicate": int(row["replicate"]),
    }
    factor = row.get("factor") or row.get("mark")
    if factor:
        sample["factor"] = factor
    for key in ["tissue", "peak_caller"]:
        if row.get(key):
            sample[key] = row[key]
    return sample


def build_config(rows, species, genome, chrom_sizes, bt2_index, black_list):
    chip_controls = []
    chip_samples = []

    for row in rows:
        sample_type = row["type"].lower()
        sample = sample_from_row(row)
        if sample_type == "control":
            chip_controls.append(sample)
        elif sample_type == "chip":
            control_id = row.get("control", "")
            if not control_id:
                raise ValueError(f"ChIP sample {row['sample_id']} is missing a control")
            sample["control"] = control_id
            chip_samples.append(sample)
        elif sample_type in {"rna", "rnaseq", "rna-seq"}:
            raise ValueError(
                f"RNA-seq sample {row['sample_id']} is not supported. "
                "This workflow is ChIP-seq differential binding only."
            )
        else:
            raise ValueError(f"Unknown type '{row['type']}' for sample {row['sample_id']}")

    return {
        "species": species,
        "references": {
            species: {
                "genome": genome,
                "chrom_sizes": chrom_sizes,
                "black_list": black_list,
                "bt2_index": bt2_index,
                "name": SPECIES_META[species]["name"],
                "gsize": SPECIES_META[species]["gsize"],
            }
        },
        "chip_controls": chip_controls,
        "chip_samples": chip_samples,
    }


def main():
    parser = argparse.ArgumentParser(description="Convert a ChIP-seq sample manifest TSV into Snakemake config.yaml")
    parser.add_argument("manifest", help="Sample manifest TSV")
    parser.add_argument("--species", choices=["human", "mouse", "rat"], required=True)
    parser.add_argument("--genome", required=True)
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
        args.chrom_sizes,
        args.bt2_index,
        args.black_list,
    )

    with open(args.output, "w") as fh:
        yaml.safe_dump(config, fh, sort_keys=False)

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
