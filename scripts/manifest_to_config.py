#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import yaml

REQUIRED_COLUMNS = [
    "sample_id",
    "type",
    "condition",
    "replicate",
    "fastq_r1",
    "fastq_r2",
    "control",
    "factor",
    "tissue",
    "peak_mode",
]
SPECIES_META = {
    "human": {"name": "hg38", "gsize": "2.7e9"},
    "mouse": {"name": "mm39", "gsize": "1.87e9"},
    "rat": {"name": "rn6", "gsize": "2.53e9"},
}


def clean(value):
    return "" if value is None else str(value).strip()


def parse_manifest(path):
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Manifest is empty")
        missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"Missing manifest columns: {', '.join(missing)}")
        for row in reader:
            rows.append({k: clean(row.get(k)) for k in reader.fieldnames})
    if not rows:
        raise ValueError("Manifest contains no sample rows")
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
    for key in ["tissue", "peak_mode"]:
        if row.get(key):
            sample[key] = row[key]
    return sample


def validate_design(chip_controls, chip_samples):
    if not chip_controls:
        raise ValueError("Manifest must contain at least one ChIP input/control sample")
    if not chip_samples:
        raise ValueError("Manifest must contain at least one ChIP sample")

    fastq_owners = defaultdict(list)
    for sample in chip_controls + chip_samples:
        fastq = sample.get("fastq")
        if (
            not isinstance(fastq, list)
            or len(fastq) != 2
            or not all(isinstance(path, str) and path for path in fastq)
        ):
            raise ValueError(f"Sample {sample['id']} requires exactly two FASTQ paths")
        if not all(path.lower().endswith((".fastq.gz", ".fq.gz")) for path in fastq):
            raise ValueError(f"Sample {sample['id']} FASTQs must end in .fastq.gz or .fq.gz")
        normalized = [str(Path(path).resolve(strict=False)) for path in fastq]
        if normalized[0] == normalized[1]:
            raise ValueError(f"Sample {sample['id']} uses the same file for R1 and R2")
        for read, path in zip(("R1", "R2"), normalized):
            fastq_owners[path].append(f"{sample['id']} {read}")

    reused_fastqs = {path: owners for path, owners in fastq_owners.items() if len(owners) > 1}
    if reused_fastqs:
        details = "; ".join(
            f"{path}: {', '.join(owners)}" for path, owners in sorted(reused_fastqs.items())
        )
        raise ValueError(f"FASTQ files are assigned more than once: {details}")

    ids = [item["id"] for item in chip_controls + chip_samples]
    duplicate_ids = sorted([sample_id for sample_id, count in Counter(ids).items() if count > 1])
    if duplicate_ids:
        raise ValueError(f"Duplicate sample IDs found: {', '.join(duplicate_ids)}")

    control_ids = {control["id"] for control in chip_controls}
    control_conditions = {}
    for control in chip_controls:
        raw_condition = control.get("condition")
        condition = raw_condition.strip() if isinstance(raw_condition, str) else ""
        if not condition:
            raise ValueError(f"Control {control['id']} is missing condition")
        control_conditions[control["id"]] = condition
    missing_controls = sorted(
        {
            sample.get("control", "")
            for sample in chip_samples
            if not sample.get("control") or sample.get("control") not in control_ids
        }
    )
    if missing_controls:
        raise ValueError(f"Unknown or missing ChIP control IDs: {', '.join(missing_controls)}")

    factor_condition_counts = Counter()
    replicates_by_factor_condition = defaultdict(list)
    modes_by_factor = defaultdict(set)
    tissues_by_factor = defaultdict(set)
    for sample in chip_samples:
        factor = sample.get("factor", "")
        condition = sample.get("condition", "")
        mode = str(sample.get("peak_mode", "")).lower()
        raw_tissue = sample.get("tissue")
        tissue = raw_tissue.strip() if isinstance(raw_tissue, str) else ""
        if not factor:
            raise ValueError(f"ChIP sample {sample['id']} is missing factor")
        if not condition:
            raise ValueError(f"ChIP sample {sample['id']} is missing condition")
        control_id = sample.get("control")
        if control_id in control_conditions and control_conditions[control_id] != condition:
            raise ValueError(
                f"ChIP sample {sample['id']} condition '{condition}' does not match "
                f"control {control_id} condition '{control_conditions[control_id]}'"
            )
        if not tissue:
            raise ValueError(f"ChIP sample {sample['id']} is missing tissue")
        replicate = sample.get("replicate")
        if isinstance(replicate, bool) or not isinstance(replicate, int) or replicate < 1:
            raise ValueError(f"ChIP sample {sample['id']} replicate must be a positive integer")
        if mode not in {"narrow", "broad"}:
            raise ValueError(f"ChIP sample {sample['id']} peak_mode must be 'narrow' or 'broad'")
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


def build_config(
    rows,
    species,
    genome,
    chrom_sizes,
    bt2_index,
    black_list,
    effective_genome_size,
    fastq_screen_conf="config/fastq_screen.conf",
    numerator_condition="treated",
    reference_condition="control",
):
    chip_controls = []
    chip_samples = []

    for row in rows:
        sample_type = row["type"].lower()
        if sample_type in {"rna", "rnaseq", "rna-seq"}:
            raise ValueError(
                f"RNA-seq sample {row['sample_id']} is not supported. "
                "This workflow is ChIP-seq differential binding only."
            )

        sample = sample_from_row(row)
        if sample_type == "control":
            chip_controls.append(sample)
        elif sample_type == "chip":
            control_id = row.get("control", "")
            if not control_id:
                raise ValueError(f"ChIP sample {row['sample_id']} is missing a control")
            sample["control"] = control_id
            chip_samples.append(sample)
        else:
            raise ValueError(f"Unknown type '{row['type']}' for sample {row['sample_id']}")

    validate_design(chip_controls, chip_samples)
    if (
        isinstance(effective_genome_size, bool)
        or not isinstance(effective_genome_size, int)
        or effective_genome_size <= 0
    ):
        raise ValueError("Effective genome size must be a positive integer")
    numerator_condition = clean(numerator_condition)
    reference_condition = clean(reference_condition)
    if (
        not numerator_condition
        or not reference_condition
        or numerator_condition == reference_condition
    ):
        raise ValueError("Numerator and reference conditions must be distinct and non-empty")
    expected_conditions = {reference_condition, numerator_condition}
    for factor in sorted({sample["factor"] for sample in chip_samples}):
        observed = {sample["condition"] for sample in chip_samples if sample["factor"] == factor}
        if observed != expected_conditions:
            raise ValueError(
                f"Factor {factor} conditions {sorted(observed)} do not match reference "
                f"'{reference_condition}' and numerator '{numerator_condition}'"
            )

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
                "effective_genome_size": effective_genome_size,
            }
        },
        "contamination": {"fastq_screen_conf": fastq_screen_conf, "subset": 100000},
        "alignment": {"min_mapq": 30, "max_insert_size": 2000},
        "peak_calling": {
            "narrow_qvalue": 0.01,
            "broad_cutoff": 0.1,
            "consensus_min_replicates": 2,
        },
        "differential_binding": {
            "numerator_condition": numerator_condition,
            "reference_condition": reference_condition,
            "narrow_summits": 200,
        },
        "motif_enrichment": {
            "window_bp": 200,
            "background_multiplier": 2,
            "seed": 1,
        },
        "deeptools": {
            "threads": 4,
            "track_bin_size": 25,
            "matrix_bin_size": 50,
            "reference_point_upstream": 3000,
            "reference_point_downstream": 3000,
            "scale_regions_upstream": 3000,
            "scale_regions_downstream": 3000,
            "scale_regions_body_length": 5000,
            "log2_pseudocount": 1,
            "log2_scale_method": "None",
            "heatmap_color_map": "RdBu_r",
            "heatmap_z_min": -2,
            "heatmap_z_max": 2,
            "plot_dpi": 200,
        },
        "chip_controls": chip_controls,
        "chip_samples": chip_samples,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Convert a ChIP-seq sample manifest TSV into Snakemake config.yaml"
    )
    parser.add_argument("manifest", help="Sample manifest TSV")
    parser.add_argument("--species", choices=["human", "mouse", "rat"], required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--bt2-index", required=True)
    parser.add_argument("--black-list", required=True)
    parser.add_argument(
        "--effective-genome-size",
        required=True,
        type=int,
        help="Positive deepTools effective genome size for the assembly, read length, and mapping policy",
    )
    parser.add_argument(
        "--fastq-screen-conf",
        default="config/fastq_screen.conf",
        help="FastQ Screen database configuration",
    )
    parser.add_argument(
        "--numerator-condition",
        default="treated",
        help="Condition whose signal is the DiffBind numerator (default: treated)",
    )
    parser.add_argument(
        "--reference-condition",
        default="control",
        help="Condition used as the DiffBind denominator/reference (default: control)",
    )
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
        args.effective_genome_size,
        args.fastq_screen_conf,
        args.numerator_condition,
        args.reference_condition,
    )

    with open(args.output, "w") as fh:
        yaml.safe_dump(config, fh, sort_keys=False)

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
