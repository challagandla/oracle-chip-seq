from copy import deepcopy
from pathlib import Path
import sys

import pytest


sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))

from build_sample_sheets import build_diffbind_sheet, validate_design  # noqa: E402
from manifest_to_config import build_config, parse_manifest  # noqa: E402


def sample(sample_id, factor, condition, replicate, control, peak_mode="narrow"):
    return {
        "id": sample_id,
        "fastq": [f"{sample_id}_R1.fastq.gz", f"{sample_id}_R2.fastq.gz"],
        "control": control,
        "condition": condition,
        "replicate": replicate,
        "factor": factor,
        "tissue": "test",
        "peak_mode": peak_mode,
    }


@pytest.fixture
def valid_config():
    controls = [
        {
            "id": f"input_{condition}_{rep}",
            "fastq": [
                f"input_{condition}_{rep}_R1.fastq.gz",
                f"input_{condition}_{rep}_R2.fastq.gz",
            ],
            "condition": condition,
        }
        for condition in ("control", "treated")
        for rep in (1, 2)
    ]
    samples = [
        sample(
            f"H3K27ac_{condition}_{rep}",
            "H3K27ac",
            condition,
            rep,
            f"input_{condition}_{rep}",
        )
        for condition in ("control", "treated")
        for rep in (1, 2)
    ]
    return {"chip_controls": controls, "chip_samples": samples}


def test_valid_design_builds_factor_specific_sheet(valid_config):
    sheet = build_diffbind_sheet(valid_config)
    assert set(sheet["Factor"]) == {"H3K27ac"}
    assert set(sheet["Peaks"].str.rsplit(".", n=1).str[-1]) == {"narrowPeak"}
    assert set(sheet["PeakMode"]) == {"narrow"}


def test_factor_condition_confounding_is_rejected(valid_config):
    cfg = deepcopy(valid_config)
    cfg["chip_samples"] = [
        sample("H3K27ac_control_1", "H3K27ac", "control", 1, "input_control_1"),
        sample("H3K27ac_control_2", "H3K27ac", "control", 2, "input_control_2"),
        sample("CTCF_treated_1", "CTCF", "treated", 1, "input_treated_1"),
        sample("CTCF_treated_2", "CTCF", "treated", 2, "input_treated_2"),
    ]
    with pytest.raises(ValueError, match="requires exactly two conditions"):
        validate_design(cfg)


def test_peak_mode_is_explicit(valid_config):
    cfg = deepcopy(valid_config)
    del cfg["chip_samples"][0]["peak_mode"]
    with pytest.raises(ValueError, match="peak_mode"):
        validate_design(cfg)


def test_duplicate_replicate_numbers_are_rejected(valid_config):
    cfg = deepcopy(valid_config)
    cfg["chip_samples"][1]["replicate"] = 1
    with pytest.raises(ValueError, match="duplicate replicate numbers"):
        validate_design(cfg)


def test_missing_factor_is_rejected_instead_of_inferred_from_id(valid_config):
    cfg = deepcopy(valid_config)
    del cfg["chip_samples"][0]["factor"]
    with pytest.raises(ValueError, match="missing factor"):
        validate_design(cfg)


def test_sample_id_whitespace_is_rejected(valid_config):
    cfg = deepcopy(valid_config)
    cfg["chip_samples"][0]["id"] = " H3K27ac_control_1"
    with pytest.raises(ValueError, match="leading or trailing whitespace"):
        validate_design(cfg)


def test_uncompressed_fastq_extension_is_rejected(valid_config):
    cfg = deepcopy(valid_config)
    cfg["chip_samples"][0]["fastq"][0] = "reads_R1.fastq"
    with pytest.raises(ValueError, match="must end in .fastq.gz or .fq.gz"):
        validate_design(cfg)


def test_same_fastq_cannot_be_used_for_both_mates(valid_config):
    cfg = deepcopy(valid_config)
    cfg["chip_samples"][0]["fastq"][1] = cfg["chip_samples"][0]["fastq"][0]
    with pytest.raises(ValueError, match="same file for R1 and R2"):
        validate_design(cfg)


def test_fastq_cannot_be_assigned_to_multiple_samples(valid_config):
    cfg = deepcopy(valid_config)
    cfg["chip_samples"][1]["fastq"][0] = cfg["chip_samples"][0]["fastq"][0]
    with pytest.raises(ValueError, match="assigned more than once"):
        validate_design(cfg)


def test_tissue_is_required(valid_config):
    cfg = deepcopy(valid_config)
    del cfg["chip_samples"][0]["tissue"]
    with pytest.raises(ValueError, match="missing tissue"):
        validate_design(cfg)


def test_control_condition_must_match_chip_condition(valid_config):
    cfg = deepcopy(valid_config)
    cfg["chip_samples"][0]["control"] = "input_treated_1"
    with pytest.raises(ValueError, match="does not match control"):
        validate_design(cfg)


def test_condition_whitespace_is_normalized_in_sheet(valid_config):
    cfg = deepcopy(valid_config)
    cfg["chip_samples"][0]["condition"] = " control "
    sheet = build_diffbind_sheet(cfg)
    assert set(sheet["Condition"]) == {"control", "treated"}


def test_multiple_tissues_require_an_explicit_extended_model(valid_config):
    cfg = deepcopy(valid_config)
    for index, record in enumerate(cfg["chip_samples"]):
        record["tissue"] = "brain" if index < 2 else "liver"
    with pytest.raises(ValueError, match="requires one tissue"):
        validate_design(cfg)


def test_manifest_config_has_integer_effective_size_and_screen_config(valid_config):
    rows = []
    for control in valid_config["chip_controls"]:
        condition, replicate = control["id"].rsplit("_", 2)[1:]
        rows.append(
            {
                "sample_id": control["id"],
                "type": "control",
                "condition": condition,
                "replicate": replicate,
                "fastq_r1": control["fastq"][0],
                "fastq_r2": control["fastq"][1],
                "control": "",
                "factor": "",
                "tissue": "cell_line",
                "peak_mode": "",
            }
        )
    for chip in valid_config["chip_samples"]:
        rows.append(
            {
                "sample_id": chip["id"],
                "type": "chip",
                "condition": chip["condition"],
                "replicate": str(chip["replicate"]),
                "fastq_r1": chip["fastq"][0],
                "fastq_r2": chip["fastq"][1],
                "control": chip["control"],
                "factor": chip["factor"],
                "tissue": chip["tissue"],
                "peak_mode": chip["peak_mode"],
            }
        )

    cfg = build_config(
        rows,
        "human",
        "hg38.fa",
        "hg38.chrom.sizes",
        "hg38",
        "hg38-blacklist.bed",
        2862010428,
        "fastq_screen.conf",
    )
    ref = cfg["references"]["human"]
    assert isinstance(ref["effective_genome_size"], int)
    assert ref["effective_genome_size"] == 2862010428
    assert cfg["contamination"]["fastq_screen_conf"] == "fastq_screen.conf"
    assert cfg["differential_binding"] == {
        "numerator_condition": "treated",
        "reference_condition": "control",
        "narrow_summits": 200,
    }


def test_shipped_manifest_builds_the_documented_k562_design():
    rows = parse_manifest(Path(__file__).parents[1] / "sample_manifest.tsv")
    cfg = build_config(
        rows,
        "human",
        "hg38.fa",
        "hg38.chrom.sizes",
        "hg38",
        "hg38-blacklist.bed",
        2862010428,
        "fastq_screen.conf",
    )

    assert {sample["tissue"] for sample in cfg["chip_samples"]} == {"K562"}
    assert {control["condition"] for control in cfg["chip_controls"]} == {
        "control",
        "treated",
    }
