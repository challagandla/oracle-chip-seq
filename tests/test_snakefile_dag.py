import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest
import yaml


ROOT = Path(__file__).parents[1]


def _record(sample_id, factor, condition, replicate, control, peak_mode, fastq_dir):
    left = fastq_dir / f"arbitrary_readset_{sample_id}_left.fq.gz"
    right = fastq_dir / f"arbitrary_readset_{sample_id}_right.fq.gz"
    left.touch()
    right.touch()
    return {
        "id": sample_id,
        "fastq": [str(left), str(right)],
        "control": control,
        "condition": condition,
        "replicate": replicate,
        "factor": factor,
        "tissue": "cell_line",
        "peak_mode": peak_mode,
    }


def _synthetic_config(tmp_path):
    fastq_dir = tmp_path / "fastq inputs"
    reference_dir = tmp_path / "reference files"
    fastq_dir.mkdir()
    reference_dir.mkdir()

    controls = []
    for condition in ("control", "treated"):
        for replicate in (1, 2):
            sample_id = f"Input_{condition}_{replicate}"
            left = fastq_dir / f"unrelated_{sample_id}_left.fq.gz"
            right = fastq_dir / f"unrelated_{sample_id}_right.fq.gz"
            left.touch()
            right.touch()
            controls.append(
                {
                    "id": sample_id,
                    "fastq": [str(left), str(right)],
                    "condition": condition,
                }
            )

    samples = []
    for factor, peak_mode in (("CTCF", "narrow"), ("H3K27me3", "broad")):
        for condition in ("control", "treated"):
            for replicate in (1, 2):
                samples.append(
                    _record(
                        f"{factor}_{condition}_{replicate}",
                        factor,
                        condition,
                        replicate,
                        f"Input_{condition}_{replicate}",
                        peak_mode,
                        fastq_dir,
                    )
                )

    genome = reference_dir / "genome.fa"
    chrom_sizes = reference_dir / "genome.chrom.sizes"
    blacklist = reference_dir / "black list.bed"
    screen = reference_dir / "fastq screen.conf"
    for path in (genome, chrom_sizes, blacklist, screen):
        path.touch()
    Path(f"{genome}.fai").touch()
    index = reference_dir / "genome index"
    for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        Path(f"{index}{suffix}").touch()

    return {
        "species": "human",
        "references": {
            "human": {
                "genome": str(genome),
                "chrom_sizes": str(chrom_sizes),
                "black_list": str(blacklist),
                "bt2_index": str(index),
                "name": "test",
                "gsize": "2.7e9",
                "effective_genome_size": 2862010428,
            }
        },
        "contamination": {"fastq_screen_conf": str(screen), "subset": 1000},
        "alignment": {"min_mapq": 30, "max_insert_size": 2000},
        "peak_calling": {
            "narrow_qvalue": 0.01,
            "broad_cutoff": 0.1,
            "consensus_min_replicates": 2,
        },
        "differential_binding": {
            "numerator_condition": "treated",
            "reference_condition": "control",
            "narrow_summits": 200,
        },
        "motif_enrichment": {
            "window_bp": 200,
            "background_multiplier": 2,
            "seed": 1,
        },
        "deeptools": {
            "threads": 2,
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
            "plot_dpi": 100,
        },
        "chip_controls": controls,
        "chip_samples": samples,
    }


def test_mixed_peak_modes_build_scoped_visualization_dag(tmp_path):
    snakemake = shutil.which("snakemake")
    if snakemake is None:
        pytest.skip("snakemake executable is not on PATH")

    config = _synthetic_config(tmp_path)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    completed = subprocess.run(
        [
            snakemake,
            "--snakefile",
            str(ROOT / "Snakefile"),
            "--directory",
            str(ROOT),
            "--configfile",
            str(config_path),
            "--cores",
            "2",
            "--dry-run",
            "--forceall",
            "--printshellcmds",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode == 0, output
    assert re.search(r"^total\s+135$", output, flags=re.MULTILINE), output

    assert re.search(r"call_peaks_narrow\s+4", output)
    assert re.search(r"call_peaks_broad\s+4", output)
    assert re.search(r"condition_consensus_peaks\s+4", output)
    assert re.search(r"motif_enrichment\s+1", output)
    assert "results/motifs/CTCF/motif_enrichment.tsv" in output
    assert "results/motifs/CTCF/motif_summary.csv" in output
    assert "results/motifs/CTCF/motif_summary.pdf" in output
    assert "results/motifs/H3K27me3/motif_summary.pdf" not in output
    assert str(Path(config["references"]["human"]["black_list"])) in output
    assert re.search(r"\s200\s+2\s+1\s+CTCF\s", output)

    assert "macs3 callpeak" in output
    assert "results/bam/raw/CTCF_control_1.bam" in output
    assert "bedtools sort -i - -g" in output
    assert "computeMatrix reference-point --referencePoint center" in output
    assert "computeMatrix scale-regions --regionBodyLength 5000" in output
    assert "--sortRegions descend --sortUsing mean" in output
    assert "--sortRegions keep --missingDataColor" in output
    assert "--colorMap RdBu_r" in output
    assert "-X 2000" in output
    assert "samtools view -b -q 30 -F 1804 -f 2" in output
    assert "samtools markdup -r -s" in output
    assert "results/qc/samtools/CTCF_control_1.markdup.txt" in output
    assert "contrast.tsv" in output
    assert " treated control 200" in output

    matrix_lines = [line for line in output.splitlines() if "computeMatrix " in line]
    assert matrix_lines
    assert all("--regionsLabel" not in line for line in matrix_lines)
    pca_lines = [line for line in output.splitlines() if "plotPCA " in line]
    assert pca_lines
    assert all("--transpose" in line and "--log2" not in line for line in pca_lines)
    assert "arbitrary_readset_CTCF_control_1_left.fq.gz" in output
    assert 'CTCF_control_1_R1.fastq.gz"' in output


def test_broad_factor_motif_target_is_rejected(tmp_path):
    snakemake = shutil.which("snakemake")
    if snakemake is None:
        pytest.skip("snakemake executable is not on PATH")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(_synthetic_config(tmp_path), sort_keys=False))
    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    completed = subprocess.run(
        [
            snakemake,
            "--snakefile",
            str(ROOT / "Snakefile"),
            "--directory",
            str(ROOT),
            "--configfile",
            str(config_path),
            "--cores",
            "2",
            "--dry-run",
            "results/motifs/H3K27me3/motif_summary.pdf",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "MissingRuleException" in output, output


def test_synthetic_workflow_passes_snakemake_lint(tmp_path):
    snakemake = shutil.which("snakemake")
    if snakemake is None:
        pytest.skip("snakemake executable is not on PATH")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(_synthetic_config(tmp_path), sort_keys=False))
    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    completed = subprocess.run(
        [
            snakemake,
            "--snakefile",
            str(ROOT / "Snakefile"),
            "--directory",
            str(ROOT),
            "--configfile",
            str(config_path),
            "--lint",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode == 0, output


def test_environment_smoke_workflow_renders_all_commands(tmp_path):
    snakemake = shutil.which("snakemake")
    if snakemake is None:
        pytest.skip("snakemake executable is not on PATH")

    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    completed = subprocess.run(
        [
            snakemake,
            "--snakefile",
            str(ROOT / "workflow" / "envs.smk"),
            "--directory",
            str(ROOT),
            "--cores",
            "1",
            "--dry-run",
            "--forceall",
            "--printshellcmds",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode == 0, output
    assert "mkdir -p .snakemake/setup-env-checks" in output
    assert "library(DiffBind)" in output
