"""Tests for the peak-mode registry.

The registry is the pipeline's central claim: that narrow and broad targets must
be treated differently, everywhere. These tests pin the classifications that are
easy to get wrong and the invariants that keep the two modes from leaking into
each other.

    pytest tests/
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from marks import MarkRegistry  # noqa: E402

REGISTRY = ROOT / "config" / "mark_registry.yaml"


@pytest.fixture(scope="module")
def reg():
    return MarkRegistry(REGISTRY)


# --- the classifications people get wrong ------------------------------------

@pytest.mark.parametrize("mark", ["H3K4me3", "H3K27ac", "H3K9ac", "H3K4me2", "H2AFZ"])
def test_punctate_histones_are_narrow(reg, mark):
    """Acetylation and H3K4me3 are punctate. Calling them broad fuses adjacent
    regulatory elements and destroys the summits downstream steps rely on."""
    assert reg.peak_mode(mark) == "narrow"
    assert reg.peak_ext(mark) == "narrowPeak"


@pytest.mark.parametrize(
    "mark", ["H3K27me3", "H3K9me3", "H3K36me3", "H3K79me2", "H4K20me1", "H3K9me2"]
)
def test_domain_histones_are_broad(reg, mark):
    assert reg.peak_mode(mark) == "broad"
    assert reg.peak_ext(mark) == "broadPeak"


def test_h3k4me1_is_broad(reg):
    """The counter-intuitive one: H3K4me1 marks enhancers but is broad in the
    ENCODE spec — its signal is a wide shoulder around the NDR, not a peak."""
    assert reg.peak_mode("H3K4me1") == "broad"


def test_tfs_are_narrow(reg):
    assert reg.peak_mode("CTCF") == "narrow"
    assert reg.get("CTCF")["class"] == "tf"


# --- consistency between mode and every downstream decision ------------------

def test_broad_never_uses_idr(reg):
    """IDR models rank consistency in a ranked peak list; broad peaks are wide with
    compressed scores and the model does not hold."""
    for mark in reg.marks:
        if reg.is_broad(mark):
            assert reg.uses_idr(mark) is False, f"{mark} is broad but routed to IDR"


def test_broad_never_runs_motifs(reg):
    """Scanning a multi-kb domain for 8-mers recovers its base composition."""
    for mark in reg.marks:
        if reg.is_broad(mark):
            assert not reg.motifs_enabled(mark), f"{mark} is broad but motifs enabled"


def test_broad_uses_background_normalisation(reg):
    """Reads-in-peaks normalisation assumes most peaks are unchanged, which fails
    for domain marks that can shift globally."""
    for mark in reg.marks:
        spec = reg.get(mark)
        if spec["peak_mode"] == "broad":
            assert spec["diffbind"]["normalize"] == "lib", mark
            assert spec["diffbind"]["summits"] is False, (
                f"{mark} is broad; a summit window would discard the domain"
            )
        else:
            assert spec["diffbind"]["normalize"] == "RLE", mark
            assert isinstance(spec["diffbind"]["summits"], int), mark


def test_broad_allows_wide_peaks_narrow_does_not(reg):
    """A 40 kb CTCF peak is an artefact; a 40 kb H3K27me3 peak is the biology."""
    assert reg.get("CTCF")["qc"]["max_peak_width"] <= 5_000
    assert reg.get("H3K27me3")["qc"]["max_peak_width"] >= 100_000


# --- MACS2 command assembly ---------------------------------------------------

def test_broad_flag_only_for_broad(reg):
    assert "--broad" in reg.macs2_args("H3K27me3", "2.7e9", True)
    assert "--broad" not in reg.macs2_args("H3K27ac", "2.7e9", True)


def test_call_summits_stripped_for_broad(reg, tmp_path):
    """--call-summits is incompatible with --broad; MACS2 exits with an error.
    Any registry entry combining them must be sanitised, not passed through."""
    import yaml

    raw = yaml.safe_load(REGISTRY.read_text())
    raw["marks"]["H3K27me3"]["macs2"]["extra"] = "--call-summits"
    p = tmp_path / "r.yaml"
    p.write_text(yaml.safe_dump(raw))

    args = MarkRegistry(p).macs2_args("H3K27me3", "2.7e9", True)
    assert "--broad" in args
    assert "--call-summits" not in args


def test_paired_uses_bampe_single_uses_extsize(reg):
    """Paired-end gives MACS2 the true fragment. Single-end has no fragment, so the
    cross-correlation estimate must be supplied as --extsize or peaks are smeared."""
    pe = reg.macs2_args("CTCF", "2.7e9", is_paired=True)
    assert "--format BAMPE" in pe
    assert "--extsize" not in pe

    se = reg.macs2_args("CTCF", "2.7e9", is_paired=False, extsize=180)
    assert "--format BAM" in se
    assert "--nomodel --extsize 180" in se


# --- fallbacks and aliases ----------------------------------------------------

def test_unknown_target_defaults_to_narrow(reg):
    """An unrecognised antibody is far more likely a TF than a domain mark, and
    over-broad calling is the more damaging error."""
    spec = reg.get("SOME_NOVEL_FACTOR")
    assert spec["known"] is False
    assert spec["peak_mode"] == "narrow"


@pytest.mark.parametrize(
    "alias,canonical",
    [("H3K27AC", "H3K27ac"), ("h3k4me3", "H3K4me3"), ("RNAPII", "POLR2A"),
     ("Pol2", "POLR2A"), ("H2A.Z", "H2AFZ")],
)
def test_aliases_resolve(reg, alias, canonical):
    """GEO/SRA metadata spells targets inconsistently; a spelling variant must not
    silently become an unknown target and fall back to the default."""
    assert reg.get(alias)["target"] == canonical


def test_polr2a_narrow_peaks_but_gene_body_profile(reg):
    """Pol II pauses in a sharp promoter peak but elongates across the gene body,
    so peak mode and profile geometry are genuinely independent fields."""
    spec = reg.get("POLR2A")
    assert spec["peak_mode"] == "narrow"
    assert spec["profile"] == "scale_regions"
