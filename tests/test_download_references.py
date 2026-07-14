import gzip
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))

import download_references as references  # noqa: E402


def _gzip_fasta(sequence=b">chr1\nACGT\n"):
    return gzip.compress(sequence)


def test_download_reuses_only_a_valid_complete_file(tmp_path, monkeypatch):
    destination = tmp_path / "genome.fa.gz"
    destination.write_bytes(_gzip_fasta())

    def unexpected_run(_cmd):
        raise AssertionError("a validated download should not be fetched again")

    monkeypatch.setattr(references, "run", unexpected_run)
    assert references.download("https://example.invalid/genome.fa.gz", destination) == destination
    assert not (tmp_path / "genome.fa.gz.part").exists()


def test_download_replaces_invalid_existing_file_atomically(tmp_path, monkeypatch):
    destination = tmp_path / "genome.fa.gz"
    destination.write_bytes(b"truncated")
    expected = _gzip_fasta()

    def fake_wget(cmd):
        Path(cmd[-1]).write_bytes(expected)

    monkeypatch.setattr(references, "run", fake_wget)
    references.download("https://example.invalid/genome.fa.gz", destination)

    assert destination.read_bytes() == expected
    assert not (tmp_path / "genome.fa.gz.part").exists()


def test_failed_download_validation_removes_corrupt_partial(tmp_path, monkeypatch):
    destination = tmp_path / "genome.fa.gz"

    def fake_wget(cmd):
        Path(cmd[-1]).write_bytes(b"not gzip data")

    monkeypatch.setattr(references, "run", fake_wget)
    with pytest.raises(subprocess.CalledProcessError):
        references.download("https://example.invalid/genome.fa.gz", destination)

    assert not destination.exists()
    assert not (tmp_path / "genome.fa.gz.part").exists()


@pytest.mark.skipif(shutil.which("gunzip") is None, reason="gunzip is not installed")
def test_decompress_replaces_invalid_fasta_via_atomic_partial(tmp_path):
    source = tmp_path / "genome.fa.gz"
    output = tmp_path / "genome.fa"
    source.write_bytes(_gzip_fasta())
    output.write_text("truncated sequence without a header\n")

    assert references.decompress(source) == output
    assert output.read_text() == ">chr1\nACGT\n"
    assert not (tmp_path / "genome.fa.part").exists()
    assert (tmp_path / "genome.fa.complete").is_file()


@pytest.mark.skipif(shutil.which("gunzip") is None, reason="gunzip is not installed")
def test_decompress_rebuilds_header_only_legacy_fasta_without_marker(tmp_path):
    source = tmp_path / "genome.fa.gz"
    output = tmp_path / "genome.fa"
    source.write_bytes(_gzip_fasta())
    output.write_text(">chr1\n")

    references.decompress(source)

    assert output.read_text() == ">chr1\nACGT\n"
    references.validate_decompression_marker(source, output)


@pytest.mark.skipif(shutil.which("gunzip") is None, reason="gunzip is not installed")
def test_decompress_rebuilds_when_completed_output_is_later_truncated(tmp_path):
    source = tmp_path / "genome.fa.gz"
    source.write_bytes(_gzip_fasta())
    output = references.decompress(source)
    output.write_text(">chr1\n")

    references.decompress(source)

    assert output.read_text() == ">chr1\nACGT\n"
    references.validate_decompression_marker(source, output)
