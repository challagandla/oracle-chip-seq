#!/usr/bin/env python3
import argparse
import os
import subprocess
from pathlib import Path

SPECIES_CONFIG = {
    "human": {
        "alias": "hg38",
        "fasta_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz",
    },
    "mouse": {
        "alias": "mm39",
        "fasta_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.primary_assembly.genome.fa.gz",
    },
    "rat": {
        "alias": "rn6",
        "fasta_url": "https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz",
    },
}


def run(cmd):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def validate_download(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Downloaded file is empty or missing: {path}")
    if path.name.endswith(".gz") or path.name.endswith(".gz.part"):
        subprocess.run(["gzip", "-t", str(path)], check=True)


def validate_fasta(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"FASTA is empty or missing: {path}")
    with path.open("rb") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                if not stripped.startswith(b">"):
                    raise ValueError(f"FASTA does not begin with a header line: {path}")
                return
    raise ValueError(f"FASTA contains no records: {path}")


def completion_marker(path):
    return path.with_name(f"{path.name}.complete")


def validate_decompression_marker(src, out):
    marker = completion_marker(out)
    try:
        source_size, output_size = [int(value) for value in marker.read_text().split()]
    except (OSError, ValueError) as error:
        raise ValueError(f"Missing or invalid completion marker: {marker}") from error
    if source_size != src.stat().st_size or output_size != out.stat().st_size:
        raise ValueError(f"Completion marker does not match reference files: {marker}")


def write_decompression_marker(src, out):
    marker = completion_marker(out)
    partial = marker.with_name(f"{marker.name}.part")
    partial.write_text(f"{src.stat().st_size}\t{out.stat().st_size}\n")
    os.replace(partial, marker)


def download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        try:
            validate_download(dest)
        except (OSError, ValueError, subprocess.CalledProcessError):
            print(f"Existing download failed validation; replacing: {dest}")
            dest.unlink()
        else:
            print(f"Using validated existing file: {dest}")
            return dest

    partial = dest.with_name(f"{dest.name}.part")
    run(["wget", "-c", url, "-O", str(partial)])
    try:
        validate_download(partial)
    except (OSError, ValueError, subprocess.CalledProcessError):
        partial.unlink(missing_ok=True)
        raise
    os.replace(partial, dest)
    return dest


def decompress(src):
    if src.suffix == ".gz":
        out = src.with_suffix("")
        marker = completion_marker(out)
        if out.exists():
            try:
                validate_fasta(out)
                validate_decompression_marker(src, out)
            except (OSError, ValueError):
                print(f"Existing FASTA failed validation; replacing: {out}")
                out.unlink()
                marker.unlink(missing_ok=True)
            else:
                print(f"Using validated existing FASTA: {out}")
                return out

        partial = out.with_name(f"{out.name}.part")
        print("Running:", " ".join(["gunzip", "-c", str(src), ">", str(partial)]))
        try:
            with partial.open("wb") as handle:
                subprocess.run(["gunzip", "-c", str(src)], stdout=handle, check=True)
            validate_fasta(partial)
        except (OSError, ValueError, subprocess.CalledProcessError):
            partial.unlink(missing_ok=True)
            raise
        os.replace(partial, out)
        write_decompression_marker(src, out)
        return out
    return src


def build_bowtie2_index(fasta, prefix):
    print(f"Building Bowtie2 index: {prefix}")
    run(["bowtie2-build", str(fasta), str(prefix)])


def make_chrom_sizes(fasta, out_path):
    print(f"Generating chromosome sizes: {out_path}")
    run(["samtools", "faidx", str(fasta)])
    with open(out_path, "w") as out:
        with open(str(fasta) + ".fai") as fai:
            for line in fai:
                parts = line.strip().split("\t")
                out.write(f"{parts[0]}\t{parts[1]}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Download genome FASTA files and build ChIP-seq alignment references."
    )
    parser.add_argument("species", choices=SPECIES_CONFIG.keys(), help="Species name")
    parser.add_argument(
        "--outdir", default="references", help="Directory to store downloaded references"
    )
    args = parser.parse_args()

    cfg = SPECIES_CONFIG[args.species]
    outdir = Path(args.outdir) / args.species
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_gz = download(cfg["fasta_url"], outdir / f"{cfg['alias']}.fa.gz")
    fasta = decompress(fasta_gz)

    build_bowtie2_index(fasta, outdir / cfg["alias"])
    chrom_sizes = outdir / f"{cfg['alias']}.chrom.sizes"
    make_chrom_sizes(fasta, chrom_sizes)

    print("Reference download and preparation complete.")
    print(f"Genome: {fasta}")
    print(f"Bowtie2 index prefix: {outdir / cfg['alias']}")
    print(f"Chrom sizes: {chrom_sizes}")


if __name__ == "__main__":
    main()
