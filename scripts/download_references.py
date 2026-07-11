#!/usr/bin/env python3
"""Fetch the reference files config.yaml expects.

Downloads the genome FASTA, chrom sizes, ENCODE blacklist, GENCODE annotation and
a prebuilt Bowtie2 index. The index is pulled prebuilt rather than built locally:
bowtie2-build on a mammalian genome takes hours, and the published no-alt analysis
set is the one ENCODE aligns against, so rolling our own would be both slower and
less standard.

    python3 scripts/download_references.py --species human --outdir references/hg38
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SPECIES = {
    "human": {
        "name": "hg38",
        "genome": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "chrom_sizes": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
        "blacklist": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
        "gtf": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz",
        "bt2_index": "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip",
        "bt2_prefix": "GRCh38_noalt_as",
    },
    "mouse": {
        "name": "mm39",
        "genome": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz",
        "chrom_sizes": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes",
        "blacklist": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
        "gtf": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.annotation.gtf.gz",
        "bt2_index": "https://genome-idx.s3.amazonaws.com/bt/mm39.zip",
        "bt2_prefix": "mm39",
    },
}


def fetch(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name}")
        return
    print(f"[get ] {dest.name}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    subprocess.run(["curl", "-fsSL", "-o", str(tmp), url], check=True)
    tmp.rename(dest)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--species", choices=sorted(SPECIES), required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--skip-index", action="store_true",
                   help="Skip the Bowtie2 index (~4 GB) if one is already available")
    args = p.parse_args()

    spec = SPECIES[args.species]
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    name = spec["name"]

    genome_gz = out / f"{name}.fa.gz"
    genome = out / f"{name}.fa"
    fetch(spec["genome"], genome_gz)
    if not genome.exists():
        print(f"[gunzip] {genome.name}")
        with open(genome, "wb") as fh:
            subprocess.run(["gunzip", "-c", str(genome_gz)], stdout=fh, check=True)
    # The motif step reads sequence straight from the FASTA and needs the index.
    if not (out / f"{name}.fa.fai").exists():
        subprocess.run(["samtools", "faidx", str(genome)], check=True)

    fetch(spec["chrom_sizes"], out / f"{name}.chrom.sizes")
    fetch(spec["blacklist"], out / f"{name}-blacklist.v2.bed.gz")
    fetch(spec["gtf"], out / "annotation.gtf.gz")

    if not args.skip_index:
        idx_dir = out / "bowtie2"
        idx_dir.mkdir(exist_ok=True)
        prefix = idx_dir / spec["bt2_prefix"]
        if not Path(f"{prefix}.1.bt2").exists():
            zipf = idx_dir / "index.zip"
            fetch(spec["bt2_index"], zipf)
            subprocess.run(["unzip", "-o", "-j", str(zipf), "-d", str(idx_dir)], check=True)
            zipf.unlink()

    print(f"\nDone. Point config.yaml at these:")
    print(f"  genome:      {genome.resolve()}")
    print(f"  chrom_sizes: {(out / f'{name}.chrom.sizes').resolve()}")
    print(f"  blacklist:   {(out / f'{name}-blacklist.v2.bed.gz').resolve()}")
    print(f"  gtf:         {(out / 'annotation.gtf.gz').resolve()}")
    if not args.skip_index:
        print(f"  bt2_index:   {(out / 'bowtie2' / spec['bt2_prefix']).resolve()}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        sys.exit(f"command failed: {e}")
