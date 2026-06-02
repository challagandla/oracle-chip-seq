#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

SPECIES_CONFIG = {
    "human": {
        "alias": "hg38",
        "fasta_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz",
        "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz",
    },
    "mouse": {
        "alias": "mm10",
        "fasta_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.primary_assembly.genome.fa.gz",
        "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.annotation.gtf.gz",
    },
    "rat": {
        "alias": "rn6",
        "fasta_url": "https://ftp.ensembl.org/pub/release-113/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.primary_assembly.fa.gz",
        "gtf_url": "https://ftp.ensembl.org/pub/release-113/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.113.gtf.gz",
    },
}


def run(cmd):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        print(f"Skipping existing file: {dest}")
        return dest
    run(["wget", "-c", url, "-O", str(dest)])
    return dest


def decompress(src):
    if src.suffix == ".gz":
        out = src.with_suffix("")
        if out.exists():
            print(f"Skipping existing decompressed file: {out}")
            return out
        run(["gunzip", "-kf", str(src)])
        return out
    return src


def build_bowtie2_index(fasta, prefix):
    print(f"Building Bowtie2 index: {prefix}")
    run(["bowtie2-build", str(fasta), str(prefix)])


def build_salmon_index(fasta, index_dir):
    print(f"Building Salmon index: {index_dir}")
    run(["salmon", "index", "-t", str(fasta), "-i", str(index_dir), "--type", "quasi"] )


def make_chrom_sizes(fasta, out_path):
    print(f"Generating chromosome sizes: {out_path}")
    run(["samtools", "faidx", str(fasta)])
    with open(out_path, "w") as out:
        with open(str(fasta) + ".fai") as fai:
            for line in fai:
                parts = line.strip().split("\t")
                out.write(f"{parts[0]}\t{parts[1]}\n")


def main():
    parser = argparse.ArgumentParser(description="Download reference data for human/mouse/rat and build indexes.")
    parser.add_argument("species", choices=SPECIES_CONFIG.keys(), help="Species name")
    parser.add_argument("--outdir", default="references", help="Directory to store downloaded references")
    parser.add_argument("--build-salmon", action="store_true", help="Build a Salmon index after download")
    args = parser.parse_args()

    cfg = SPECIES_CONFIG[args.species]
    outdir = Path(args.outdir) / args.species
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_gz = download(cfg["fasta_url"], outdir / f"{cfg['alias']}.fa.gz")
    gtf_gz = download(cfg["gtf_url"], outdir / f"{cfg['alias']}.gtf.gz")

    fasta = decompress(fasta_gz)
    gtf = decompress(gtf_gz)

    build_bowtie2_index(fasta, outdir / cfg["alias"])
    make_chrom_sizes(fasta, outdir / f"{cfg['alias']}.chrom.sizes")

    if args.build_salmon:
        build_salmon_index(fasta, outdir / "salmon_index")

    chrom_sizes = outdir / f"{cfg['alias']}.chrom.sizes"
    print("Reference download and preparation complete.")
    print(f"Genome: {fasta}")
    print(f"Annotation: {gtf}")
    print(f"Bowtie2 index prefix: {outdir / cfg['alias']}")
    print(f"Chrom sizes: {chrom_sizes}")


if __name__ == "__main__":
    main()
