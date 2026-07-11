#!/usr/bin/env python3
"""Fetch FASTQ files for the accessions listed in a sample table.

Runs prefetch + fasterq-dump, then compresses. Output names follow the
convention the workflow expects:

    single:  data/raw/{sample_id}.fastq.gz
    paired:  data/raw/{sample_id}_R1.fastq.gz, data/raw/{sample_id}_R2.fastq.gz

The declared layout in the sample table is checked against what fasterq-dump
actually produces; a mismatch is an error rather than a silent single-end run.
"""
import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def compress(src: Path, dest: Path, threads: int) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    pigz = shutil.which("pigz")
    with open(tmp, "wb") as out:
        if pigz:
            run([pigz, "-c", "-p", str(threads), str(src)], stdout=out)
        else:
            run(["gzip", "-c", str(src)], stdout=out)
    tmp.rename(dest)
    src.unlink()


def expected_outputs(sample_id: str, layout: str, outdir: Path):
    if layout == "paired":
        return [outdir / f"{sample_id}_R1.fastq.gz", outdir / f"{sample_id}_R2.fastq.gz"]
    return [outdir / f"{sample_id}.fastq.gz"]


def fetch(row, outdir: Path, threads: int, max_spots: int | None) -> None:
    sample_id, srr, layout = row["sample_id"], row["srr"], row["layout"].lower()
    targets = expected_outputs(sample_id, layout, outdir)
    if all(t.exists() and t.stat().st_size > 0 for t in targets):
        print(f"[skip] {sample_id} ({srr}) already present", flush=True)
        return

    print(f"[fetch] {sample_id} <- {srr} ({layout})", flush=True)
    with tempfile.TemporaryDirectory(prefix=f"sra_{srr}_", dir=outdir) as tmp:
        tmpdir = Path(tmp)
        run(["prefetch", "--max-size", "100G", "--output-directory", str(tmpdir), srr])

        cmd = ["fasterq-dump", "--split-3", "--threads", str(threads), "--outdir", str(tmpdir)]
        if max_spots:
            # fasterq-dump has no spot cap, so bound the read via fastq-dump instead.
            cmd = ["fastq-dump", "--split-3", "--maxSpotId", str(max_spots), "--outdir", str(tmpdir)]
        cmd.append(str(tmpdir / srr))
        run(cmd)

        r1, r2 = tmpdir / f"{srr}_1.fastq", tmpdir / f"{srr}_2.fastq"
        single = tmpdir / f"{srr}.fastq"

        if r1.exists() and r2.exists():
            if layout != "paired":
                raise SystemExit(
                    f"{sample_id} ({srr}) is declared '{layout}' but SRA returned paired reads"
                )
            compress(r1, targets[0], threads)
            compress(r2, targets[1], threads)
        elif single.exists():
            if layout != "single":
                raise SystemExit(
                    f"{sample_id} ({srr}) is declared '{layout}' but SRA returned single reads"
                )
            compress(single, targets[0], threads)
        else:
            raise SystemExit(f"{sample_id} ({srr}): fasterq-dump produced no FASTQ")

    for t in targets:
        print(f"  -> {t} ({t.stat().st_size / 1e9:.2f} GB)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", default="samples.tsv")
    p.add_argument("--outdir", default="data/raw")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument(
        "--max-spots",
        type=int,
        default=None,
        help="Cap reads per run. Use for a fast smoke test; omit for the real analysis.",
    )
    p.add_argument("--only", nargs="*", help="Restrict to these sample_ids")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.samples, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("sample_id")]
    if args.only:
        rows = [r for r in rows if r["sample_id"] in set(args.only)]
    if not rows:
        sys.exit("No samples selected")

    for row in rows:
        fetch(row, outdir, args.threads, args.max_spots)
    print(f"[done] {len(rows)} samples in {outdir}", flush=True)


if __name__ == "__main__":
    main()
