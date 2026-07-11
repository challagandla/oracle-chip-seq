#!/usr/bin/env python3
"""featureCounts wrapper that handles a mixed single-end / paired-end cohort.

`-p --countReadPairs` is a per-invocation flag, so a cohort containing both
layouts cannot be counted in one call: counting a PE library without `-p` counts
each mate separately and doubles its counts relative to the SE libraries. If the
layouts happen to correlate with the experimental groups — as they do whenever a
study is topped up with a second, paired-end batch — that 2x lands directly on the
contrast and manufactures differential binding out of nothing.

So: count each layout separately, both as fragments, then join on region.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


def run_featurecounts(saf: str, bams: list[str], paired: bool, threads: int, tmp: Path) -> pd.DataFrame:
    out = tmp / ("pe.txt" if paired else "se.txt")
    cmd = [
        "featureCounts",
        "-a", saf, "-F", "SAF",
        "-o", str(out),
        "-T", str(threads),
        "-Q", "30",
        "--fracOverlap", "0.2",
    ]
    if paired:
        # Count fragments, not mates. -B requires both ends mapped; the BAMs are
        # already filtered to proper pairs so this is a no-op safety net.
        cmd += ["-p", "--countReadPairs", "-B"]
    cmd += bams
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stdout + res.stderr)
        sys.exit(f"featureCounts failed ({'PE' if paired else 'SE'})")
    sys.stderr.write(res.stderr[-2000:])

    df = pd.read_csv(out, sep="\t", comment="#")
    df = df.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
    df = df.rename(columns={"Geneid": "region"})
    # featureCounts labels columns with the BAM path; recover the sample id.
    df.columns = ["region"] + [Path(c).name.replace(".filtered.bam", "") for c in df.columns[1:]]
    return df


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--saf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--se", nargs="*", default=[])
    p.add_argument("--pe", nargs="*", default=[])
    args = p.parse_args()

    se = [b for b in args.se if b]
    pe = [b for b in args.pe if b]
    if not se and not pe:
        sys.exit("No BAMs given")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        frames = []
        if se:
            frames.append(run_featurecounts(args.saf, se, False, args.threads, tmp))
        if pe:
            frames.append(run_featurecounts(args.saf, pe, True, args.threads, tmp))

        merged = frames[0]
        for f in frames[1:]:
            merged = merged.merge(f, on="region", how="outer")

    merged = merged.fillna(0)
    for c in merged.columns[1:]:
        merged[c] = merged[c].astype(int)
    merged.to_csv(args.out, sep="\t", index=False)

    print(
        f"{len(merged)} regions x {len(merged.columns)-1} samples "
        f"({len(se)} single-end, {len(pe)} paired-end counted as fragments)"
    )


if __name__ == "__main__":
    main()
