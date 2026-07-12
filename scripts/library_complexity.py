#!/usr/bin/env python3
"""ENCODE library-complexity metrics: NRF, PBC1, PBC2.

A ChIP library sequenced from too little starting material is dominated by PCR
duplicates of a few fragments. Deduplication hides this — the deduplicated BAM
looks fine while the underlying library has almost no distinct information. These
three numbers expose it:

    NRF  = distinct fragment positions / total fragments
    PBC1 = positions seen exactly once / distinct positions
    PBC2 = positions seen exactly once / positions seen exactly twice

A "fragment position" is (chrom, 5' start, strand) for single-end, and the
(chrom, start, mate start) tuple for paired-end, since for PE the true fragment
is known and strand is redundant.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from collections import Counter


def positions(bam: str, mapq: int, paired: bool):
    if paired:
        # -f 0x42: read paired + first in pair. Restricting to read1 counts each
        # fragment once instead of twice.
        flags = ["-f", "0x42", "-F", "0x904"]
    else:
        flags = ["-F", "0x904"]

    cmd = ["samtools", "view", "-q", str(mapq), *flags, bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    for line in proc.stdout:
        f = line.split("\t", 9)
        flag = int(f[1])
        chrom, pos = f[2], f[3]
        if paired:
            yield (chrom, pos, f[7])  # RNEXT position of the mate
        else:
            strand = "-" if flag & 16 else "+"
            yield (chrom, pos, strand)
    proc.stdout.close()
    if proc.wait() != 0:
        sys.exit(f"samtools view failed on {bam}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--mapq", type=int, default=30)
    p.add_argument("--paired", default="0")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    paired = args.paired in ("1", "true", "True", "paired")
    counts = Counter(positions(args.bam, args.mapq, paired))

    total = sum(counts.values())
    distinct = len(counts)
    once = sum(1 for c in counts.values() if c == 1)
    twice = sum(1 for c in counts.values() if c == 2)

    nrf = distinct / total if total else 0.0
    pbc1 = once / distinct if distinct else 0.0
    pbc2 = (once / twice) if twice else float("inf")

    with open(args.out, "w") as fh:
        fh.write("sample\ttotal_fragments\tdistinct_positions\tNRF\tPBC1\tPBC2\n")
        pbc2_s = "inf" if pbc2 == float("inf") else f"{pbc2:.4f}"
        fh.write(
            f"{args.sample}\t{total}\t{distinct}\t{nrf:.4f}\t{pbc1:.4f}\t{pbc2_s}\n"
        )
    print(f"{args.sample}: NRF={nrf:.3f} PBC1={pbc1:.3f} PBC2={pbc2_s} (n={total})")


if __name__ == "__main__":
    main()
