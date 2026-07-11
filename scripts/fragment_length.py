#!/usr/bin/env python3
"""Estimate fragment length, and for single-end libraries the ChIP enrichment
quality metrics that depend on it.

Paired-end: the fragment is observed, so take the median template length.

Single-end: the fragment is not observed. Reads pile up on the + strand upstream
of a binding site and on the - strand downstream, so the strand cross-correlation
peaks at the fragment length. MACS2's `predictd` implements this, and the value it
returns is what we hand back as --extsize. Getting it wrong smears peaks (too
large) or splits each peak in two (too small), so this is a load-bearing number,
not a QC nicety.

Falls back to the configured default only if prediction fails outright — which
itself is a red flag, because a library with no cross-correlation signal has no
ChIP enrichment.
"""
from __future__ import annotations

import argparse
import re
import statistics
import subprocess
import sys

FALLBACK_FRAGMENT = 200


def paired_fragment(bam: str, limit: int = 2_000_000) -> tuple[int, str]:
    cmd = ["samtools", "view", "-f", "0x42", "-F", "0x904", bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    lengths = []
    for line in proc.stdout:
        tlen = abs(int(line.split("\t", 9)[8]))
        if 0 < tlen < 2000:
            lengths.append(tlen)
        if len(lengths) >= limit:
            break
    proc.stdout.close()
    proc.wait()
    if not lengths:
        return FALLBACK_FRAGMENT, "fallback (no proper pairs with usable TLEN)"
    return int(statistics.median(lengths)), "observed insert size (paired-end)"


def single_fragment(bam: str, gsize: str) -> tuple[int, str]:
    cmd = ["macs2", "predictd", "-i", bam, "-g", gsize, "--outdir", "/tmp"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    text = res.stderr + res.stdout
    m = re.search(r"predicted fragment length is (\d+)", text)
    if m:
        return int(m.group(1)), "strand cross-correlation (macs2 predictd)"
    m = re.search(r"alternative fragment length\(s\) may be (\d+)", text)
    if m:
        return int(m.group(1)), "cross-correlation, alternative peak"
    print(text[-1500:], file=sys.stderr)
    return FALLBACK_FRAGMENT, "FALLBACK — cross-correlation failed; suspect poor enrichment"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--paired", default="0")
    p.add_argument("--gsize", default="2.7e9")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    paired = args.paired in ("1", "true", "True", "paired")
    if paired:
        frag, method = paired_fragment(args.bam)
    else:
        frag, method = single_fragment(args.bam, args.gsize)

    with open(args.out, "w") as fh:
        fh.write("sample\tlayout\tfragment_length\tmethod\n")
        fh.write(f"{args.sample}\t{'paired' if paired else 'single'}\t{frag}\t{method}\n")
    print(f"{args.sample}: fragment={frag} bp via {method}")


if __name__ == "__main__":
    main()
