#!/usr/bin/env python3
"""Build the intervals that reads are counted in for differential binding.

This is the summit rule from the registry, applied:

  narrow (summits: N)   Re-centre each consensus peak on its strongest MACS2 summit
                        and take +/- N bp. Called narrow peaks vary in width from
                        ~150 bp to several kb, and a wide peak accumulates counts
                        simply for being wide. Fixed windows on the summit make the
                        count comparable across peaks and concentrate the signal
                        where the factor actually binds, which is what gives narrow
                        marks their sensitivity.

  broad  (summits: false)
                        Count over the whole called interval. A Polycomb domain has
                        no meaningful summit — the biology is the breadth — and
                        clipping it to a 500 bp window around the highest point
                        would discard almost all of the signal.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from marks import MarkRegistry


def load_summits(paths: list[str]) -> dict[str, list[tuple[int, float]]]:
    """chrom -> [(pos, score), ...] from MACS2 *_summits.bed."""
    by_chrom: dict[str, list[tuple[int, float]]] = {}
    for path in paths:
        p = Path(path)
        if not p.exists():
            continue
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            f = line.split("\t")
            try:
                by_chrom.setdefault(f[0], []).append((int(f[1]), float(f[4])))
            except (IndexError, ValueError):
                continue
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def best_summit(summits: list[tuple[int, float]], start: int, end: int) -> int | None:
    import bisect

    lo = bisect.bisect_left(summits, (start, -1.0))
    hi = bisect.bisect_right(summits, (end, float("inf")))
    window = summits[lo:hi]
    if not window:
        return None
    return max(window, key=lambda s: s[1])[0]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--peaks", required=True, help="Consensus peak BED")
    p.add_argument("--target", required=True)
    p.add_argument("--summits", nargs="*", default=[], help="MACS2 *_summits.bed files")
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--registry", default="config/mark_registry.yaml")
    p.add_argument("--out-bed", required=True)
    p.add_argument("--out-saf", required=True)
    args = p.parse_args()

    reg = MarkRegistry(args.registry)
    spec = reg.get(args.target)
    half = spec["diffbind"]["summits"]

    limits = {}
    for line in Path(args.chrom_sizes).read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            limits[f[0]] = int(f[1])

    summit_index = load_summits(args.summits) if half else {}

    rows = []
    n_recentred = 0
    for line in Path(args.peaks).read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        chrom, start, end = f[0], int(f[1]), int(f[2])
        name = f[3] if len(f) > 3 else f"{chrom}:{start}-{end}"

        if half:
            s = best_summit(summit_index.get(chrom, []), start, end)
            if s is None:
                # No summit inside this consensus peak (it came from a replicate
                # whose summit fell just outside the merged interval). Fall back to
                # the peak midpoint rather than dropping the region.
                s = (start + end) // 2
            else:
                n_recentred += 1
            new_start = max(0, s - int(half))
            new_end = min(limits.get(chrom, s + int(half)), s + int(half))
        else:
            new_start, new_end = start, min(limits.get(chrom, end), end)

        if new_end <= new_start:
            continue
        rows.append((chrom, new_start, new_end, name))

    rows.sort(key=lambda r: (r[0], r[1]))

    with open(args.out_bed, "w") as bed, open(args.out_saf, "w") as saf:
        saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for chrom, start, end, name in rows:
            bed.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t.\n")
            # SAF is 1-based inclusive; BED is 0-based half-open.
            saf.write(f"{name}\t{chrom}\t{start + 1}\t{end}\t+\n")

    mode = f"summit +/-{half} bp" if half else "full called interval"
    print(
        f"{args.target} [{spec['peak_mode']}]: {len(rows)} count regions ({mode})"
        + (f", {n_recentred} re-centred on a summit" if half else ""),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
