#!/usr/bin/env python3
"""Genome-browser snapshots with pyGenomeTracks.

One panel per ChIP target, replicate tracks overlaid, conditions in contrasting
colours, with the consensus peaks and gene models underneath.

The y-axis is shared within a target and free across targets, which is the only
honest choice: H3K27me3 and CTCF differ in dynamic range by an order of magnitude,
so a common scale would flatten one of them to a line. Within a target the scale
must be shared, or a visual "increase" could just be autoscaling.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

from marks import MarkRegistry

# Resting = blue, stimulated = vermilion. Same encoding as the ggplot figures.
COND_COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]


def build_ini(target: str, rows: pd.DataFrame, conditions: list[str], peaks: Path,
              gtf: Path, ini: Path, ymax: float | None) -> None:
    colour = {c: COND_COLOURS[i % len(COND_COLOURS)] for i, c in enumerate(conditions)}
    lines: list[str] = []

    for _, r in rows.iterrows():
        bw = Path("results/bigwig") / f"{r.sample_id}.log2ratio.bw"
        if not bw.exists():
            continue
        lines += [
            f"[{r.sample_id}]",
            f"file = {bw}",
            f"title = {r.target} {r.condition} r{r.replicate}",
            "height = 1.6",
            f"color = {colour[r.condition]}",
            "min_value = 0",
            *( [f"max_value = {ymax:.2f}"] if ymax else [] ),
            "number_of_bins = 700",
            "nans_to_zeros = true",
            "summary_method = mean",
            "show_data_range = true",
            "file_type = bigwig",
            "",
        ]

    lines += [
        "[spacer]",
        "height = 0.1",
        "",
        f"[{target} peaks]",
        f"file = {peaks}",
        f"title = {target} peaks",
        "height = 0.4",
        "color = #444444",
        "border_color = none",
        "display = collapsed",
        "labels = false",
        "file_type = bed",
        "",
        "[genes]",
        f"file = {gtf}",
        "title = genes",
        "height = 2",
        "color = black",
        "prefered_name = gene_name",
        "merge_transcripts = true",
        "labels = true",
        "style = UCSC",
        "fontsize = 8",
        "file_type = gtf",
        "",
        "[x-axis]",
        "fontsize = 8",
        "",
    ]
    ini.write_text("\n".join(lines))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True)
    p.add_argument("--registry", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--loci", nargs="+", required=True, help="NAME:chrom:start-end")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    reg = MarkRegistry(args.registry)

    samples = pd.read_csv(args.samples, sep="\t", dtype=str).fillna("")
    chip = samples[samples.assay == "chip"]
    targets = sorted(chip.target.unique())
    conditions = list(dict.fromkeys(chip.condition))

    # pyGenomeTracks wants an uncompressed, sorted GTF.
    gtf = Path(args.gtf)
    plain = outdir / "genes.gtf"
    if not plain.exists():
        with open(plain, "w") as fh:
            src = ["zcat", str(gtf)] if gtf.suffix == ".gz" else ["cat", str(gtf)]
            zcat = subprocess.Popen(src, stdout=subprocess.PIPE)
            subprocess.run(
                ["awk", '$3=="gene" || $3=="transcript" || $3=="exon"'],
                stdin=zcat.stdout, stdout=fh, check=True,
            )
            zcat.wait()

    made, failed = 0, 0
    for locus in args.loci:
        name, coords = locus.split(":", 1)
        region = coords  # chrom:start-end

        for target in targets:
            rows = chip[chip.target == target].sort_values(["condition", "replicate"])
            peaks = Path("results/peaks/consensus") / f"{target}.bed"
            if not peaks.exists():
                continue

            ini = outdir / f"{name}_{target}.ini"
            out = outdir / f"{name}_{target}.pdf"
            # A shared y-max within the target; broad marks need headroom.
            ymax = 6.0 if reg.is_broad(target) else 10.0
            build_ini(target, rows, conditions, peaks, plain, ini, ymax)

            res = subprocess.run(
                ["pyGenomeTracks", "--tracks", str(ini), "--region", region,
                 "--outFileName", str(out), "--dpi", "300",
                 "--title", f"{name} — {target} ({reg.peak_mode(target)})",
                 "--width", "18"],
                capture_output=True, text=True,
            )
            if res.returncode != 0:
                failed += 1
                print(f"[warn] {name}/{target}: {res.stderr.strip()[-300:]}", file=sys.stderr)
            else:
                made += 1
                print(f"  {out}")

    print(f"[done] {made} browser panels ({failed} failed)")
    if made == 0:
        sys.exit("No browser panels were produced")


if __name__ == "__main__":
    main()
