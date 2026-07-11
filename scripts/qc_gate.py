#!/usr/bin/env python3
"""Mark-aware ChIP-seq QC gate.

Aggregates the per-sample metrics and judges each against the threshold for that
particular target, taken from the mark registry. The point of doing it per mark is
that no single cutoff is defensible across targets:

    FRiP 2%   -> a failure for H3K4me3, entirely normal for H3K9me3
    NSC 1.05  -> reasonable for CTCF, unreachable for a broad repressive domain
    45M reads -> the ENCODE floor for broad marks, overkill for a sharp TF

Emits PASS / WARN / FAIL per sample with the reason attached. Nothing is dropped
automatically — a flagged sample stays in the analysis, but the flag travels with
it into the summary so a marginal result is never presented as a clean one.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

from marks import MarkRegistry


def _read_tsvs(paths: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(p, sep="\t") for p in paths if p.exists()]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _flagstat_mapped(path: Path) -> int:
    """Reads surviving the filtering step (the 'usable' count in ENCODE terms)."""
    if not path.exists():
        return 0
    for line in path.read_text().splitlines():
        if " in total " in line:
            return int(line.split()[0])
    return 0


def _alignment_rate(path: Path) -> float | None:
    """Overall alignment rate from the bowtie2 summary.

    Worth checking on its own rather than inferring it from the usable-read count:
    a library can look adequately deep and still be mostly non-human. A rate far
    below the ~90-98% typical of a clean human ChIP means the reads are dominated
    by adapter dimer, another organism, or degraded input — and no amount of
    downstream normalisation repairs that.
    """
    if not path.exists():
        return None
    for line in path.read_text().splitlines():
        if "overall alignment rate" in line:
            try:
                return float(line.split("%")[0].strip()) / 100.0
            except ValueError:
                return None
    return None


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True)
    p.add_argument("--registry", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--qcdir", default="results/qc")
    p.add_argument("--tsv", required=True)
    p.add_argument("--md", required=True)
    args = p.parse_args()

    cfg = yaml.safe_load(open(args.config))
    reg = MarkRegistry(args.registry)
    qc = Path(args.qcdir)

    samples = pd.read_csv(args.samples, sep="\t", dtype=str).fillna("")
    frip = _read_tsvs(sorted((qc / "frip").glob("*.tsv")))
    cplx = _read_tsvs(sorted((qc / "complexity").glob("*.tsv")))
    frag = _read_tsvs(sorted((qc / "fragment").glob("*.tsv")))

    df = samples[["sample_id", "target", "assay", "condition", "replicate", "layout"]].copy()
    for extra in (frip, cplx, frag):
        if extra.empty:
            continue
        # The fragment table repeats `layout`; merging it in unrenamed would produce
        # layout_x/layout_y and silently break every later reference to `layout`.
        dup = [c for c in extra.columns if c in df.columns and c != "sample"]
        extra = extra.drop(columns=dup)
        df = df.merge(extra, left_on="sample_id", right_on="sample", how="left")
        df = df.drop(columns=["sample"])

    df["usable_reads"] = [
        _flagstat_mapped(qc / "flagstat" / f"{s}.txt") for s in df.sample_id
    ]
    logdir = qc.parent / "logs" / "bowtie2"
    df["alignment_rate"] = [_alignment_rate(logdir / f"{s}.log") for s in df.sample_id]

    qcfg = cfg["qc"]
    rows = []
    for _, r in df.iterrows():
        flags: list[str] = []
        status = "PASS"

        def fail(msg: str) -> None:
            nonlocal status
            flags.append(msg)
            status = "FAIL"

        def warn(msg: str) -> None:
            nonlocal status
            flags.append(msg)
            if status != "FAIL":
                status = "WARN"

        is_input = r.assay == "input"
        spec = reg.get(r.target) if not is_input else None
        broad = (not is_input) and reg.is_broad(r.target)

        # --- depth. ENCODE floors differ by peak mode; a broad domain needs far
        # more reads than a punctate TF to reach the same power.
        floor = (
            qcfg["min_usable_reads_broad"] if broad else qcfg["min_usable_reads_narrow"]
        )
        if not is_input and r.usable_reads < floor:
            warn(
                f"depth {r.usable_reads/1e6:.1f}M < {floor/1e6:.0f}M ENCODE floor "
                f"for {'broad' if broad else 'narrow'} marks (underpowered)"
            )

        # --- alignment rate. A clean human ChIP aligns at ~90-98%. Well below that
        # means most reads are not human: adapter dimer, contamination, or degraded
        # material. Deduplication and normalisation cannot recover from it.
        ar = r.get("alignment_rate")
        if pd.notna(ar):
            if ar < 0.50:
                fail(f"alignment rate {ar:.0%}: most reads are not human (contamination or adapter dimer)")
            elif ar < 0.80:
                warn(f"alignment rate {ar:.0%} < 80%")

        # --- library complexity
        nrf = r.get("NRF")
        pbc1 = r.get("PBC1")
        if pd.notna(nrf) and nrf < qcfg["min_nrf"]:
            warn(f"NRF {nrf:.2f} < {qcfg['min_nrf']} (PCR bottleneck)")
        if pd.notna(pbc1) and pbc1 < qcfg["min_pbc1"]:
            warn(f"PBC1 {pbc1:.2f} < {qcfg['min_pbc1']} (PCR bottleneck)")

        # --- enrichment, judged against this mark's own threshold
        if not is_input:
            fr = r.get("FRiP")
            min_frip = spec["qc"]["min_frip"]
            if pd.isna(fr):
                warn("FRiP not computed")
            elif fr < min_frip:
                fail(f"FRiP {fr:.3f} < {min_frip} required for {r.target}")

            npk = r.get("n_peaks")
            min_peaks = spec["qc"]["min_peaks"]
            if pd.notna(npk) and npk < min_peaks:
                warn(f"{int(npk)} peaks < {min_peaks} typical for {r.target}")

        rows.append(
            {
                "sample": r.sample_id,
                "target": r.target,
                "assay": r.assay,
                "condition": r.condition,
                "replicate": r.replicate,
                "layout": r.layout,
                "peak_mode": "-" if is_input else spec["peak_mode"],
                "usable_reads": int(r.usable_reads),
                "alignment_rate": r.get("alignment_rate"),
                "fragment_length": r.get("fragment_length"),
                "NRF": r.get("NRF"),
                "PBC1": r.get("PBC1"),
                "FRiP": r.get("FRiP"),
                "n_peaks": r.get("n_peaks"),
                "status": status,
                "flags": "; ".join(flags) if flags else "",
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(args.tsv, sep="\t", index=False)

    n_fail = int((out.status == "FAIL").sum())
    n_warn = int((out.status == "WARN").sum())

    with open(args.md, "w") as fh:
        fh.write("# QC gate\n\n")
        fh.write(
            f"{len(out)} libraries: {int((out.status=='PASS').sum())} PASS, "
            f"{n_warn} WARN, {n_fail} FAIL.\n\n"
        )
        fh.write(
            "Thresholds are per target, from `config/mark_registry.yaml`. "
            "A FRiP that fails H3K4me3 is normal for H3K9me3, so a single global "
            "cutoff would be meaningless.\n\n"
        )
        cols = [
            "sample", "target", "peak_mode", "usable_reads", "alignment_rate",
            "fragment_length", "NRF", "FRiP", "n_peaks", "status",
        ]
        fh.write("| " + " | ".join(cols) + " |\n")
        fh.write("|" + "---|" * len(cols) + "\n")
        for _, r in out.iterrows():
            vals = []
            for c in cols:
                v = r[c]
                if c == "usable_reads":
                    v = f"{int(v)/1e6:.1f}M"
                elif c in ("NRF", "FRiP") and pd.notna(v):
                    v = f"{float(v):.3f}"
                elif c == "alignment_rate" and pd.notna(v):
                    v = f"{float(v):.0%}"
                elif c == "n_peaks" and pd.notna(v):
                    v = f"{int(v):,}"
                vals.append("" if pd.isna(v) else str(v))
            fh.write("| " + " | ".join(vals) + " |\n")

        flagged = out[out.flags != ""]
        if not flagged.empty:
            fh.write("\n## Flags\n\n")
            for _, r in flagged.iterrows():
                fh.write(f"- **{r['sample']}** ({r['status']}): {r['flags']}\n")

    print(out.to_string(index=False))
    print(f"\n{n_fail} FAIL, {n_warn} WARN -> {args.tsv}")
    if n_fail:
        # Visible, but not fatal: the downstream analysis still runs and the
        # summary carries the flag forward.
        print("NOTE: FAIL libraries remain in the analysis; see the summary.", file=sys.stderr)


if __name__ == "__main__":
    main()
