#!/usr/bin/env python3
"""Write the run's analysis summary: what was found, and what to distrust.

Deliberately reports the caveats alongside the results. A summary that lists
differential peak counts without saying that one arm of the contrast is 5x
shallower than the other is worse than no summary, because it invites the reader
to believe a number the data cannot support.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from marks import MarkRegistry


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", default="results")
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--registry", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    R = Path(args.results)
    cfg = yaml.safe_load(open(args.config))
    reg = MarkRegistry(args.registry)
    samples = pd.read_csv(args.samples, sep="\t", dtype=str).fillna("")
    chip = samples[samples.assay == "chip"]
    targets = sorted(chip.target.unique())

    ref, trt = cfg["contrast"]["reference"], cfg["contrast"]["treatment"]
    fdr, lfc = cfg["differential"]["fdr"], cfg["differential"]["min_lfc"]

    gate = read_tsv(R / "qc" / "qc_gate.tsv")

    L: list[str] = []
    L.append(f"# ChIP-seq differential binding: {trt} vs {ref}\n")
    L.append(
        f"{len(samples)} libraries, {len(chip)} ChIP and {len(samples)-len(chip)} input, "
        f"across {len(targets)} targets.\n"
    )

    # ---- peak mode table
    L.append("## Peak-mode assignment\n")
    L.append(
        "Assigned from `config/mark_registry.yaml`. The mode selects the peak "
        "caller, the reproducibility test, the profile geometry, the counting "
        "window, the normalisation and whether motif analysis runs.\n"
    )
    L.append("| target | class | mode | reproducibility | normalisation | motifs |")
    L.append("|---|---|---|---|---|---|")
    for t in targets:
        s = reg.get(t)
        L.append(
            f"| {t} | {s['class']} | **{s['peak_mode']}** | {s['reproducibility']} | "
            f"{s['diffbind']['normalize']} | {'yes' if s['motifs']['enabled'] else 'no'} |"
        )
    L.append("")

    # ---- QC
    if not gate.empty:
        n_fail = int((gate.status == "FAIL").sum())
        n_warn = int((gate.status == "WARN").sum())
        L.append("## Library QC\n")
        L.append(
            f"{int((gate.status=='PASS').sum())} PASS, {n_warn} WARN, {n_fail} FAIL "
            f"(thresholds are per mark; see `results/qc/qc_gate.md`).\n"
        )
        flagged = gate[gate.flags.fillna("") != ""]
        if not flagged.empty:
            L.append("Flagged libraries:\n")
            for _, r in flagged.iterrows():
                L.append(f"- `{r['sample']}` **{r['status']}** — {r['flags']}")
            L.append("")

        # Depth imbalance between the arms of the contrast is the failure mode that
        # most often produces spurious differential binding, so check for it by name.
        L.append("### Depth balance across the contrast\n")
        L.append("| target | " + f"{ref} (M) | {trt} (M) | ratio |")
        L.append("|---|---|---|---|")
        for t in targets:
            sub = gate[(gate.target == t) & (gate.assay == "chip")]
            a = sub[sub.condition == ref].usable_reads.mean() / 1e6
            b = sub[sub.condition == trt].usable_reads.mean() / 1e6
            ratio = (max(a, b) / min(a, b)) if min(a, b) > 0 else float("nan")
            warn = " ⚠️" if ratio >= 2 else ""
            L.append(f"| {t} | {a:.1f} | {b:.1f} | {ratio:.1f}×{warn} |")
        L.append("")
        L.append(
            "A ratio above ~2x means the two arms are not equally powered. "
            "Normalisation rescales the counts but cannot restore the information "
            "that was never sequenced: the shallower arm has genuinely fewer reads "
            "per peak, so peaks that are real but weak will be detected in the deep "
            "arm and missed in the shallow one. That asymmetry looks exactly like "
            "biology and is not.\n"
        )

    # ---- differential binding
    L.append("## Differential binding\n")
    L.append(f"DESeq2, FDR < {fdr} and |log2FC| >= {lfc}, design `~ replicate + condition`.\n")
    L.append(f"| target | mode | regions tested | up in {trt} | down in {trt} | normalisation |")
    L.append("|---|---|---|---|---|---|")
    for t in targets:
        d = read_tsv(R / "differential" / t / "results.tsv")
        s = reg.get(t)
        if d.empty or "direction" not in d:
            L.append(f"| {t} | {s['peak_mode']} | — | — | — | {s['diffbind']['normalize']} |")
            continue
        up = int((d.direction == "up").sum())
        dn = int((d.direction == "down").sum())
        L.append(
            f"| {t} | {s['peak_mode']} | {len(d):,} | {up:,} | {dn:,} | "
            f"{s['diffbind']['normalize']} |"
        )
    L.append("")

    L.append(
        "Broad marks are normalised on genome-wide background bins rather than on "
        "reads-in-peaks. Median-of-ratios on peak counts assumes most peaks are "
        "unchanged; a repressive domain mark can move globally, and when it does "
        "that assumption silently divides the effect out. Where a truly global shift "
        "is expected, an exogenous spike-in is the only unbiased normaliser — "
        "background bins are a proxy.\n"
    )

    # ---- top differential regions
    for t in targets:
        d = read_tsv(R / "differential" / t / "results.tsv")
        genes = read_tsv(R / "annotation" / t / "differential_genes.tsv")
        if d.empty or "direction" not in d:
            continue
        sig = d[d.direction != "ns"]
        if sig.empty:
            continue
        L.append(f"### {t} — strongest changes\n")
        top = sig.reindex(sig.log2FoldChange.abs().sort_values(ascending=False).index).head(8)
        L.append("| region | log2FC | padj | direction |")
        L.append("|---|---|---|---|")
        for _, r in top.iterrows():
            L.append(
                f"| {r['chr']}:{int(r['start']):,}-{int(r['end']):,} | "
                f"{r['log2FoldChange']:.2f} | {r['padj']:.2e} | {r['direction']} |"
            )
        L.append("")
        if not genes.empty and "gene" in genes:
            for direction in ("up", "down"):
                g = genes[genes.direction == direction].gene.dropna().unique()
                if len(g):
                    L.append(
                        f"- **{direction}** ({len(g)} genes): "
                        + ", ".join(f"`{x}`" for x in list(g)[:15])
                        + ("…" if len(g) > 15 else "")
                    )
            L.append("")

    # ---- GO
    go_any = False
    for t in targets:
        go = read_tsv(R / "annotation" / t / "go_enrichment.tsv")
        if go.empty or "Description" not in go:
            continue
        if not go_any:
            L.append("## GO enrichment (biological process)\n")
            L.append(
                "Universe: genes near the consensus peaks for that mark, not the whole "
                "genome. Against a whole-genome universe an H3K27ac set would simply "
                "rediscover that H3K27ac marks expressed genes.\n"
            )
            go_any = True
        for direction in go.direction.unique():
            sub = go[go.direction == direction].nsmallest(5, "p.adjust")
            if sub.empty:
                continue
            L.append(f"**{t} / {direction}**\n")
            for _, r in sub.iterrows():
                L.append(f"- {r['Description']} (q = {r['p.adjust']:.1e})")
            L.append("")

    # ---- motifs
    mot_targets = [t for t in targets if reg.motifs_enabled(t)]
    skipped = [t for t in targets if not reg.motifs_enabled(t)]
    if mot_targets:
        L.append("## Motif enrichment (JASPAR / monaLisa)\n")
        L.append(
            "Foreground: peaks that changed. Background: peaks of the same mark that "
            "did not. Against a random genomic background every enhancer set is "
            "enriched for every enhancer factor, which answers nothing.\n"
        )
        for t in mot_targets:
            for direction in ("up", "down"):
                d = read_tsv(R / "motifs" / t / direction / "motif_enrichment.tsv")
                if d.empty or "negLog10Padj" not in d:
                    continue
                d = d[(d.log2enr > 0) & d.negLog10Padj.notna()]
                d = d.nlargest(5, "negLog10Padj")
                if d.empty:
                    continue
                L.append(f"**{t} / {direction}**\n")
                for _, r in d.iterrows():
                    L.append(
                        f"- {r['name']} (log2 enrichment {r['log2enr']:.2f}, "
                        f"q = {10 ** -r['negLog10Padj']:.1e})"
                    )
                L.append("")
    if skipped:
        L.append(
            f"Motif analysis was not run for {', '.join(skipped)}: these are broad "
            "marks, and scanning a multi-kb domain for 8-mers recovers its base "
            "composition rather than any bound factor.\n"
        )

    L.append("## Figures\n")
    for name, desc in [
        ("fig1_qc.pdf", "QC: FRiP against each mark's own threshold, depth, complexity, fragment size"),
        ("fig2_peak_landscape.pdf", "Peak widths by registry class — the direct test of the narrow/broad rule"),
        ("fig3_sample_relationships.pdf", "Spearman correlation and PCA on log2(ChIP/Input)"),
        ("fig4_differential.pdf", "Volcano per mark and differential peak counts"),
        ("fig5_annotation_go.pdf", "Genomic distribution and GO enrichment"),
        ("fig6_motifs.pdf", "Known-motif enrichment in differential peaks"),
        ("browser/", "Genome-browser panels at the activation loci"),
    ]:
        if (R / "figures" / name).exists():
            L.append(f"- `results/figures/{name}` — {desc}")
    L.append("")

    Path(args.out).write_text("\n".join(L) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
