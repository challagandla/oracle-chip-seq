#!/usr/bin/env python3
"""Resolve a ChIP target to its peak mode and downstream parameters.

This is the only place that decides "narrow or broad". The Snakefile, the QC
gate, the DiffBind driver, the profile plots and the motif step all resolve
their settings through here, so a target cannot be called narrow by the peak
caller and treated as broad by the differential analysis.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import yaml

DEFAULT_REGISTRY = "config/mark_registry.yaml"


def _deep_merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


class MarkRegistry:
    def __init__(self, path: str | Path = DEFAULT_REGISTRY):
        self.path = Path(path)
        with open(self.path) as fh:
            self.raw = yaml.safe_load(fh)
        self.defaults = self.raw["defaults"]
        self.marks = self.raw.get("marks", {}) or {}
        self.peak_modes = self.raw["peak_modes"]
        self.aliases = {k.upper(): v for k, v in (self.raw.get("aliases") or {}).items()}
        # Case-insensitive lookup of canonical names (H3K27AC -> H3K27ac).
        self._canonical = {name.upper(): name for name in self.marks}

    def canonical(self, target: str) -> str | None:
        key = target.strip().upper()
        if key in self._canonical:
            return self._canonical[key]
        if key in self.aliases:
            return self.aliases[key]
        return None

    def is_known(self, target: str) -> bool:
        return self.canonical(target) is not None

    def get(self, target: str) -> dict[str, Any]:
        """Full resolved parameter set for a target. Unknown targets fall back to
        the defaults (narrow), which is deliberate: an unrecognised antibody is far
        more likely to be a TF than a domain-forming histone mark, and mis-calling a
        broad mark as narrow loses breadth, whereas mis-calling a punctate mark as
        broad silently fuses neighbouring regulatory elements."""
        name = self.canonical(target)
        entry = self.marks.get(name, {}) if name else {}
        resolved = _deep_merge(self.defaults, entry)
        resolved["target"] = name or target
        resolved["known"] = name is not None
        mode = resolved["peak_mode"]
        if mode not in self.peak_modes:
            raise ValueError(f"{target}: unknown peak_mode '{mode}'")
        resolved["peak_ext"] = self.peak_modes[mode]["peak_ext"]
        resolved["macs2_flag"] = self.peak_modes[mode]["macs2_flag"]
        return resolved

    # -- convenience accessors used throughout the workflow ------------------

    def peak_mode(self, target: str) -> str:
        return self.get(target)["peak_mode"]

    def peak_ext(self, target: str) -> str:
        return self.get(target)["peak_ext"]

    def is_broad(self, target: str) -> bool:
        return self.peak_mode(target) == "broad"

    def uses_idr(self, target: str) -> bool:
        return self.get(target)["reproducibility"] == "idr"

    def motifs_enabled(self, target: str) -> bool:
        return bool(self.get(target)["motifs"]["enabled"])

    def macs2_args(self, target: str, gsize: str, is_paired: bool, extsize: int | None = None) -> str:
        """Assemble the MACS2 command line for this target.

        Two conditions are injected here:
          * peak mode  -> --broad (+ --broad-cutoff) or narrow with a tighter q
          * library layout -> BAMPE lets MACS2 use the true fragment length;
            single-end has no fragment, so we supply the cross-correlation
            fragment estimate via --nomodel --extsize.
        """
        cfg = self.get(target)
        macs2 = cfg["macs2"]
        args = [f"--gsize {gsize}", f"--qvalue {macs2.get('qvalue', 0.05)}"]

        if cfg["peak_mode"] == "broad":
            args.append("--broad")
            args.append(f"--broad-cutoff {macs2.get('broad_cutoff', 0.1)}")

        if is_paired:
            args.append("--format BAMPE")
        else:
            args.append("--format BAM")
            if extsize:
                args.append(f"--nomodel --extsize {int(extsize)}")

        extra = (macs2.get("extra") or "").strip()
        if extra:
            # --call-summits is incompatible with --broad; MACS2 errors out.
            if cfg["peak_mode"] == "broad":
                extra = extra.replace("--call-summits", "").strip()
            if extra:
                args.append(extra)
        return " ".join(args)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("target", nargs="?", help="Target to resolve; omit to dump the whole table")
    p.add_argument("--registry", default=DEFAULT_REGISTRY)
    p.add_argument("--field", help="Print a single top-level field")
    args = p.parse_args()

    reg = MarkRegistry(args.registry)
    if not args.target:
        rows = [reg.get(m) for m in reg.marks]
        width = max(len(r["target"]) for r in rows)
        print(f'{"target".ljust(width)}  {"mode":7} {"class":15} {"repro":8} {"profile":16} {"summits":8} {"norm":5} motifs')
        for r in sorted(rows, key=lambda r: (r["peak_mode"], r["target"])):
            print(
                f'{r["target"].ljust(width)}  {r["peak_mode"]:7} {r["class"]:15} '
                f'{r["reproducibility"]:8} {r["profile"]:16} '
                f'{str(r["diffbind"]["summits"]):8} {r["diffbind"]["normalize"]:5} '
                f'{"yes" if r["motifs"]["enabled"] else "no"}'
            )
        return

    resolved = reg.get(args.target)
    if args.field:
        print(resolved[args.field])
    else:
        print(json.dumps(resolved, indent=2, default=str))


if __name__ == "__main__":
    main()
