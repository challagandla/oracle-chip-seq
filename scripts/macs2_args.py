#!/usr/bin/env python3
"""Print the MACS2 command-line flags for a target.

Exists so the shell rules never assemble peak-caller flags themselves. The peak
mode, q-value and single-end fragment size all come from the registry, which
means "is this mark narrow or broad" is answered in exactly one place.
"""
import argparse

from marks import MarkRegistry


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", required=True)
    p.add_argument("--gsize", required=True)
    p.add_argument("--paired", type=int, default=0)
    p.add_argument("--extsize", type=int, default=None)
    p.add_argument("--registry", default="config/mark_registry.yaml")
    args = p.parse_args()

    reg = MarkRegistry(args.registry)
    print(reg.macs2_args(args.target, args.gsize, bool(args.paired), args.extsize))


if __name__ == "__main__":
    main()
