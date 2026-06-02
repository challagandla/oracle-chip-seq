#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Generate a Snakemake HTML report for the current workflow.")
    parser.add_argument("output", help="Path to the generated report HTML")
    parser.add_argument("--snakefile", default="Snakefile", help="Snakefile path")
    parser.add_argument("--configfile", default="config.yaml", help="Config file path")
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "snakemake",
        "--snakefile",
        args.snakefile,
        "--configfile",
        args.configfile,
        "--report",
        str(output),
        "--nolock",
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"Created report: {output}")


if __name__ == "__main__":
    main()
