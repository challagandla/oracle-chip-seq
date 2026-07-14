#!/usr/bin/env bash
# Run the reference helper with the ChIP rule environment on PATH.
set -Eeuo pipefail

ENV_NAME="oracle-chip-runner"
ENV_SNAKEFILE="workflow/envs.smk"
REQUIRED_ENVS=("envs/chipseq.yaml")

usage() {
    cat <<'EOF'
Download and index a ChIP-seq reference using packages installed by setup.sh.

Usage:
  bash prepare_references.sh human --outdir references
  bash prepare_references.sh mouse --outdir references
  bash prepare_references.sh rat --outdir references
EOF
}

die() { printf '[ERROR] %s\n' "$1" >&2; exit 1; }
[[ "${1:-}" != "--help" && "${1:-}" != "-h" ]] || { usage; exit 0; }
[[ $# -gt 0 ]] || { usage >&2; exit 1; }

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$ROOT_DIR/.snakemake/cache}"
mkdir -p "$XDG_CACHE_HOME"

CONDA_BIN=""
for candidate in \
    "${CONDA_EXE:-}" \
    "${MINIFORGE_HOME:+${MINIFORGE_HOME%/}/bin/conda}" \
    "$(type -P conda 2>/dev/null || true)" \
    "$HOME/miniforge3/bin/conda" \
    "$HOME/miniconda3/bin/conda"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then CONDA_BIN="$candidate"; break; fi
done
[[ -n "$CONDA_BIN" ]] || die "Conda is missing; run: bash setup.sh"
"$CONDA_BIN" run --name "$ENV_NAME" python -c 'import sys' >/dev/null 2>&1 || \
    die "Runner $ENV_NAME is missing; run: bash setup.sh"

version="$("$CONDA_BIN" run --name "$ENV_NAME" snakemake --version | tail -n 1)"
major="${version%%.*}"
[[ "$major" =~ ^[0-9]+$ ]] || die "Could not parse Snakemake version: $version"
if ((major >= 8)); then deployment=(--software-deployment-method conda); else deployment=(--use-conda); fi

listing="$("$CONDA_BIN" run --name "$ENV_NAME" snakemake \
    --snakefile "$ENV_SNAKEFILE" --cores 1 "${deployment[@]}" --list-conda-envs 2>&1)" || \
    die "Could not locate rule environments; run: bash setup.sh"

PATH_PREFIX=""
for spec in "${REQUIRED_ENVS[@]}"; do
    location="$(awk -v wanted="$spec" '$1 == wanted {print $NF; exit}' <<< "$listing")"
    [[ -n "$location" ]] || die "Environment $spec is not registered; run: bash setup.sh"
    prefix="$ROOT_DIR/$location"
    [[ -d "$prefix/conda-meta" ]] || die "Environment $spec is not installed; run: bash setup.sh"
    PATH_PREFIX="${PATH_PREFIX:+$PATH_PREFIX:}$prefix/bin"
done

export PATH="$PATH_PREFIX:$PATH"
for tool in python wget gunzip gzip bowtie2-build samtools; do
    command -v "$tool" >/dev/null 2>&1 || die "$tool is unavailable; rerun: bash setup.sh"
done

exec env -u PYTHONPATH -u R_LIBS -u R_LIBS_USER \
    python scripts/download_references.py "$@"
