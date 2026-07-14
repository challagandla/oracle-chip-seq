#!/usr/bin/env bash
# Run the manifest converter in the declared runner environment.
set -Eeuo pipefail

ENV_NAME="oracle-chip-runner"

usage() {
    cat <<'EOF'
Build config.yaml from a tab-separated ChIP-seq sample manifest.

Usage:
  bash make_config.sh sample_manifest.tsv \
    --species human \
    --genome /path/to/hg38.fa \
    --chrom-sizes /path/to/hg38.chrom.sizes \
    --bt2-index /path/to/hg38 \
    --black-list /path/to/hg38-blacklist.v2.bed \
    --effective-genome-size 2862010428 \
    --fastq-screen-conf config/fastq_screen.conf \
    --output config.yaml

Run bash setup.sh before using this command. All arguments after the script
name are passed to scripts/manifest_to_config.py.
EOF
}

die() { printf '[ERROR] %s\n' "$1" >&2; exit 1; }
[[ $# -gt 0 ]] || { usage >&2; exit 1; }

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

CONDA_BIN=""
for candidate in \
    "${CONDA_EXE:-}" \
    "${MINIFORGE_HOME:+${MINIFORGE_HOME%/}/bin/conda}" \
    "$(type -P conda 2>/dev/null || true)" \
    "$HOME/miniforge3/bin/conda" \
    "$HOME/miniconda3/bin/conda"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then
        CONDA_BIN="$candidate"
        break
    fi
done

[[ -n "$CONDA_BIN" ]] || die "Conda is missing; run: bash setup.sh"
"$CONDA_BIN" run --name "$ENV_NAME" python -c 'import yaml' >/dev/null 2>&1 || \
    die "Runner $ENV_NAME is missing; run: bash setup.sh"

exec env -u PYTHONPATH -u R_LIBS -u R_LIBS_USER \
    R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null \
    "$CONDA_BIN" run --no-capture-output --name "$ENV_NAME" \
    python scripts/manifest_to_config.py "$@"
