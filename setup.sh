#!/usr/bin/env bash
# Reproducible, beginner-friendly ORACLE pipeline installer.
set -Eeuo pipefail

PIPELINE_NAME="oracle-chip-seq"
ENV_NAME="oracle-chip-runner"
RUNNER_FILE="environment.runner.yml"
ENV_SNAKEFILE="workflow/envs.smk"
EXPECTED_ENV_COUNT=2
MINIFORGE_VERSION="${MINIFORGE_VERSION:-26.3.2-2}"
MIN_FREE_GB="${MIN_FREE_GB:-20}"

CHECK_ONLY=0
INSTALL_RULE_ENVS=1

usage() {
    cat <<EOF
Install and verify ${PIPELINE_NAME}.

Usage:
  bash setup.sh                 Install/update the runner and every rule environment
  bash setup.sh --runner-only   Install/update only the small runner environment
  bash setup.sh --check         Check the existing installation without changing it
  bash setup.sh --help          Show this help

The installer uses an existing Conda/Miniforge installation when available. If
Conda is absent, it installs a pinned, checksummed Miniforge release under
\${MINIFORGE_HOME:-\$HOME/miniforge3}. Set MINIFORGE_HOME to use a custom
prefix. It never modifies your shell startup files. The environment
specifications are solver-validated for Linux x86_64 and Linux aarch64;
package smoke checks run on the current host during installation.
EOF
}

log() { printf '[%s] %s\n' "$1" "$2"; }
die() { log ERROR "$1" >&2; exit 1; }

while (($#)); do
    case "$1" in
        --check) CHECK_ONLY=1 ;;
        --runner-only) INSTALL_RULE_ENVS=0 ;;
        --help|-h) usage; exit 0 ;;
        *) die "Unknown option: $1 (use --help)" ;;
    esac
    shift
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$ROOT_DIR/.snakemake/cache}"
mkdir -p "$XDG_CACHE_HOME"
[[ -f "$RUNNER_FILE" ]] || die "Missing $RUNNER_FILE"
[[ -f "$ENV_SNAKEFILE" ]] || die "Missing $ENV_SNAKEFILE"

SYSTEM="$(uname -s)"
ARCHITECTURE="$(uname -m)"
[[ "$SYSTEM" == "Linux" ]] || \
    die "The declared rule environments support Linux only; found $SYSTEM"
case "$ARCHITECTURE" in
    x86_64) ;;
    aarch64|arm64) ARCHITECTURE="aarch64" ;;
    *) die "The declared rule environments support x86_64 or aarch64; found $ARCHITECTURE" ;;
esac

available_kb="$(df -Pk "$ROOT_DIR" | awk 'NR == 2 {print $4}')"
if [[ "$available_kb" =~ ^[0-9]+$ ]]; then
    available_gb=$((available_kb / 1024 / 1024))
    if ((available_gb < MIN_FREE_GB)); then
        log WARN "Only ${available_gb} GB is free; at least ${MIN_FREE_GB} GB is recommended"
    else
        log INFO "Free disk space: ${available_gb} GB"
    fi
fi

CONDA_BIN=""
SOLVER_BIN=""
MINIFORGE_TEMP_DIR=""

cleanup() {
    if [[ -n "$MINIFORGE_TEMP_DIR" && -d "$MINIFORGE_TEMP_DIR" ]]; then
        rm -rf -- "$MINIFORGE_TEMP_DIR"
    fi
}
trap cleanup EXIT

locate_conda() {
    local candidate
    for candidate in \
        "${CONDA_EXE:-}" \
        "${MINIFORGE_HOME:+${MINIFORGE_HOME%/}/bin/conda}" \
        "$(type -P conda 2>/dev/null || true)" \
        "$HOME/miniforge3/bin/conda" \
        "$HOME/miniconda3/bin/conda"; do
        if [[ -n "$candidate" && -x "$candidate" ]]; then
            CONDA_BIN="$candidate"
            return 0
        fi
    done
    return 1
}

download_file() {
    local url="$1" destination="$2"
    if command -v curl >/dev/null 2>&1; then
        curl --fail --location --silent --show-error "$url" --output "$destination"
    elif command -v wget >/dev/null 2>&1; then
        wget --quiet "$url" --output-document "$destination"
    else
        die "curl or wget is required to download Miniforge"
    fi
}

install_miniforge() {
    local asset base_url temp_dir expected actual prefix
    asset="Miniforge3-${MINIFORGE_VERSION}-Linux-${ARCHITECTURE}.sh"
    base_url="https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}"
    temp_dir="$(mktemp -d)"
    MINIFORGE_TEMP_DIR="$temp_dir"
    log INFO "Downloading checksummed Miniforge ${MINIFORGE_VERSION}"
    download_file "$base_url/$asset" "$temp_dir/$asset"
    download_file "$base_url/$asset.sha256" "$temp_dir/$asset.sha256"
    expected="$(awk 'NR == 1 {print $1}' "$temp_dir/$asset.sha256")"
    if command -v sha256sum >/dev/null 2>&1; then
        actual="$(sha256sum "$temp_dir/$asset" | awk '{print $1}')"
    elif command -v shasum >/dev/null 2>&1; then
        actual="$(shasum -a 256 "$temp_dir/$asset" | awk '{print $1}')"
    else
        die "sha256sum or shasum is required to verify Miniforge"
    fi
    [[ -n "$expected" && "$actual" == "$expected" ]] || die "Miniforge checksum verification failed"

    prefix="${MINIFORGE_HOME:-$HOME/miniforge3}"
    [[ ! -e "$prefix" ]] || die "$prefix already exists but no usable conda was found"
    bash "$temp_dir/$asset" -b -p "$prefix"
    CONDA_BIN="$prefix/bin/conda"
    log OK "Installed Miniforge at $prefix"
}

choose_solver() {
    local conda_dir candidate
    conda_dir="$(dirname "$CONDA_BIN")"
    for candidate in "$conda_dir/mamba" "$(type -P mamba 2>/dev/null || true)"; do
        if [[ -n "$candidate" && -x "$candidate" ]]; then
            SOLVER_BIN="$candidate"
            return
        fi
    done
    SOLVER_BIN="$CONDA_BIN"
}

env_exists() {
    "$CONDA_BIN" run --name "$ENV_NAME" python -c 'import sys' >/dev/null 2>&1
}

clean_run() {
    env -u R_LIBS -u R_LIBS_USER -u PYTHONPATH \
        CONDA_ALWAYS_YES=true \
        R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null \
        "$@"
}

set_deployment_args() {
    local version major
    version="$(clean_run "$CONDA_BIN" run --name "$ENV_NAME" snakemake --version | tail -n 1)"
    major="${version%%.*}"
    [[ "$major" =~ ^[0-9]+$ ]] || die "Could not parse Snakemake version: $version"
    SNAKEMAKE_VERSION="$version"
    if ((major >= 8)); then
        DEPLOYMENT_ARGS=(--software-deployment-method conda)
    else
        DEPLOYMENT_ARGS=(--use-conda)
    fi
}

verify_runner() {
    clean_run "$CONDA_BIN" run --name "$ENV_NAME" python -c \
        'import pandas, yaml; print("Python imports: pandas + PyYAML OK")'
    set_deployment_args
    clean_run "$CONDA_BIN" run --name "$ENV_NAME" dot -V 2>&1 | head -n 1
    log OK "Runner $ENV_NAME: Snakemake $SNAKEMAKE_VERSION"
}

list_rule_envs() {
    local listing listed installed spec env_path
    listing="$(clean_run "$CONDA_BIN" run --name "$ENV_NAME" \
        snakemake --snakefile "$ENV_SNAKEFILE" --cores 1 \
        "${DEPLOYMENT_ARGS[@]}" --list-conda-envs 2>&1)" || {
            printf '%s\n' "$listing" >&2
            die "Could not inspect rule environments"
        }
    listed="$(printf '%s\n' "$listing" | grep -Ec '\.ya?ml([[:space:]]|$)' || true)"
    installed=0
    while read -r spec env_path _; do
        if [[ "$spec" =~ \.ya?ml$ && -n "${env_path:-}" && -d "$env_path/conda-meta" ]]; then
            installed=$((installed + 1))
        fi
    done <<< "$listing"
    printf '%s\n' "$listing"
    ((listed == EXPECTED_ENV_COUNT)) || die "Expected $EXPECTED_ENV_COUNT rule specs, found $listed"
    log INFO "Rule environments present: $installed/$listed"
    [[ "$installed" -eq "$listed" ]]
}

run_rule_env_smoke_checks() {
    clean_run "$CONDA_BIN" run --no-capture-output --name "$ENV_NAME" \
        snakemake --snakefile "$ENV_SNAKEFILE" --cores 1 \
        "${DEPLOYMENT_ARGS[@]}" --forceall --rerun-incomplete
}

if ! locate_conda; then
    if ((CHECK_ONLY)); then
        die "Conda is not installed; run: bash setup.sh"
    fi
    install_miniforge
fi
choose_solver
log INFO "Conda: $CONDA_BIN"
log INFO "Solver: $SOLVER_BIN"

if ((CHECK_ONLY)); then
    env_exists || die "Runner $ENV_NAME is not installed; run: bash setup.sh"
    verify_runner
    list_rule_envs || die "One or more rule environments are missing; run: bash setup.sh"
    log INFO "Running package smoke checks inside every rule environment"
    run_rule_env_smoke_checks
    log OK "${PIPELINE_NAME} installation is complete"
    exit 0
fi

if env_exists; then
    log INFO "Updating runner environment $ENV_NAME"
    clean_run "$SOLVER_BIN" env update --name "$ENV_NAME" --file "$RUNNER_FILE" --prune
else
    log INFO "Creating runner environment $ENV_NAME"
    clean_run "$SOLVER_BIN" env create --yes --name "$ENV_NAME" --file "$RUNNER_FILE"
fi
verify_runner

if ((INSTALL_RULE_ENVS)); then
    log INFO "Installing $EXPECTED_ENV_COUNT Snakemake rule environments"
    clean_run "$CONDA_BIN" run --no-capture-output --name "$ENV_NAME" \
        snakemake --snakefile "$ENV_SNAKEFILE" --cores 1 \
        "${DEPLOYMENT_ARGS[@]}" --conda-create-envs-only
    list_rule_envs || die "Rule-environment verification failed"
    log INFO "Running package smoke checks inside every rule environment"
    run_rule_env_smoke_checks
else
    log INFO "Skipped rule environments (--runner-only); the first run will create them automatically"
fi

log OK "${PIPELINE_NAME} is ready"
printf '\nNext steps:\n  bash run.sh --dry-run\n  bash run.sh --cores 8\n'
