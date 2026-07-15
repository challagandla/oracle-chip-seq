# Validation

This file records repository checks run on 2026-07-15. It is not evidence that
the workflow has been biologically
validated on a real ChIP-seq cohort. The GitHub Actions workflow reruns the
portable static and synthetic checks for each pushed commit and pull request.

## Automated checks

The following checks were run from the repository root:

```bash
ruff check scripts tests
ruff format --check scripts tests
python -m py_compile scripts/*.py tests/*.py
bash -n setup.sh run.sh make_config.sh prepare_references.sh
Rscript -e 'parse(file="analysis/diffbind_analysis.R"); parse(file="analysis/motif_enrichment.R")'
python -m pytest -q
git diff --check
```

Results on 2026-07-15:

- Ruff 0.15.17 reported no lint or formatting findings.
- All shell entry points and Python modules parsed successfully.
- Both R analysis scripts parsed with R 4.6.1.
- Pytest reported `30 passed` with Python 3.13.13 and Snakemake 9.22.0.
- The tests include Snakemake lint and a forced dry-run of a 135-job synthetic
  mixed narrow/broad design.
- The tracked diff contained no whitespace errors.

The synthetic design creates empty placeholder inputs. It verifies workflow
construction, declared paths, rule selection, shell-command composition, and
validation failures; it does not execute the scientific tools.

## Environment resolution

Each version-constrained environment was dry-solved without installation for
both declared Linux platforms:

```bash
MAMBA_ROOT_PREFIX=/tmp/oracle-chip-mamba-root \
  mamba create --dry-run --yes --platform linux-64 \
  --prefix /tmp/oracle-chip-runner-linux-64 --file environment.runner.yml
MAMBA_ROOT_PREFIX=/tmp/oracle-chip-mamba-root \
  mamba create --dry-run --yes --platform linux-aarch64 \
  --prefix /tmp/oracle-chip-runner-linux-aarch64 --file environment.runner.yml

MAMBA_ROOT_PREFIX=/tmp/oracle-chip-mamba-root \
  mamba create --dry-run --yes --platform linux-64 \
  --prefix /tmp/oracle-chip-tools-linux-64 --file envs/chipseq.yaml
MAMBA_ROOT_PREFIX=/tmp/oracle-chip-mamba-root \
  mamba create --dry-run --yes --platform linux-aarch64 \
  --prefix /tmp/oracle-chip-tools-linux-aarch64 --file envs/chipseq.yaml

MAMBA_ROOT_PREFIX=/tmp/oracle-chip-mamba-root \
  mamba create --dry-run --yes --platform linux-64 \
  --prefix /tmp/oracle-chip-r-linux-64 --file envs/r_analysis.yaml
MAMBA_ROOT_PREFIX=/tmp/oracle-chip-mamba-root \
  mamba create --dry-run --yes --platform linux-aarch64 \
  --prefix /tmp/oracle-chip-r-linux-aarch64 --file envs/r_analysis.yaml
```

All six original solves completed successfully on 2026-07-14. After the motif
environment reconciliation on 2026-07-15, both current rule environments were
dry-solved again for Linux x86_64 and Linux aarch64. A cross-platform solve
checks package availability and dependency compatibility; it is not native
aarch64 runtime testing. These environment files constrain direct dependency
versions but are not lockfiles, so a later solve may select newer compatible
transitive dependencies.

The reconciled R environment was also installed on Linux x86_64 with R 4.3.3,
and all 15 direct R/Bioconductor libraries passed an isolated import smoke test.
A real synthetic motif run then executed `monaLisa` 1.8.0 with three foreground
and six deterministic, GC-matched background windows. It wrote the declared
TSV, CSV, and PDF; a second run with the same seed produced byte-identical TSV
and CSV files. This is a software-path smoke test, not biological validation.

The corrected ChIP rule environment was installed on Linux x86_64 with
deepTools 3.5.5 and NumPy 1.26.4. A synthetic BigWig and BED then completed a
real `computeMatrix reference-point` run and `plotHeatmap`, producing a matrix,
plain-text matrix export, and PNG. This specifically exercises the NumPy code
path that is incompatible with NumPy 2 in deepTools 3.5.x.

## Reference endpoints

The three FASTA URLs declared in `scripts/download_references.py` returned HTTP
200 to `curl --fail --location --head` on 2026-07-14. Endpoint availability is
time-dependent; the helper additionally rejects empty or corrupt compressed
downloads before moving them into place.

## Validation boundaries

The following remain operator responsibilities:

- confirm sample identity, biological replication, antibody suitability, and
  the correct matched input for each ChIP library;
- use reference FASTA, chromosome sizes, blacklist, Bowtie2, and FastQ Screen
  indexes from the same assembly;
- choose an effective genome size appropriate for the assembly, read length,
  and mapping/filtering policy;
- inspect FastQC, FastQ Screen, alignment, library-complexity, fingerprint,
  correlation, PCA, peak, and browser-track outputs before interpretation;
- use an experimental normalization strategy such as spike-in controls when
  the biological question concerns global occupancy shifts; and
- perform a real-data end-to-end run on each production architecture before
  treating that architecture as operationally qualified.
