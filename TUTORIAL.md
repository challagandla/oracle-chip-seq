# Oracle ChIP-seq: complete beginner tutorial

This tutorial takes a paired-end ChIP-seq experiment from raw FASTQ files to
quality-control reports, aligned reads, peaks, signal tracks, heatmaps,
profiles, differential binding results, and motif enrichment.

Run every command from the root of the `oracle-chip-seq` repository unless the
tutorial says otherwise.

Use the guide in order for a first run:

- Sections 1-5 explain the experiment, pipeline, files, and tools.
- Sections 6-12 install the software and prepare the data and configuration.
- Sections 13-15 run the workflow and show where the files are written.
- Sections 16-23 explain how to review the results.
- Sections 24-29 cover limits, troubleshooting, a checklist, and a glossary.

## 1. Before you start

Use this workflow only when these statements describe the project:

- Your libraries are ChIP-seq libraries, not RNA-seq or ATAC-seq libraries.
- Every library was sequenced as paired-end reads.
- Your FASTQ files are gzip-compressed and end in `.fastq.gz` or `.fq.gz`.
- Each factor has exactly two biological conditions.
- Each factor has at least two biological ChIP replicates in each condition.
- All samples for one factor come from one tissue or cell-type label.
- You are running Linux on an x86_64 or aarch64 computer.

The validator enforces paired compressed FASTQ paths, two conditions, replicate
counts, one tissue label per factor, and matching ChIP/input condition labels.
You remain responsible for confirming the assay type, that declared replicates
are truly biological replicates, and that every input is appropriate for its
ChIP library. Software cannot verify laboratory provenance or biological
suitability. A donor-paired, blocked, batch-aware, multi-tissue, longitudinal,
or three-condition experiment needs a different DiffBind design.

Allow at least 20 GB of free space before installation. A real ChIP-seq project
usually needs much more space for FASTQs, the reference genome, BAM files, and
results. Keep at least as much free space as the combined size of the FASTQs,
and preferably several times that amount.

## 2. What ChIP-seq measures

ChIP-seq is used to find genomic DNA associated with a protein or histone
modification. A typical experiment has four main laboratory steps:

1. DNA and proteins are cross-linked or otherwise preserved together.
2. Chromatin is broken into fragments.
3. An antibody enriches fragments associated with the target.
4. The enriched DNA is sequenced.

The resulting reads do not directly say that a protein binds one exact base.
They show where fragments became enriched after immunoprecipitation.

### ChIP and input libraries

A ChIP library contains the antibody-enriched fragments. An input library is
prepared from chromatin before immunoprecipitation. The input captures effects
such as chromatin accessibility, fragmentation bias, copy-number differences,
and sequencing bias. The pipeline uses the input both for MACS3 peak calling
and for log2 ChIP/input signal tracks.

Match every ChIP sample to an input from the same biological material,
condition, preparation, and batch as closely as possible. The configuration
can reuse one input for several ChIP samples, but reuse is only sound when that
input is genuinely matched to all of them.

An input control and a biological control condition are different things:

| Term | Meaning |
| --- | --- |
| Input control | DNA from the same chromatin preparation without the target-specific pull-down. It measures technical background. |
| Control condition | The biological group used as the reference for a comparison, such as untreated cells. It usually has both ChIP and input libraries. |
| Treated condition | The biological group being compared with the control condition. It also usually has both ChIP and input libraries. |

### Biological replicates

A biological replicate comes from an independently prepared biological sample.
Splitting one library across two sequencing lanes creates technical replicates,
not biological replicates. Technical lanes should normally be combined before
they are entered as biological replicates.

This workflow requires at least two biological ChIP replicates per condition.
Three or more are preferable when the study design and budget allow it.

A minimum two-condition example contains four ChIP libraries and their matched
inputs:

| Condition | ChIP library | Biological replicate | Matched input |
| --- | --- | ---: | --- |
| control | H3K27ac control rep1 | 1 | Input control rep1 |
| control | H3K27ac control rep2 | 2 | Input control rep2 |
| treated | H3K27ac treated rep1 | 1 | Input treated rep1 |
| treated | H3K27ac treated rep2 | 2 | Input treated rep2 |

This table describes eight sequencing libraries: the four ChIP libraries in
the second column and the four input libraries in the last column.

### Narrow and broad enrichment

Some targets produce compact peaks. Most sequence-specific transcription
factors are treated as narrow, and punctate histone signals such as H3K27ac are
usually analyzed as narrow peaks. Other marks form domains that cover much
larger regions; H3K27me3 and H3K9me3 are common broad examples.

The pipeline never guesses the geometry from the factor name. Every ChIP
sample must set `peak_mode` to either `narrow` or `broad`, and all samples for
one factor must use the same mode.

Ask the person who designed the experiment or check the accepted analysis
practice for the target if you are unsure. Choosing a mode only because a
target is a histone mark is not sufficient.

## 3. What the pipeline does

The workflow follows this path:

```text
paired FASTQ files
        |
        +--> FastQC and FastQ Screen on raw reads
        |
        +--> Trim Galore
                |
                +--> FastQC on trimmed reads
                |
                +--> Bowtie2 alignment
                        |
                        +--> samtools filtering, duplicate removal and QC
                                |
                                +--> MACS3 peaks against matched input
                                |       |
                                |       +--> blacklist filtering
                                |       +--> replicate-supported peaks
                                |
                                +--> RPGC BigWig tracks
                                |       |
                                |       +--> correlation and PCA
                                |
                                +--> deepTools fingerprint
                                |
                                +--> log2 ChIP/input BigWig tracks
                                        |
                                        +--> deepTools heatmaps and profiles

replicate peak files + BAM files --> DiffBind differential binding
narrow consensus peaks          --> monaLisa/JASPAR motif enrichment
all main outputs                 --> MultiQC
```

Snakemake controls this process. It works out which files are needed, runs the
steps in the right order, and resumes from completed files after an
interruption.

## 4. Repository layout

The important files in a fresh checkout are:

```text
oracle-chip-seq/
├── Snakefile                         main workflow
├── config.yaml                       configuration used by default
├── config.sample.yaml                clean example configuration
├── sample_manifest.tsv               example sample table
├── setup.sh                          software installer and checker
├── run.sh                            activation-free workflow runner
├── make_config.sh                    activation-free manifest converter
├── prepare_references.sh             reference download/index wrapper
├── environment.runner.yml            Snakemake controller environment
├── envs/
│   ├── chipseq.yaml                  command-line analysis tools
│   └── r_analysis.yaml               R and DiffBind tools
├── config/
│   └── fastq_screen.conf.example     contamination-screen template
├── scripts/
│   ├── manifest_to_config.py         converts the sample table to YAML
│   ├── build_sample_sheets.py        builds the DiffBind sample sheet
│   └── download_references.py        reference helper used by the wrapper
├── analysis/
│   ├── diffbind_analysis.R           differential-binding analysis
│   └── motif_enrichment.R            motif enrichment, table, and figure
├── docs/
│   └── ANALYSIS_AND_VISUALIZATION.md detailed method notes
├── VALIDATION.md                     repository checks and validation limits
└── tests/                             workflow and configuration tests
```

You will normally add these project directories:

```text
data/raw/          paired FASTQ files
references/        genome, Bowtie2 index, chromosome sizes and blacklist
results/           files produced by the workflow
```

Do not put two unrelated analyses into the same `results/` directory. Use a
separate checkout or archive the first result directory before changing the
reference assembly, conditions, or sample design.

## 5. Tool map

You do not need to run these tools by hand. Snakemake runs each one in the
correct software environment.

| Tool | Purpose in this pipeline |
| --- | --- |
| Snakemake | Connects all inputs, commands, and outputs into one reproducible workflow. |
| Conda/Miniforge | Installs version-constrained, isolated software environments. |
| FastQC | Reports read quality, adapter content, duplication, and other FASTQ checks. |
| FastQ Screen | Maps a subset of reads to selected databases to look for unexpected sources. |
| Trim Galore | Removes adapters and low-quality sequence from paired reads. |
| Bowtie2 | Aligns the trimmed paired reads to the selected reference genome. |
| samtools | Filters alignments, repairs mate information, removes duplicates, indexes BAMs, and reports alignment counts. |
| MACS3 | Calls enriched narrow peaks or broad domains using the matched input. |
| bedtools | Removes blacklisted regions and combines replicate-supported intervals. |
| deepTools | Creates normalized BigWigs, heatmaps, average profiles, fingerprints, correlations, and PCA plots. |
| DiffBind | Counts reads over peak regions and tests for differential binding between the two conditions. |
| monaLisa + JASPAR | Tests known vertebrate motifs in narrow consensus peaks against deterministic, blacklist-excluded genomic background windows. |
| MultiQC | Collects many quality-control reports into one HTML page. |

## 6. Get the repository

If you do not already have it:

```bash
git clone https://github.com/challagandla/oracle-chip-seq.git
cd oracle-chip-seq
```

If the repository is already on your computer, open a terminal and change to
its root directory. Confirm your location:

```bash
pwd
ls
```

You should see `Snakefile`, `config.yaml`, `setup.sh`, and `run.sh`.

## 7. Install the software

Run:

```bash
bash setup.sh
```

The installer looks for an existing Conda installation. If none is available,
it downloads a pinned and checksum-verified Miniforge installer and installs it
under `$HOME/miniforge3`. It does not change your shell startup files.

The installer then creates:

- `oracle-chip-runner`, which contains Snakemake and the small Python tools
  needed to control the workflow;
- a ChIP-seq environment from `envs/chipseq.yaml`; and
- an R environment from `envs/r_analysis.yaml`.

The first installation can take tens of minutes depending on the computer,
network, and package cache. Interrupted environment creation can normally be
resumed by running `bash setup.sh` again. An interrupted Miniforge installer may
leave an incomplete prefix; see the troubleshooting section before removing it.

No `conda activate` command is required. `run.sh` locates the installation and
runs Snakemake in the correct environment.

Check an existing installation with:

```bash
bash setup.sh --check
```

Install only the small runner, leaving the scientific environments to be built
on the first run, with:

```bash
bash setup.sh --runner-only
```

For a custom Miniforge location, use the same value for setup and later runs:

```bash
export MINIFORGE_HOME=/opt/miniforge3
bash setup.sh
bash run.sh --dry-run
```

## 8. Organize the FASTQ files

Create a data directory if needed:

```bash
mkdir -p data/raw
```

Copy or link each pair of FASTQ files into it. A simple naming scheme is:

```text
data/raw/Input_control_rep1_R1.fastq.gz
data/raw/Input_control_rep1_R2.fastq.gz
data/raw/H3K27ac_control_rep1_R1.fastq.gz
data/raw/H3K27ac_control_rep1_R2.fastq.gz
```

The names on disk do not have to follow this pattern. The paths in the
configuration are what matter. `R1` must be listed first and `R2` second.

Check that a file exists and is not empty:

```bash
test -s data/raw/H3K27ac_control_rep1_R1.fastq.gz && echo "file is present"
```

Check gzip integrity without unpacking the file:

```bash
gzip -t data/raw/H3K27ac_control_rep1_R1.fastq.gz
```

No output means the gzip check passed. Repeat this check for all FASTQs. Make
sure paired files contain the same library and that they were not swapped
between samples.

## 9. Prepare the reference genome

All reference files must describe the same assembly. For example, do not mix an
hg19 FASTA with an hg38 blacklist or Bowtie2 index.

The workflow needs:

1. a genome FASTA;
2. a samtools FASTA index beside it (`genome.fa.fai`);
3. a chromosome-size table;
4. a complete Bowtie2 index prefix;
5. an assembly-matched blacklist BED;
6. a MACS3 genome size; and
7. a deepTools effective genome size written as a positive integer.

### Use the reference helper

After `bash setup.sh` has completed, the helper can download a primary-assembly
FASTA, create its samtools index, build a Bowtie2 index, and create chromosome
sizes:

```bash
bash prepare_references.sh human --outdir references
```

Replace `human` with `mouse` or `rat` when appropriate.

For human, the command prints paths similar to:

```text
Genome: references/human/hg38.fa
FASTA index: references/human/hg38.fa.fai
Bowtie2 index prefix: references/human/hg38
Chrom sizes: references/human/hg38.chrom.sizes
```

If you use an existing FASTA instead of the helper, activate an environment
that provides `samtools` and create the required index once:

```bash
GENOME=/absolute/path/to/genome.fa
samtools faidx "$GENOME"
```

This must create exactly `$GENOME.fai`. On a shared read-only reference, ask
the reference administrator to provide that adjacent index.

Building a whole-genome Bowtie2 index can take a long time. The helper safely
reuses a complete download, but it rebuilds the index when it is run again, so
do not call it repeatedly without a reason.

The helper does not download a blacklist. Obtain a blacklist made for the exact
assembly from an authoritative source such as ENCODE. Save it under
`references/` and record its path. Check that chromosome names agree between
the FASTA-derived chromosome-size file and the blacklist:

```bash
cut -f1 references/human/hg38.chrom.sizes | head
cut -f1 references/human/hg38-blacklist.v2.bed | head
```

Both should use the same naming style, such as `chr1` in both files.

### Use an existing institutional reference

You may use a reference maintained by your laboratory or computing center. In
that case, find the Bowtie2 prefix rather than one individual index file. For a
prefix such as `/refs/hg38/bowtie2/hg38`, the directory should contain either
all six `.bt2` files or all six `.bt2l` files.

Check the prefix with:

```bash
ls /refs/hg38/bowtie2/hg38*.bt2*
```

The workflow stops if it finds a partial index.

### Genome-size values

`gsize` is passed to MACS3. `effective_genome_size` is passed to deepTools for
RPGC normalization. They are related but are not the same setting.

The manifest converter supplies the assembly label and MACS3 `gsize` defaults:

| Species | Assembly label | MACS3 `gsize` |
| --- | --- | --- |
| human | hg38 | `2.7e9` |
| mouse | mm39 | `1.87e9` |
| rat | rn6 | `2.53e9` |

It deliberately does not guess the deepTools value. Pass
`--effective-genome-size` explicitly because the appropriate value depends on
the assembly, read length, and mapping/filtering policy. This workflow removes
low-MAPQ alignments, so use a read-length-specific uniquely mappable value, not
the non-N genome length intended for analyses that retain multimappers. The
checked human and mouse examples use the published deepTools 150-bp values; the rat
example fails closed until an appropriate value is supplied. In YAML, write the
value as an unquoted positive integer, not text such as
`effective_genome_size: "2.7e9"`. See the
[deepTools 3.5.6 effective-genome-size guidance](https://deeptools.readthedocs.io/en/3.5.6/content/feature/effectiveGenomeSize.html)
for the calculation methods and published examples.

## 10. Configure FastQ Screen

FastQ Screen is part of the default workflow, so its configuration must be
valid even if contamination screening is not your main analysis goal.

Copy the template:

```bash
cp config/fastq_screen.conf.example config/fastq_screen.conf
```

Open `config/fastq_screen.conf` in a text editor. Each `DATABASE` line must end
with a Bowtie2 index prefix, not an individual `.bt2` file. Remove databases
you do not have rather than leaving `/path/to/...` placeholders.

A minimal human example is:

```text
BOWTIE2 bowtie2

DATABASE Human /absolute/path/to/references/human/hg38
```

This checks mapping to the host genome. To look for contamination, add real
indexes for relevant sources such as PhiX, mouse, yeast, or bacteria. Only add
databases that are meaningful for the experiment and available on disk.

The path in `config.yaml` should then be:

```yaml
contamination:
  fastq_screen_conf: "config/fastq_screen.conf"
  subset: 100000
```

The default checks a subset of 100,000 reads from each FASTQ. This is usually
enough for a practical screen without aligning the complete file to every
database.

## 11. Describe the samples

There are two ways to create `config.yaml`:

1. edit `sample_manifest.tsv` and convert it; or
2. edit `config.sample.yaml` directly.

The manifest route is easier to check in a spreadsheet and is recommended for
a first project.

### Option A: use the sample manifest

Keep the supplied example and make a project copy:

```bash
cp sample_manifest.tsv sample_manifest.project.tsv
```

Open `sample_manifest.project.tsv` in a plain-text editor or spreadsheet
program. It is tab-separated. Keep the header unchanged. If a spreadsheet
program asks for a separator, choose tab, and save the file as tab-separated
text rather than as an Excel workbook.

| Column | Meaning |
| --- | --- |
| `sample_id` | Unique filesystem-safe name. Use letters, numbers, `.`, `_`, or `-`; do not use spaces. |
| `type` | `control` for input libraries or `chip` for ChIP libraries. |
| `condition` | Biological condition, for example `control` or `treated`. |
| `replicate` | Positive biological replicate number within the condition. |
| `fastq_r1` | Path to paired-end read 1. |
| `fastq_r2` | Path to paired-end read 2. |
| `control` | For a ChIP row, the exact `sample_id` of its matched input. Leave blank for an input row. |
| `factor` | Target or histone mark, such as `CTCF` or `H3K27ac`. Leave blank for an input row. |
| `tissue` | One tissue or cell-type label for the factor. |
| `peak_mode` | `narrow` or `broad` for a ChIP row. Leave blank for an input row. |

The supplied manifest already shows a valid design: two control-condition ChIP
replicates, two treated-condition ChIP replicates, and matched inputs. Replace
the example names and paths with your own data.

For each factor, check the following before conversion:

- exactly two condition names are present;
- both condition names use exactly the same spelling on every row;
- each condition has at least ChIP replicate 1 and ChIP replicate 2;
- replicate numbers are not repeated within a factor and condition;
- every ChIP `control` value matches an input `sample_id` exactly;
- every ChIP row has a factor, tissue, and peak mode; and
- one factor does not mix `narrow` and `broad`.

Set the real reference paths. This human example assumes the reference helper
was used and that you added the blacklist:

```bash
ROOT="$(pwd)"
GENOME="$ROOT/references/human/hg38.fa"
CHROM_SIZES="$ROOT/references/human/hg38.chrom.sizes"
BT2_INDEX="$ROOT/references/human/hg38"
BLACKLIST="$ROOT/references/human/hg38-blacklist.v2.bed"
```

Generate the configuration:

```bash
bash make_config.sh \
  sample_manifest.project.tsv \
  --species human \
  --genome "$GENOME" \
  --chrom-sizes "$CHROM_SIZES" \
  --bt2-index "$BT2_INDEX" \
  --black-list "$BLACKLIST" \
  --effective-genome-size 2862010428 \
  --fastq-screen-conf config/fastq_screen.conf \
  --numerator-condition treated \
  --reference-condition control \
  --output config.yaml
```

This command replaces `config.yaml`. The original example remains available as
`config.sample.yaml`.

The effective size in this example is for uniquely mappable 150-bp GRCh38
reads. Choose a value for the actual assembly, read length, and mapping policy.
Change the species and paths for mouse or rat. Change `treated` and `control`
when your two condition labels are different.

The numerator and reference settings determine the direction of the DiffBind
contrast. For example, numerator `treated` and reference `control` means the
reported contrast is treated versus control.

### Option B: edit YAML directly

Start from the example:

```bash
cp config.sample.yaml config.yaml
```

Edit `config.yaml` carefully. YAML uses spaces for indentation. Do not use tabs.
The main sections are:

```yaml
species: human

references:
  human:
    genome: "/absolute/path/to/hg38.fa"
    chrom_sizes: "/absolute/path/to/hg38.chrom.sizes"
    black_list: "/absolute/path/to/hg38-blacklist.v2.bed"
    bt2_index: "/absolute/path/to/hg38"
    name: "hg38"
    gsize: "2.7e9"
    # Example for uniquely mappable 150-bp GRCh38 reads.
    effective_genome_size: 2862010428
```

Then set the screen, alignment, peaks, contrast, and plot options:

```yaml
contamination:
  fastq_screen_conf: "config/fastq_screen.conf"
  subset: 100000

alignment:
  min_mapq: 30
  max_insert_size: 2000

peak_calling:
  narrow_qvalue: 0.01
  broad_cutoff: 0.1
  consensus_min_replicates: 2

differential_binding:
  numerator_condition: treated
  reference_condition: control
  narrow_summits: 200

motif_enrichment:
  window_bp: 200
  background_multiplier: 2
  seed: 1

deeptools:
  threads: 4
  track_bin_size: 25
  matrix_bin_size: 50
  reference_point_upstream: 3000
  reference_point_downstream: 3000
  scale_regions_upstream: 3000
  scale_regions_downstream: 3000
  scale_regions_body_length: 5000
  log2_pseudocount: 1
  log2_scale_method: "None"
  heatmap_color_map: "RdBu_r"
  heatmap_z_min: -2
  heatmap_z_max: 2
  plot_dpi: 200
```

`motif_enrichment.window_bp` sets the fixed midpoint-centered sequence width.
`background_multiplier: 2` asks for two non-overlapping background windows per
foreground window, stratified by chromosome. `seed` makes those sampled windows
reproducible. Keep these values unchanged for a first run.

Finally, list input libraries under `chip_controls` and antibody-enriched
libraries under `chip_samples`:

```yaml
chip_controls:
  - id: Input_control_rep1
    condition: control
    fastq:
      - data/raw/Input_control_rep1_R1.fastq.gz
      - data/raw/Input_control_rep1_R2.fastq.gz

chip_samples:
  - id: H3K27ac_control_rep1
    fastq:
      - data/raw/H3K27ac_control_rep1_R1.fastq.gz
      - data/raw/H3K27ac_control_rep1_R2.fastq.gz
    control: Input_control_rep1
    condition: control
    replicate: 1
    factor: H3K27ac
    tissue: my_cell_type
    peak_mode: narrow
```

This single sample is only an illustration. A valid factor needs at least two
such ChIP entries for each of the two conditions.

### What the main settings mean

- `min_mapq: 30` retains confidently mapped proper pairs. Lowering it admits
  more ambiguous alignments and should have a biological reason.
- `max_insert_size: 2000` sets Bowtie2's maximum valid paired-end fragment
  length explicitly. Lower values can silently discard long ChIP fragments.
- `effective_genome_size` controls RPGC normalization. With a MAPQ filter, use
  a read-length-specific uniquely mappable value rather than a non-N genome
  length.
- `narrow_qvalue` is the MACS3 q-value cutoff for narrow peaks.
- `broad_cutoff` is the broad-region cutoff passed to MACS3.
- `consensus_min_replicates: 2` keeps intervals present in at least two peak
  files within each condition. It must not exceed the number of replicates in
  either condition.
- `narrow_summits: 200` tells DiffBind to recenter narrow intervals around
  summits with 200 bases on each side. Broad intervals retain their widths.
- `track_bin_size` controls BigWig resolution. Smaller bins create larger files
  and take longer to compute.
- `matrix_bin_size` controls the deepTools heatmap/profile bins. Window and body
  lengths must be exact multiples of this value.
- `log2_pseudocount` prevents division by zero in ChIP/input ratios.
- `log2_scale_method: "None"` independently CPM-normalizes ChIP and input before
  their log2 ratio. The accepted alternatives are `readCount` and `SES`.
- `heatmap_z_min` and `heatmap_z_max` clip the displayed color scale. They must
  lie below and above zero.

The defaults are a reasonable starting point, not a substitute for checking
the target, antibody, fragment size, sequencing depth, and experimental design.

## 12. Run the preflight checks

First confirm the installation:

```bash
bash setup.sh --check
```

Search for placeholders that should have been replaced:

```bash
grep -R "/path/to" config.yaml config/fastq_screen.conf
```

No output should be returned.

Confirm the main files exist:

```bash
test -s config.yaml
test -s config/fastq_screen.conf
test -s "$GENOME"
test -s "$GENOME.fai"
test -s "$CHROM_SIZES"
test -s "$BLACKLIST"
```

If you edited YAML directly and did not set the shell variables from the
manifest example, replace `$GENOME`, `$CHROM_SIZES`, and `$BLACKLIST` with your
real paths.

Now ask Snakemake to build the plan without running any analysis:

```bash
bash run.sh --dry-run
```

A dry-run reads the configuration, validates the design, checks declared input
paths, and prints the jobs that would run. It also rejects FastQ Screen files
that still contain `/path/to/` placeholders, but it cannot prove that every
configured database index is complete. It does not trim, align, or analyze
reads.

Read every error. The workflow tries to report all design problems together,
so one dry-run may list several items that need correction.

## 13. Run the analysis

When the dry-run succeeds, start the full workflow:

```bash
bash run.sh --cores 8
```

Choose a core count appropriate for the computer. Eight or twelve cores is a
reasonable start on a workstation. On a shared computer, follow local resource
rules. A high core count does not remove disk or memory limits.

The default run creates all main outputs, including QC, peaks, BigWigs,
deepTools plots, DiffBind results, narrow-factor motifs, provenance, and
MultiQC.

Snakemake prints each command and the job progress. Rule logs are written under:

```text
results/logs/
```

If the terminal closes or a job fails, correct the cause and run the same
command again:

```bash
bash run.sh --cores 8 --rerun-incomplete
```

Completed outputs are reused. Do not delete the `.snakemake/` directory while
a workflow is running.

To watch a particular log in another terminal:

```bash
tail -f results/logs/bowtie2_H3K27ac_control_rep1.log
```

Press `Ctrl-C` to stop watching the log. This does not stop the workflow in the
other terminal.

## 14. Run one result at a time

You can give `run.sh` a target path. Snakemake builds that file and everything
it needs.

```bash
# One alignment QC file
bash run.sh --cores 8 results/qc/samtools/H3K27ac_control_rep1.flagstat.txt

# One blacklist-filtered peak file
bash run.sh --cores 8 results/peaks/H3K27ac_control_rep1_peaks.narrowPeak

# One factor heatmap
bash run.sh --cores 8 results/deeptools/H3K27ac/heatmap.png

# One profile table
bash run.sh --cores 8 results/deeptools/H3K27ac/profile.tsv

# Differential binding
bash run.sh --cores 8 results/diffbind/H3K27ac/diffbind_summary.csv

# Motifs for a narrow factor
bash run.sh --cores 8 results/motifs/H3K27ac/motif_summary.pdf

# Full workflow plus the final combined report
bash run.sh --cores 8 results/multiqc/multiqc_report.html

# Optional Snakemake report
bash run.sh --cores 4 results/report/snakemake_report.html
```

Replace `H3K27ac` with the factor's output slug. Spaces and unsafe punctuation
in factor names are converted to hyphens. Sample IDs themselves may not contain
spaces.

Motif enrichment is only available for factors configured as `narrow`.

## 15. Results layout

After a complete run, `results/` has this general structure:

```text
results/
├── fastqc/
│   ├── raw/                         raw-read FastQC HTML and ZIP files
│   └── trimmed/                     post-trimming FastQC files
├── contamination/fastq_screen/      FastQ Screen HTML and text files
├── trimmed/                         paired trimmed FASTQs
├── bam/                             duplicate-removed BAMs and indexes
├── qc/samtools/                     flagstat and duplicate-removal metrics
├── peaks/
│   ├── raw/                         direct MACS3 outputs
│   ├── consensus/<factor>/          condition-level supported intervals
│   ├── consensus/<factor>.bed       union across the two conditions
│   └── <sample>_peaks.*Peak         blacklist-filtered sample peaks
├── bigwig/                          RPGC and log2 ChIP/input tracks
├── deeptools/<factor>/              matrices, heatmaps, profiles, and QC
├── diffbind/
│   ├── sample_sheet.csv             generated DiffBind design table
│   └── <factor>/                    contrast, table, plots, and R object
├── motifs/<factor>/                 full enrichment table and motif summary
├── provenance/resolved_config.yaml  configuration recorded for this run
├── multiqc/multiqc_report.html      combined QC report
├── logs/                            one log per workflow step
└── report/                          optional Snakemake HTML report
```

The best way to review a run is not to start with the final differential table.
Review the experiment in the order shown below.

## 16. Review raw and trimmed reads

Open:

```text
results/multiqc/multiqc_report.html
```

You can also open the individual FastQC reports under `results/fastqc/`.

Check:

- whether base quality falls sharply near the read ends;
- whether adapter sequence is present before trimming and reduced afterward;
- whether read lengths after trimming are still useful;
- whether one sample differs strongly from all others;
- whether duplication is unusually high; and
- whether the number of reads is comparable to the planned sequencing depth.

ChIP-seq can have genuine duplicated fragments when enrichment is strong, so a
FastQC duplication warning is not automatically a failed library. Compare it
with the samtools duplicate-removal count, peak enrichment, and other samples.

Open the FastQ Screen HTML files under:

```text
results/contamination/fastq_screen/
```

For a human experiment, many reads should map to the human database. A large
unexpected signal from another organism or PhiX needs investigation. The exact
expected pattern depends on the databases in your FastQ Screen configuration.

## 17. Review alignment and duplicate removal

Bowtie2 logs are written to:

```text
results/logs/bowtie2_<sample>.log
```

The pipeline keeps concordant proper pairs with the configured minimum mapping
quality and removes secondary, unmapped, mate-unmapped, QC-failed, and
duplicate-marked alignments as applicable. It then removes duplicate pairs with
samtools.

Review:

```text
results/qc/samtools/<sample>.flagstat.txt
results/qc/samtools/<sample>.markdup.txt
```

Look for:

- a reasonable number of mapped proper pairs;
- a severe loss of reads in one sample compared with its peers;
- very high duplicate removal; and
- large differences between ChIP replicates that were prepared together.

There is no universal mapping or duplication threshold that makes every
ChIP-seq experiment valid. Poor reference choice, contamination, low library
complexity, under-sequencing, and a failed immunoprecipitation can produce
different combinations of symptoms.

## 18. Review peaks

MACS3 calls each ChIP sample against its configured input. Direct outputs are
under:

```text
results/peaks/raw/
```

The workflow removes intervals overlapping the configured blacklist and sorts
the remaining intervals into:

```text
results/peaks/<sample>_peaks.narrowPeak
results/peaks/<sample>_peaks.broadPeak
```

Only one extension exists for a given sample, according to its `peak_mode`.

Count filtered intervals with:

```bash
wc -l results/peaks/H3K27ac_control_rep1_peaks.narrowPeak
```

Peak count alone is not a quality score. Very many weak peaks can indicate a
permissive call or noisy background, while few peaks can reflect a focused
factor, low sequencing depth, or a failed ChIP. Inspect the tracks in a genome
browser and compare biological replicates.

### Replicate-supported consensus peaks

Within each factor and condition, bedtools keeps intervals supported by at
least `consensus_min_replicates` sample peak files:

```text
results/peaks/consensus/H3K27ac/control.bed
results/peaks/consensus/H3K27ac/treated.bed
```

The two condition files are then merged into:

```text
results/peaks/consensus/H3K27ac.bed
```

This factor file contains regions reproducible within either condition. It is
not the intersection of both conditions. A treatment-specific region can be
included if it is supported by enough treated replicates.

This overlap procedure is not IDR. Do not describe the result as an IDR peak
set or as ENCODE certification.

## 19. Review signal tracks in a genome browser

The workflow makes two kinds of BigWig:

```text
results/bigwig/<sample>.rpgc.bw
results/bigwig/<chip-sample>.log2ratio.bw
```

RPGC tracks are produced for ChIP and input libraries. They are useful for
viewing library signal and for the factor correlation/PCA analysis.

The log2-ratio track compares each ChIP with its matched input. With the default
settings:

- zero means equal normalized ChIP and input signal in that bin;
- a positive value means relatively greater ChIP signal; and
- a negative value means relatively greater input signal.

Open the BigWigs together with the sample and consensus peak BED files in a
genome browser such as IGV. Check several known positive loci, quiet regions,
and regions that appear condition-specific. A convincing figure at one locus
does not replace whole-experiment QC.

## 20. Read the deepTools plots

All factor-level plots are under:

```text
results/deeptools/<factor>/
```

The matrix uses the factor consensus union from both conditions and one
matched-control log2 track for each ChIP sample. It therefore shows every
replicate over one common region set, including regions supported in only one
of the two conditions.

### Heatmap

Open `heatmap.png` for quick viewing and `heatmap.pdf` for a scalable figure.

For a narrow factor:

- each row is one factor consensus interval;
- the center column is the interval center;
- the default window covers 3 kb before and 3 kb after the center; and
- rows are ordered from higher to lower mean log2 ChIP/input signal.

For a broad factor:

- each row is one broad consensus interval;
- the interval body is scaled to the configured common length;
- the default includes 3 kb flanks on both sides; and
- scaling makes different domains comparable but hides their true width
  differences in the heatmap.

With the default `RdBu_r` scale, red represents positive log2 ChIP/input and
blue represents negative values. The default display limits are -2 and 2.
Values beyond those limits are shown with the end color; they are not changed
in the underlying matrix. Gray marks missing data.

The reusable data are:

```text
matrix.gz       deepTools binary matrix
matrix.tsv      numeric matrix
regions.bed     plotted regions in heatmap row order
```

Use `regions.bed` when you need to identify a particular heatmap row.

### Average profile

Open `profile.png` or `profile.pdf`. The profile uses the same matrix as the
heatmap and shows the mean signal across regions for every ChIP sample.

A sharp central signal is common for a successful narrow-factor analysis. A
broad-mark profile should be read across the scaled body and its flanks.

An average can hide important subgroups. A smooth profile does not prove that
every region behaves the same way. Always inspect the heatmap and the profile
together. Numeric profile values are in `profile.tsv`.

### Fingerprint

Open:

```text
results/deeptools/<factor>/qc/fingerprint.png
```

The plot compares read concentration in the ChIP libraries and their matched
inputs. An enriched ChIP often bends farther from the diagonal because a
larger fraction of reads is concentrated in a smaller fraction of genomic
bins. Interpret the curve together with the table in
`fingerprint_metrics.tsv`, not as a pass/fail test by itself.

### Correlation

Open:

```text
results/deeptools/<factor>/qc/spearman_heatmap.png
```

The values summarize each ChIP sample's RPGC signal over the factor consensus
peaks; inputs are not included in this matrix. Biological replicates usually
correlate well, but there is no universal acceptable number. A strong treatment
may also make conditions differ. An isolated replicate should be checked
against FastQC, alignment, duplicate, peak, and metadata evidence before it is
excluded.

The numeric matrix is `spearman.tsv`.

### PCA

Open:

```text
results/deeptools/<factor>/qc/pca.png
```

Each point is a ChIP sample; input controls are not included. Points that are
close have similar peak-level RPGC signal. Replicates should usually be closer
to each other than to unrelated samples, and condition separation can support
a biological effect. PCA does not prove an effect and does not diagnose its
cause. This PCA is calculated over a peak set derived from the same experiment,
so confirm the pattern with other QC and biology.

Coordinates and variance information are in `pca.tsv`.

## 21. Read the DiffBind results

Before interpreting differential binding, open:

```text
results/diffbind/sample_sheet.csv
results/diffbind/<factor>/contrast.tsv
```

Confirm sample names, BAMs, controls, peak files, conditions, replicates, peak
mode, numerator, and reference.

The main outputs are:

```text
results/diffbind/<factor>/diffbind_summary.csv
results/diffbind/<factor>/diffbind_plots.pdf
results/diffbind/<factor>/diffbind.rds
results/diffbind/<factor>/contrast.tsv
```

DiffBind counts reads over a common set of intervals supported by at least two
sample peak sets, creates the configured numerator-versus-reference contrast,
and reports sites passing an FDR threshold of 0.05.

In a treated-versus-control contrast, a positive `Fold` value indicates greater
binding in treated and a negative value indicates greater binding in control.
This follows the numerator-minus-reference direction recorded in
`contrast.tsv`. Confirm that file before labeling any site as gained or lost.

The CSV can contain no rows. That is a valid result when no site passes the
0.05 reporting threshold; it is not necessarily a software error.

The PDF contains overview, PCA, heatmap, and MA plots. Use these to check sample
behavior and effect-size patterns. The saved `.rds` object is for further work
in R by someone familiar with DiffBind.

The supplied model does not account for batches, paired donors, or other
covariates. Do not use it to make adjusted claims about a design that contains
those effects.

## 22. Read the motif results

For a narrow factor, the workflow centers a 200 bp sequence window on each
factor-consensus interval midpoint and tests JASPAR vertebrate CORE motifs with
`monaLisa`. These midpoints are not retained MACS summit coordinates.

That consensus BED is the union of both condition-level consensus files. The
default motif result therefore describes the combined factor region set, not a
treated-only or control-only motif comparison.

Open:

```text
results/motifs/<factor>/motif_summary.pdf
results/motifs/<factor>/motif_summary.csv
```

The PDF shows up to 20 motifs ranked by adjusted significance. The complete
table, including effect size and raw and adjusted significance, is:

```text
results/motifs/<factor>/motif_enrichment.tsv
```

The background sequences are sampled reproducibly in the same chromosome
proportions as the foreground. They are non-overlapping and exclude both the
blacklist and the complete factor-consensus intervals. `monaLisa` then weights
them to reduce GC and short k-mer composition differences. No precomputed
background BED is required, but the configured assembly-matched FASTA, index,
and blacklist are required. This background does not automatically match assay
accessibility or mappability.

An enriched motif is an association. It does not prove that its named
transcription factor bound those regions. Related factors can share motifs,
sequence composition affects enrichment, and the background choice matters.
Compare the result with the ChIP target, known partners, expression evidence,
and independent experiments.

This is enrichment in the combined factor region set, not differential motif
analysis between conditions. The default workflow does not run fixed-window
motif analysis for broad domains, where the biological interpretation is
usually unclear.

## 23. Review the combined reports and provenance

The MultiQC report is part of the default run:

```text
results/multiqc/multiqc_report.html
```

It is the quickest entry point for read and alignment QC, but it does not
replace the factor-specific deepTools, DiffBind, motif, and browser review.

The exact resolved configuration used by Snakemake is saved at:

```text
results/provenance/resolved_config.yaml
```

Keep this file with the results. Also record the repository commit, sample
metadata, reference source and assembly, antibody information, and sequencing
run details outside the pipeline.

The optional Snakemake report is built only when requested:

```bash
bash run.sh --cores 4 results/report/snakemake_report.html
```

After the full run, repeat the dry-run:

```bash
bash run.sh --dry-run
```

If the inputs and configuration have not changed, Snakemake should report that
nothing needs to be done. Also confirm the final report and recorded config are
not empty:

```bash
test -s results/multiqc/multiqc_report.html && echo "MultiQC report is present"
test -s results/provenance/resolved_config.yaml && echo "Run config is present"
```

## 24. Relative normalization and global changes

The RPGC, CPM, log2 ChIP/input, and default DiffBind results are relative to
the sequenced libraries. They assume that a useful reference distribution
remains across samples.

They cannot establish a genome-wide gain or loss of chromatin occupancy when
most of the genome changes in the same direction. A sample with a true global
increase can still be rescaled to look similar overall.

If the biological question concerns a global shift, use an appropriate
experimental spike-in and a validated spike-in-aware analysis. Adding spike-in
reads after sequencing or changing only a plotting scale does not solve this
problem. This workflow does not implement spike-in normalization.

Phrase conclusions as relative enrichment or relative differential binding
unless the experimental design supports an absolute or global claim.

## 25. What this workflow does not prove

A completed run does not by itself prove that:

- the antibody was specific;
- the ChIP experiment succeeded;
- every peak is a direct binding event;
- consensus peaks passed IDR;
- the experiment meets all ENCODE standards;
- a motif identifies the protein that bound the DNA;
- a PCA separation is biological rather than technical;
- a differential site causes a change in gene expression; or
- library-normalized signal measures a global gain or loss.

Good analysis combines the pipeline outputs with experimental metadata,
antibody validation, known positive and negative loci, genome-browser review,
independent assays, and an appropriate statistical design.

## 26. Troubleshooting

### `Conda is missing` or the runner is missing

Run:

```bash
bash setup.sh
bash setup.sh --check
```

If you set `MINIFORGE_HOME` during installation, export the same value before
running the workflow.

### Setup was interrupted

Run `bash setup.sh` again. Completed downloads and complete environments are
reused. If Miniforge installation itself stopped and its prefix exists without
an executable `bin/conda`, inspect that path carefully, remove or move only the
incomplete Miniforge directory, and rerun setup. The default prefix is
`$HOME/miniforge3`; a custom prefix comes from `MINIFORGE_HOME`.

### The dry-run reports missing FASTQs

Check each path exactly as it appears in `config.yaml`:

```bash
ls -lh /the/path/from/config.fastq.gz
```

Relative paths are interpreted from the repository root. Check capitalization,
file extensions, and R1/R2 order.

### The Bowtie2 index is missing or incomplete

`bt2_index` must be a prefix, not a `.1.bt2` filename and not just a directory.
All six small-index files or all six large-index files must exist. Rebuild an
incomplete index.

### Chromosome names do not match

The FASTA, chromosome-size file, blacklist, BAMs, and BED files must use the
same names. `chr1` and `1` are not interchangeable. Compare the first column
of the chromosome-size and blacklist files. Obtain or convert the correct
assembly-matched resource before running the analysis. Do not simply add or
remove `chr` without confirming every contig.

### FastQ Screen fails immediately

Check for remaining placeholders:

```bash
grep "/path/to" config/fastq_screen.conf
```

Then confirm every `DATABASE` value is a real Bowtie2 prefix. Remove unused
database lines rather than leaving broken paths.

### YAML cannot be parsed

Use spaces, not tabs. Keep list dashes aligned. Quote paths containing special
characters. Starting again from `config.sample.yaml` is often faster than
repairing badly indented YAML.

### The configuration reports an invalid design

Read the complete list in the error. Common causes are:

- only one condition;
- only one biological ChIP replicate in a condition;
- condition spelling that does not match the numerator/reference settings;
- a repeated replicate number;
- an unknown input control ID;
- mixed narrow and broad modes for one factor; or
- more than one tissue label for one factor.

Correct the metadata rather than renaming samples just to bypass the check.

### No consensus peaks remain

The workflow stops when a factor/condition has no intervals supported by the
required number of replicate peak files. Inspect the individual filtered peak
files, MACS3 logs, alignment QC, controls, and browser tracks.

Do not lower `consensus_min_replicates` below two. If two genuine biological
replicates have almost no peak overlap, the disagreement is a result that needs
investigation.

### The heatmap is blank or nearly flat

Check:

- whether the factor consensus BED contains intervals;
- whether log2-ratio BigWigs contain signal in a genome browser;
- whether chromosome names match across BAMs, BigWigs, BEDs, and references;
- whether ChIP and input were assigned correctly;
- whether the chosen narrow/broad mode is suitable; and
- whether the fixed display limits hide a much smaller signal range.

Inspect `matrix.tsv` before changing colors or limits.

### DiffBind returns no significant sites

First check that the run completed and the CSV exists. An empty table can mean
that no interval passed FDR 0.05. Review the contrast, replicate QC, sample
sheet, effect sizes, and sequencing depth. Do not raise the FDR threshold only
to obtain a desired result.

### A workflow lock remains after a crash

Make sure no Snakemake process is still running. Only then unlock the working
directory:

```bash
bash run.sh --cores 1 --unlock
```

Never unlock a directory while another workflow is using it.

### A rule failed

Snakemake prints the rule name and log path. Read the end of that log:

```bash
tail -n 100 results/logs/the_failed_rule.log
```

Fix the input, reference, disk-space, memory, or configuration problem, then
rerun with `--rerun-incomplete`.

### The disk is full

Check free space:

```bash
df -h .
```

Do not delete files while Snakemake is running. After the run has stopped,
move unrelated data or extend the filesystem, then rerun. BAMs, trimmed FASTQs,
BigWigs, Conda packages, and `.snakemake/conda/` can all use substantial space.

## 27. Final checklist

Before the full run:

- [ ] Linux x86_64 or aarch64 is available.
- [ ] `bash setup.sh --check` passes.
- [ ] Every library has distinct R1 and R2 gzip-compressed FASTQs, and no FASTQ is assigned twice.
- [ ] Gzip integrity checks pass.
- [ ] Genome, FASTA index, chromosome sizes, blacklist, and Bowtie2 index use one assembly.
- [ ] The Bowtie2 prefix has all six index files.
- [ ] FastQ Screen contains no placeholder database paths.
- [ ] Every ChIP has a biologically matched input ID with the same condition label.
- [ ] Each factor has exactly two conditions.
- [ ] Each condition has at least two biological ChIP replicates.
- [ ] Replicate numbers are unique within each factor and condition.
- [ ] Each factor uses one tissue label and one peak mode.
- [ ] DiffBind numerator and reference labels match the sample conditions.
- [ ] `effective_genome_size` is a positive, unquoted integer chosen for the assembly, read length, and mapping policy.
- [ ] `bash run.sh --dry-run` succeeds.

After the run:

- [ ] Review MultiQC and individual FastQC/FastQ Screen reports.
- [ ] Review Bowtie2, flagstat, and duplicate-removal metrics.
- [ ] Compare replicate peak calls and condition consensus peaks.
- [ ] Inspect RPGC and log2 ChIP/input tracks in a genome browser.
- [ ] Read the heatmap and profile together.
- [ ] Check fingerprint, correlation, and PCA for unexpected samples.
- [ ] Confirm DiffBind contrast direction before naming gains or losses.
- [ ] Treat motifs as enrichment evidence, not direct binding proof.
- [ ] Keep `resolved_config.yaml` and experimental metadata with the results.
- [ ] State the limits of relative normalization in the conclusions.

## 28. Short glossary

- **BAM:** A compressed file of reads aligned to a reference genome.
- **BED:** A tab-separated file describing genomic intervals.
- **BigWig:** An indexed signal track designed for fast genome-browser display.
- **Blacklist:** A set of assembly-specific regions known to produce
  artifactual signal in many sequencing experiments.
- **Consensus peak:** In this workflow, an interval supported by a configured
  minimum number of replicate peak files within one condition.
- **Control or input:** Chromatin sequenced without the target-specific
  immunoprecipitation and used to model background.
- **FDR:** False discovery rate, an estimate of the proportion of reported
  discoveries expected to be false under the statistical model.
- **MAPQ:** Mapping quality, a score describing confidence in an alignment
  location.
- **Peak:** A genomic interval with ChIP signal enriched relative to background.
- **RPGC:** Reads per genomic content, a deepTools coverage normalization that
  uses an effective genome size.
- **Technical replicate:** Repeated measurement or sequencing of the same
  biological material.
- **Biological replicate:** An independently prepared biological sample
  representing the same condition.

## 29. Further reading

- [Analysis and visualization notes](docs/ANALYSIS_AND_VISUALIZATION.md)
- [deepTools 3.5.5 documentation](https://deeptools.readthedocs.io/en/3.5.5/)
- [MACS3 callpeak documentation](https://macs3-project.github.io/MACS/docs/callpeak.html)
- [DiffBind reference manual](https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf)
- [ENCODE histone ChIP-seq standards](https://www.encodeproject.org/chip-seq/histone-encode4/)
- [ENCODE transcription-factor ChIP-seq standards](https://www.encodeproject.org/chip-seq/transcription-factor-encode4/)
