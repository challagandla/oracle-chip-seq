def config_int(value, path):
    """Read an actual YAML integer and report a configuration error clearly."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise WorkflowError(f"{path} must be a YAML integer, got {value!r}")
    return value


def config_float(value, path):
    """Read a finite YAML number without accepting booleans or numeric text."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkflowError(f"{path} must be a YAML number, got {value!r}")
    value = float(value)
    if not value == value or value in {float("inf"), float("-inf")}:
        raise WorkflowError(f"{path} must be finite, got {value!r}")
    return value


def config_mapping(value, path):
    if not isinstance(value, dict):
        raise WorkflowError(f"{path} must be a YAML mapping, got {type(value).__name__}")
    return value


SPECIES = config.get("species")
REFERENCES = config_mapping(config.get("references", {}), "references")
REF = config_mapping(REFERENCES.get(SPECIES, {}), f"references.{SPECIES}")
CONTAMINATION = config_mapping(config.get("contamination", {}), "contamination")
FASTQ_SCREEN_CONF = CONTAMINATION.get("fastq_screen_conf", "config/fastq_screen.conf.example")
FASTQ_SCREEN_SUBSET = config_int(CONTAMINATION.get("subset", 100000), "contamination.subset")

ALIGNMENT = config_mapping(config.get("alignment", {}), "alignment")
MIN_MAPQ = config_int(ALIGNMENT.get("min_mapq", 30), "alignment.min_mapq")
MAX_INSERT_SIZE = config_int(
    ALIGNMENT.get("max_insert_size", 2000), "alignment.max_insert_size"
)

PEAK_CALLING = config_mapping(config.get("peak_calling", {}), "peak_calling")
NARROW_QVALUE = config_float(PEAK_CALLING.get("narrow_qvalue", 0.01), "peak_calling.narrow_qvalue")
BROAD_CUTOFF = config_float(PEAK_CALLING.get("broad_cutoff", 0.1), "peak_calling.broad_cutoff")
CONSENSUS_MIN_REPLICATES = config_int(
    PEAK_CALLING.get("consensus_min_replicates", 2),
    "peak_calling.consensus_min_replicates",
)

DIFFERENTIAL_BINDING = config_mapping(
    config.get("differential_binding", {}), "differential_binding"
)
DIFFBIND_NUMERATOR_RAW = DIFFERENTIAL_BINDING.get("numerator_condition", "treated")
DIFFBIND_REFERENCE_RAW = DIFFERENTIAL_BINDING.get("reference_condition", "control")
DIFFBIND_NUMERATOR = str(DIFFBIND_NUMERATOR_RAW).strip()
DIFFBIND_REFERENCE = str(DIFFBIND_REFERENCE_RAW).strip()
DIFFBIND_NARROW_SUMMITS = config_int(
    DIFFERENTIAL_BINDING.get("narrow_summits", 200),
    "differential_binding.narrow_summits",
)

DEEPTOOLS = config_mapping(config.get("deeptools", {}), "deeptools")
DEEPTOOLS_THREADS = config_int(DEEPTOOLS.get("threads", 4), "deeptools.threads")
TRACK_BIN_SIZE = config_int(DEEPTOOLS.get("track_bin_size", 25), "deeptools.track_bin_size")
MATRIX_BIN_SIZE = config_int(DEEPTOOLS.get("matrix_bin_size", 50), "deeptools.matrix_bin_size")
REFERENCE_POINT_UPSTREAM = config_int(
    DEEPTOOLS.get("reference_point_upstream", 3000),
    "deeptools.reference_point_upstream",
)
REFERENCE_POINT_DOWNSTREAM = config_int(
    DEEPTOOLS.get("reference_point_downstream", 3000),
    "deeptools.reference_point_downstream",
)
SCALE_REGIONS_UPSTREAM = config_int(
    DEEPTOOLS.get("scale_regions_upstream", 3000),
    "deeptools.scale_regions_upstream",
)
SCALE_REGIONS_DOWNSTREAM = config_int(
    DEEPTOOLS.get("scale_regions_downstream", 3000),
    "deeptools.scale_regions_downstream",
)
SCALE_REGIONS_BODY_LENGTH = config_int(
    DEEPTOOLS.get("scale_regions_body_length", 5000),
    "deeptools.scale_regions_body_length",
)
LOG2_PSEUDOCOUNT = config_float(
    DEEPTOOLS.get("log2_pseudocount", 1), "deeptools.log2_pseudocount"
)
LOG2_SCALE_METHOD_RAW = DEEPTOOLS.get("log2_scale_method", "None")
LOG2_SCALE_METHOD = str(LOG2_SCALE_METHOD_RAW)
HEATMAP_ZMIN = config_float(DEEPTOOLS.get("heatmap_z_min", -2), "deeptools.heatmap_z_min")
HEATMAP_ZMAX = config_float(DEEPTOOLS.get("heatmap_z_max", 2), "deeptools.heatmap_z_max")
HEATMAP_COLOR_MAP_RAW = DEEPTOOLS.get("heatmap_color_map", "RdBu_r")
HEATMAP_COLOR_MAP = str(HEATMAP_COLOR_MAP_RAW)
PLOT_DPI = config_int(DEEPTOOLS.get("plot_dpi", 200), "deeptools.plot_dpi")

MOTIF_ENRICHMENT = config_mapping(
    config.get("motif_enrichment", {}), "motif_enrichment"
)
MOTIF_WINDOW_BP = config_int(
    MOTIF_ENRICHMENT.get("window_bp", 200), "motif_enrichment.window_bp"
)
MOTIF_BACKGROUND_MULTIPLIER = config_int(
    MOTIF_ENRICHMENT.get("background_multiplier", 2),
    "motif_enrichment.background_multiplier",
)
MOTIF_SEED = config_int(
    MOTIF_ENRICHMENT.get("seed", 1), "motif_enrichment.seed"
)
if MOTIF_WINDOW_BP <= 0:
    raise WorkflowError("motif_enrichment.window_bp must be a positive integer")
if MOTIF_BACKGROUND_MULTIPLIER <= 0:
    raise WorkflowError(
        "motif_enrichment.background_multiplier must be a positive integer"
    )
if MOTIF_SEED < 0:
    raise WorkflowError("motif_enrichment.seed must be a non-negative integer")


def sample_factor(sample):
    return str(sample.get("factor") or sample.get("mark") or "").strip()


def factor_slug(value):
    slug = re.sub(r"[^A-Za-z0-9._-]+", "-", str(value).strip()).strip("-._")
    return slug or "factor"


def peak_mode(sample):
    """Return the explicit biological peak geometry for this sample."""
    return str(sample.get("peak_mode") or "").strip().lower()


def validate_config(cfg):
    """Fail before DAG construction when a run would be ambiguous or invalid."""
    problems = []
    if SPECIES not in cfg.get("references", {}):
        problems.append(f"references.{SPECIES} is missing")
    for key in (
        "genome", "chrom_sizes", "black_list", "bt2_index", "gsize",
        "effective_genome_size",
    ):
        if key not in REF or REF[key] in (None, ""):
            problems.append(f"references.{SPECIES}.{key} is missing")
    for key in ("genome", "chrom_sizes", "black_list", "bt2_index"):
        if key in REF and (not isinstance(REF[key], str) or not REF[key].strip()):
            problems.append(f"references.{SPECIES}.{key} must be a non-empty path")
    effective_size = REF.get("effective_genome_size")
    if isinstance(effective_size, bool) or not isinstance(effective_size, int) or effective_size <= 0:
        problems.append(
            f"references.{SPECIES}.effective_genome_size must be a positive integer; "
            "scientific-notation text such as '2.7e9' is invalid for deepTools"
        )
    gsize = REF.get("gsize")
    gsize_text = str(gsize).strip()
    numeric_gsize = re.fullmatch(
        r"(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?",
        gsize_text,
    )
    if (
        isinstance(gsize, bool)
        or not isinstance(gsize, (str, int, float))
        or (gsize_text not in {"hs", "mm", "ce", "dm"} and numeric_gsize is None)
        or (
            numeric_gsize is not None
            and (
                float(gsize_text) <= 0
                or float(gsize_text) == float("inf")
            )
        )
    ):
        problems.append(
            f"references.{SPECIES}.gsize must be a positive MACS genome size or "
            f"one of hs/mm/ce/dm; got {gsize!r}"
        )

    controls = cfg.get("chip_controls", [])
    chip_samples = cfg.get("chip_samples", [])
    if not isinstance(controls, list):
        problems.append("chip_controls must be a YAML list")
        controls = []
    if not isinstance(chip_samples, list):
        problems.append("chip_samples must be a YAML list")
        chip_samples = []
    invalid_records = [
        f"{collection}[{index}]"
        for collection, records in (("chip_controls", controls), ("chip_samples", chip_samples))
        for index, record in enumerate(records)
        if not isinstance(record, dict)
    ]
    if invalid_records:
        problems.append("sample records must be YAML mappings: " + ", ".join(invalid_records))
        controls = [record for record in controls if isinstance(record, dict)]
        chip_samples = [record for record in chip_samples if isinstance(record, dict)]
    if not controls:
        problems.append("chip_controls is empty")
    if not chip_samples:
        problems.append("chip_samples is empty")

    all_records = controls + chip_samples
    ids = [str(item.get("id", "")).strip() for item in all_records]
    duplicate_ids = sorted(sample_id for sample_id, count in Counter(ids).items() if count > 1)
    if duplicate_ids:
        problems.append(f"duplicate sample IDs: {', '.join(duplicate_ids)}")
    safe_id = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
    fastq_owners = defaultdict(list)
    for record in all_records:
        raw_id = record.get("id", "")
        sid = str(raw_id).strip()
        if not isinstance(raw_id, str):
            problems.append(f"sample ID {raw_id!r} must be text")
        elif raw_id != sid:
            problems.append(
                f"sample ID {raw_id!r} has leading or trailing whitespace; remove it"
            )
        if not safe_id.fullmatch(sid):
            problems.append(
                f"sample ID '{sid}' is unsafe; use only letters, numbers, '.', '_' and '-'"
            )
        fastq = record.get("fastq")
        if (
            not isinstance(fastq, list)
            or len(fastq) != 2
            or not all(isinstance(path, str) and path.strip() for path in fastq)
        ):
            problems.append(
                f"{sid or '<unnamed>'}: this workflow requires exactly two paired-end FASTQ paths"
            )
        elif not all(path.lower().endswith((".fastq.gz", ".fq.gz")) for path in fastq):
            problems.append(
                f"{sid}: FASTQs must be gzip-compressed .fastq.gz or .fq.gz files"
            )
        else:
            normalized = [str(Path(path).resolve(strict=False)) for path in fastq]
            if normalized[0] == normalized[1]:
                problems.append(f"{sid}: R1 and R2 resolve to the same FASTQ file")
            for read, path in zip(("R1", "R2"), normalized):
                fastq_owners[path].append(f"{sid} {read}")

    for path, owners in sorted(fastq_owners.items()):
        if len(owners) > 1:
            problems.append(
                f"FASTQ file is assigned more than once ({path}): {', '.join(owners)}"
            )

    control_ids = {control.get("id") for control in controls}
    control_conditions = {}
    for control in controls:
        raw_condition = control.get("condition")
        condition = raw_condition.strip() if isinstance(raw_condition, str) else ""
        if not condition:
            problems.append(f"{control.get('id', '<unnamed>')}: control condition is missing")
        else:
            control_conditions[control.get("id")] = condition
    factor_condition_counts = Counter()
    replicates_by_factor_condition = defaultdict(list)
    modes_by_factor = defaultdict(set)
    tissues_by_factor = defaultdict(set)
    slugs = defaultdict(set)
    condition_slugs_by_factor = defaultdict(lambda: defaultdict(set))
    for sample in chip_samples:
        sid = sample.get("id", "")
        control = sample.get("control")
        if control not in control_ids:
            problems.append(f"{sid}: unknown or missing control '{control or ''}'")
        raw_factor = sample.get("factor") or sample.get("mark")
        raw_condition = sample.get("condition", "")
        raw_tissue = sample.get("tissue")
        factor = sample_factor(sample)
        condition = str(raw_condition).strip()
        tissue = raw_tissue.strip() if isinstance(raw_tissue, str) else ""
        mode = peak_mode(sample)
        if not isinstance(raw_factor, str) or not factor:
            problems.append(f"{sid}: factor is missing")
        if not isinstance(raw_condition, str) or not condition:
            problems.append(f"{sid}: condition is missing")
        if control in control_conditions and control_conditions[control] != condition:
            problems.append(
                f"{sid}: condition '{condition}' does not match control {control} "
                f"condition '{control_conditions[control]}'"
            )
        if not tissue:
            problems.append(f"{sid}: tissue is missing")
        replicate = None
        raw_replicate = sample.get("replicate", 0)
        if isinstance(raw_replicate, bool) or not isinstance(raw_replicate, int) or raw_replicate < 1:
            problems.append(f"{sid}: replicate must be a positive integer")
        else:
            replicate = raw_replicate
        if mode not in {"narrow", "broad"}:
            problems.append(f"{sid}: peak_mode must be 'narrow' or 'broad', got '{mode}'")
        factor_condition_counts[(factor, condition)] += 1
        if replicate is not None:
            replicates_by_factor_condition[(factor, condition)].append(replicate)
        modes_by_factor[factor].add(mode)
        tissues_by_factor[factor].add(tissue)
        slugs[factor_slug(factor)].add(factor)
        condition_slugs_by_factor[factor][factor_slug(condition)].add(condition)

    for factor, modes in modes_by_factor.items():
        if len(modes) > 1:
            problems.append(f"{factor}: samples mix peak modes: {', '.join(sorted(modes))}")
        if len(tissues_by_factor[factor]) > 1:
            problems.append(
                f"{factor}: the supplied simple DiffBind model requires one tissue per "
                f"factor; found {', '.join(sorted(tissues_by_factor[factor]))}"
            )
        conditions = sorted(
            condition for f, condition in factor_condition_counts if f == factor
        )
        if len(conditions) != 2:
            problems.append(
                f"{factor}: this workflow requires exactly two conditions for one "
                f"unambiguous DiffBind contrast; found {len(conditions)}"
            )
        if set(conditions) != {DIFFBIND_REFERENCE, DIFFBIND_NUMERATOR}:
            problems.append(
                f"{factor}: conditions must match differential_binding.reference_condition="
                f"'{DIFFBIND_REFERENCE}' and numerator_condition='{DIFFBIND_NUMERATOR}'; "
                f"found {', '.join(conditions)}"
            )
        for condition in conditions:
            count = factor_condition_counts[(factor, condition)]
            if count < 2:
                problems.append(
                    f"{factor}/{condition}: differential binding needs at least two replicates; found {count}"
                )
            if CONSENSUS_MIN_REPLICATES > count:
                problems.append(
                    f"{factor}/{condition}: consensus_min_replicates="
                    f"{CONSENSUS_MIN_REPLICATES} exceeds the {count} configured replicates"
                )
            replicates = replicates_by_factor_condition[(factor, condition)]
            if len(set(replicates)) != len(replicates):
                problems.append(
                    f"{factor}/{condition}: replicate numbers must be unique; got {replicates}"
                )
        for condition_slug, condition_names in condition_slugs_by_factor[factor].items():
            if len(condition_names) > 1:
                problems.append(
                    f"{factor}: condition names collapse to output slug '{condition_slug}': "
                    + ", ".join(sorted(condition_names))
                )
    for slug, factor_names in slugs.items():
        if len(factor_names) > 1:
            problems.append(
                f"factor names collapse to the same output slug '{slug}': "
                + ", ".join(sorted(factor_names))
            )

    if not 0 <= MIN_MAPQ <= 255:
        problems.append("alignment.min_mapq must be between 0 and 255")
    if MAX_INSERT_SIZE <= 0:
        problems.append("alignment.max_insert_size must be positive")
    if "require_proper_pairs" in ALIGNMENT and ALIGNMENT["require_proper_pairs"] is not True:
        problems.append(
            "alignment.require_proper_pairs may only be true; concordant proper pairs are "
            "an invariant of this paired-end workflow"
        )
    if FASTQ_SCREEN_SUBSET <= 0:
        problems.append("contamination.subset must be positive")
    if not isinstance(FASTQ_SCREEN_CONF, str) or not FASTQ_SCREEN_CONF.strip():
        problems.append("contamination.fastq_screen_conf must be a non-empty path")
    elif FASTQ_SCREEN_CONF.endswith(".example"):
        problems.append(
            "contamination.fastq_screen_conf points to a template; copy it to "
            "config/fastq_screen.conf and replace every database path"
        )
    elif Path(FASTQ_SCREEN_CONF).is_file():
        try:
            screen_config = Path(FASTQ_SCREEN_CONF).read_text(encoding="utf-8")
        except OSError as error:
            problems.append(f"cannot read contamination.fastq_screen_conf: {error}")
        else:
            if "/path/to/" in screen_config:
                problems.append(
                    "contamination.fastq_screen_conf still contains /path/to/ placeholders"
                )
    if not 0 < NARROW_QVALUE <= 1:
        problems.append("peak_calling.narrow_qvalue must be greater than 0 and at most 1")
    if not 0 < BROAD_CUTOFF <= 1:
        problems.append("peak_calling.broad_cutoff must be greater than 0 and at most 1")
    if CONSENSUS_MIN_REPLICATES < 2:
        problems.append("peak_calling.consensus_min_replicates must be at least 2")
    if not isinstance(DIFFBIND_REFERENCE_RAW, str) or not isinstance(DIFFBIND_NUMERATOR_RAW, str):
        problems.append("differential_binding condition names must be text")
    if not DIFFBIND_REFERENCE or not DIFFBIND_NUMERATOR:
        problems.append("differential_binding condition names must be non-empty")
    if DIFFBIND_REFERENCE == DIFFBIND_NUMERATOR:
        problems.append("differential_binding numerator and reference conditions must differ")
    if DIFFBIND_NARROW_SUMMITS <= 0:
        problems.append("differential_binding.narrow_summits must be positive")
    for name, value in {
        "peak_calling.consensus_min_replicates": CONSENSUS_MIN_REPLICATES,
        "deeptools.threads": DEEPTOOLS_THREADS,
        "deeptools.track_bin_size": TRACK_BIN_SIZE,
        "deeptools.matrix_bin_size": MATRIX_BIN_SIZE,
        "deeptools.reference_point_upstream": REFERENCE_POINT_UPSTREAM,
        "deeptools.reference_point_downstream": REFERENCE_POINT_DOWNSTREAM,
        "deeptools.scale_regions_upstream": SCALE_REGIONS_UPSTREAM,
        "deeptools.scale_regions_downstream": SCALE_REGIONS_DOWNSTREAM,
        "deeptools.scale_regions_body_length": SCALE_REGIONS_BODY_LENGTH,
        "deeptools.plot_dpi": PLOT_DPI,
    }.items():
        if value <= 0:
            problems.append(f"{name} must be positive")
    for name, value in {
        "deeptools.reference_point_upstream": REFERENCE_POINT_UPSTREAM,
        "deeptools.reference_point_downstream": REFERENCE_POINT_DOWNSTREAM,
        "deeptools.scale_regions_upstream": SCALE_REGIONS_UPSTREAM,
        "deeptools.scale_regions_downstream": SCALE_REGIONS_DOWNSTREAM,
        "deeptools.scale_regions_body_length": SCALE_REGIONS_BODY_LENGTH,
    }.items():
        if value > 0 and MATRIX_BIN_SIZE > 0 and value % MATRIX_BIN_SIZE != 0:
            problems.append(
                f"{name}={value} must be an exact multiple of "
                f"deeptools.matrix_bin_size={MATRIX_BIN_SIZE}"
            )
    if not isinstance(LOG2_SCALE_METHOD_RAW, str):
        problems.append("deeptools.log2_scale_method must be text")
    if LOG2_SCALE_METHOD not in {"readCount", "SES", "None"}:
        problems.append("deeptools.log2_scale_method must be readCount, SES, or None")
    factor_modes = DEEPTOOLS.get("factor_modes", {})
    if not isinstance(factor_modes, dict):
        problems.append("deeptools.factor_modes must be a mapping from factor to matrix mode")
        factor_modes = {}
    for factor, mode in factor_modes.items():
        if factor not in modes_by_factor:
            problems.append(f"deeptools.factor_modes contains unknown factor '{factor}'")
        if not isinstance(mode, str) or mode not in {
            "narrow", "broad", "reference_point", "scale_regions"
        }:
            problems.append(
                f"deeptools.factor_modes.{factor} must be narrow, broad, "
                "reference_point, or scale_regions"
            )
    if LOG2_PSEUDOCOUNT <= 0:
        problems.append("deeptools.log2_pseudocount must be positive")
    if not HEATMAP_ZMIN < 0 < HEATMAP_ZMAX:
        problems.append("deeptools heatmap limits must straddle zero")
    if not isinstance(HEATMAP_COLOR_MAP_RAW, str):
        problems.append("deeptools.heatmap_color_map must be text")
    if not re.fullmatch(r"[A-Za-z0-9_]+", HEATMAP_COLOR_MAP):
        problems.append("deeptools.heatmap_color_map must be a safe Matplotlib colormap name")
    if problems:
        raise WorkflowError("Invalid ChIP-seq configuration:\n  - " + "\n  - ".join(problems))


validate_config(config)

CHIP_SAMPLES = [s["id"] for s in config["chip_samples"]]
CHIP_CONTROLS = {c["id"]: c for c in config["chip_controls"]}
ALL_SAMPLES = CHIP_SAMPLES + list(CHIP_CONTROLS.keys())

FASTQ = {s["id"]: s["fastq"] for s in config["chip_samples"] + config["chip_controls"]}
CONTROL_MAP = {s["id"]: s["control"] for s in config["chip_samples"]}

# Peak geometry is explicit: histone/TF identity is not a substitute for
# broad/narrow biology (for example, H3K27ac is usually punctate).
PEAK_MODE = {s["id"]: peak_mode(s) for s in config["chip_samples"]}
NARROW_SAMPLES = [sid for sid in CHIP_SAMPLES if PEAK_MODE[sid] == "narrow"]
BROAD_SAMPLES = [sid for sid in CHIP_SAMPLES if PEAK_MODE[sid] == "broad"]
FACTOR = {s["id"]: sample_factor(s) for s in config["chip_samples"]}
FACTOR_TO_SAMPLES = defaultdict(list)
FACTOR_TO_CONDITION_SAMPLES = defaultdict(lambda: defaultdict(list))
for sample_id in CHIP_SAMPLES:
    FACTOR_TO_SAMPLES[FACTOR[sample_id]].append(sample_id)
for sample in config["chip_samples"]:
    FACTOR_TO_CONDITION_SAMPLES[sample_factor(sample)][str(sample["condition"]).strip()].append(
        sample["id"]
    )
FACTOR_NAMES = list(FACTOR_TO_SAMPLES)
FACTOR_SLUG = {factor: factor_slug(factor) for factor in FACTOR_NAMES}
SLUG_TO_FACTOR = {slug: factor for factor, slug in FACTOR_SLUG.items()}
FACTOR_SLUGS = list(SLUG_TO_FACTOR)
FACTOR_PEAK_MODE = {
    factor: PEAK_MODE[sample_ids[0]] for factor, sample_ids in FACTOR_TO_SAMPLES.items()
}
NARROW_FACTOR_SLUGS = [
    slug for slug in FACTOR_SLUGS
    if FACTOR_PEAK_MODE[SLUG_TO_FACTOR[slug]] == "narrow"
]
FACTOR_CONDITION_SLUG_TO_NAME = {
    factor: {factor_slug(condition): condition for condition in condition_samples}
    for factor, condition_samples in FACTOR_TO_CONDITION_SAMPLES.items()
}
CONDITION_SLUGS = sorted(
    {
        condition_slug
        for mapping in FACTOR_CONDITION_SLUG_TO_NAME.values()
        for condition_slug in mapping
    }
)


def peak_ext(sample_id):
    """MACS emits broadPeak with --broad and narrowPeak without it."""
    return "broadPeak" if PEAK_MODE[sample_id] == "broad" else "narrowPeak"


def peak_file(sample_id):
    return f"results/peaks/{sample_id}_peaks.{peak_ext(sample_id)}"


def bowtie2_index_files(_wildcards):
    prefix = str(REF["bt2_index"])
    small = [
        f"{prefix}.1.bt2", f"{prefix}.2.bt2", f"{prefix}.3.bt2", f"{prefix}.4.bt2",
        f"{prefix}.rev.1.bt2", f"{prefix}.rev.2.bt2",
    ]
    large = [f"{path}l" for path in small]
    large_exists = [Path(path).exists() for path in large]
    small_exists = [Path(path).exists() for path in small]
    if all(large_exists):
        return large
    if all(small_exists):
        return small
    if any(large_exists) or any(small_exists):
        missing_small = [path for path, exists in zip(small, small_exists) if not exists]
        missing_large = [path for path, exists in zip(large, large_exists) if not exists]
        raise WorkflowError(
            "Bowtie2 index is incomplete. Provide all six .bt2 shards or all six "
            ".bt2l shards. Missing small-index files: "
            + ", ".join(missing_small)
            + "; missing large-index files: "
            + ", ".join(missing_large)
        )
    return small


def factor_samples(wildcards):
    return FACTOR_TO_SAMPLES[SLUG_TO_FACTOR[wildcards.factor]]


def factor_peaks(wildcards):
    return [peak_file(sample) for sample in factor_samples(wildcards)]


def factor_condition_samples(wildcards):
    factor = SLUG_TO_FACTOR[wildcards.factor]
    condition = FACTOR_CONDITION_SLUG_TO_NAME[factor][wildcards.condition]
    return FACTOR_TO_CONDITION_SAMPLES[factor][condition]


def factor_condition_peaks(wildcards):
    return [peak_file(sample) for sample in factor_condition_samples(wildcards)]


def factor_condition_consensus(wildcards):
    factor = SLUG_TO_FACTOR[wildcards.factor]
    return [
        f"results/peaks/consensus/{wildcards.factor}/{factor_slug(condition)}.bed"
        for condition in FACTOR_TO_CONDITION_SAMPLES[factor]
    ]


def factor_rpgc_bigwigs(wildcards):
    return [f"results/bigwig/{sample}.rpgc.bw" for sample in factor_samples(wildcards)]


def factor_log2_bigwigs(wildcards):
    return [f"results/bigwig/{sample}.log2ratio.bw" for sample in factor_samples(wildcards)]


def factor_labels(wildcards):
    return " ".join(quote(sample) for sample in factor_samples(wildcards))


def factor_bam_ids(wildcards):
    ordered = []
    for sample in factor_samples(wildcards):
        for sample_id in (sample, CONTROL_MAP[sample]):
            if sample_id not in ordered:
                ordered.append(sample_id)
    return ordered


def factor_bams(wildcards):
    return [f"results/bam/{sample}.sorted.bam" for sample in factor_bam_ids(wildcards)]


def factor_bais(wildcards):
    return [f"results/bam/{sample}.sorted.bam.bai" for sample in factor_bam_ids(wildcards)]


def factor_bam_labels(wildcards):
    return " ".join(quote(sample) for sample in factor_bam_ids(wildcards))


def factor_label(wildcards):
    return quote(SLUG_TO_FACTOR[wildcards.factor])


def matrix_mode(wildcards):
    factor = SLUG_TO_FACTOR[wildcards.factor]
    return DEEPTOOLS.get("factor_modes", {}).get(factor, FACTOR_PEAK_MODE[factor])


def matrix_args(wildcards):
    mode = matrix_mode(wildcards)
    if mode == "narrow" or mode == "reference_point":
        return (
            "reference-point --referencePoint center "
            f"--beforeRegionStartLength {REFERENCE_POINT_UPSTREAM} "
            f"--afterRegionStartLength {REFERENCE_POINT_DOWNSTREAM} "
            f"--binSize {MATRIX_BIN_SIZE}"
        )
    if mode == "broad" or mode == "scale_regions":
        return (
            "scale-regions "
            f"--regionBodyLength {SCALE_REGIONS_BODY_LENGTH} "
            f"--beforeRegionStartLength {SCALE_REGIONS_UPSTREAM} "
            f"--afterRegionStartLength {SCALE_REGIONS_DOWNSTREAM} "
            f"--binSize {MATRIX_BIN_SIZE}"
        )
    raise WorkflowError(f"Unsupported deepTools matrix mode for {wildcards.factor}: {mode}")


def plot_axis_args(wildcards):
    return (
        "--refPointLabel 'peak center'"
        if matrix_mode(wildcards) in {"narrow", "reference_point"}
        else "--startLabel 'peak start' --endLabel 'peak end'"
    )


def bamcompare_normalization_args(_wildcards):
    if LOG2_SCALE_METHOD == "None":
        return "--scaleFactorsMethod None --normalizeUsing CPM"
    return f"--scaleFactorsMethod {LOG2_SCALE_METHOD}"


CONFIG_DIGEST = hashlib.sha256(
    json.dumps(config, sort_keys=True, default=str).encode("utf-8")
).hexdigest()

RAW_FASTQC_HTML = expand("results/fastqc/raw/{sample}_{read}_fastqc.html", sample=ALL_SAMPLES, read=["R1", "R2"])
RAW_FASTQC_ZIP = expand("results/fastqc/raw/{sample}_{read}_fastqc.zip", sample=ALL_SAMPLES, read=["R1", "R2"])
TRIMMED_FASTQC_HTML = expand("results/fastqc/trimmed/{sample}_R1_val_1_fastqc.html", sample=ALL_SAMPLES) + expand("results/fastqc/trimmed/{sample}_R2_val_2_fastqc.html", sample=ALL_SAMPLES)
TRIMMED_FASTQC_ZIP = expand("results/fastqc/trimmed/{sample}_R1_val_1_fastqc.zip", sample=ALL_SAMPLES) + expand("results/fastqc/trimmed/{sample}_R2_val_2_fastqc.zip", sample=ALL_SAMPLES)
FASTQ_SCREEN_TEXT = expand("results/contamination/fastq_screen/{sample}_{read}_screen.txt", sample=ALL_SAMPLES, read=["R1", "R2"])
FASTQ_SCREEN_HTML = expand("results/contamination/fastq_screen/{sample}_{read}_screen.html", sample=ALL_SAMPLES, read=["R1", "R2"])
PEAKS = [peak_file(s) for s in CHIP_SAMPLES]
RPGC_BIGWIGS = expand("results/bigwig/{sample}.rpgc.bw", sample=ALL_SAMPLES)
LOG2_BIGWIGS = expand("results/bigwig/{sample}.log2ratio.bw", sample=CHIP_SAMPLES)
FACTOR_CONSENSUS = expand("results/peaks/consensus/{factor}.bed", factor=FACTOR_SLUGS)
DEEPTOOLS_OUTPUTS = []
for factor in FACTOR_SLUGS:
    DEEPTOOLS_OUTPUTS.extend(
        [
            f"results/deeptools/{factor}/matrix.gz",
            f"results/deeptools/{factor}/matrix.tsv",
            f"results/deeptools/{factor}/regions.bed",
            f"results/deeptools/{factor}/heatmap.png",
            f"results/deeptools/{factor}/heatmap.pdf",
            f"results/deeptools/{factor}/profile.png",
            f"results/deeptools/{factor}/profile.pdf",
            f"results/deeptools/{factor}/profile.tsv",
            f"results/deeptools/{factor}/qc/fingerprint.png",
            f"results/deeptools/{factor}/qc/fingerprint.tsv",
            f"results/deeptools/{factor}/qc/fingerprint_metrics.tsv",
            f"results/deeptools/{factor}/qc/peak_summary.npz",
            f"results/deeptools/{factor}/qc/peak_signal.tsv",
            f"results/deeptools/{factor}/qc/spearman_heatmap.png",
            f"results/deeptools/{factor}/qc/spearman.tsv",
            f"results/deeptools/{factor}/qc/pca.png",
            f"results/deeptools/{factor}/qc/pca.tsv",
        ]
    )
DIFFBIND_OUTPUTS = []
MOTIF_OUTPUTS = []
for factor in FACTOR_SLUGS:
    DIFFBIND_OUTPUTS.extend(
        [
            f"results/diffbind/{factor}/diffbind_summary.csv",
            f"results/diffbind/{factor}/diffbind_plots.pdf",
            f"results/diffbind/{factor}/diffbind.rds",
            f"results/diffbind/{factor}/contrast.tsv",
        ]
    )
    if FACTOR_PEAK_MODE[SLUG_TO_FACTOR[factor]] == "narrow":
        MOTIF_OUTPUTS.extend(
            [
                f"results/motifs/{factor}/motif_enrichment.tsv",
                f"results/motifs/{factor}/motif_summary.csv",
                f"results/motifs/{factor}/motif_summary.pdf",
            ]
        )
MULTIQC_REPORT = "results/multiqc/multiqc_report.html"


# Pin sample, factor, and condition wildcards to the validated design labels.
wildcard_constraints:
    sample="|".join(re.escape(s) for s in ALL_SAMPLES),
    factor="|".join(re.escape(s) for s in FACTOR_SLUGS),
    condition="|".join(re.escape(s) for s in CONDITION_SLUGS),
    ext="broadPeak|narrowPeak"
