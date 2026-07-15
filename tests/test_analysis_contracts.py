from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_motif_enrichment_uses_seeded_blacklist_excluded_background():
    script = (ROOT / "analysis" / "motif_enrichment.R").read_text()
    assert 'background = "otherBins"' in script
    assert "set.seed(seed)" in script
    assert "SerialParam(RNGseed = seed)" in script
    assert "blacklist_regions" in script
    assert "IRanges::setdiff(IRanges::IRanges(1L, chromosome_length), excluded)" in script
    assert "width(foreground_regions) == window_bp" in script
    assert "anyDuplicated(background_ranges)" in script
    assert 'match("foreground", colnames(enrichment))' in script
    assert "P.adjust = adjusted_p" in script
    assert "filter(is.finite" not in script
    assert "filter(!is.na(negLog10Padj)" in script
    assert "if_else(is.infinite(negLog10Padj)" in script


def test_motif_rule_is_factor_scoped_and_declares_complete_outputs():
    rules = (ROOT / "workflow" / "analysis.smk").read_text()
    motif = rules.split("rule motif_enrichment:", 1)[1]
    assert 'peaks="results/peaks/consensus/{factor}.bed"' in motif
    assert "genome_index=lambda wc: f\"{REF['genome']}.fai\"" in motif
    assert 'blacklist=lambda wc: REF["black_list"]' in motif
    assert 'table="results/motifs/{factor}/motif_enrichment.tsv"' in motif
    assert 'csv="results/motifs/{factor}/motif_summary.csv"' in motif
    assert 'pdf="results/motifs/{factor}/motif_summary.pdf"' in motif
    assert "NARROW_FACTOR_SLUGS" in motif
    assert "MOTIF_BACKGROUND_MULTIPLIER" in motif
    assert "MOTIF_SEED" in motif


def test_diffbind_tracks_bam_indexes_as_inputs():
    rules = (ROOT / "workflow" / "analysis.smk").read_text()
    run_diffbind = rules.split("rule run_diffbind:", 1)[1].split("rule motif_enrichment:", 1)[0]
    assert "bams=factor_bams" in run_diffbind
    assert "bais=factor_bais" in run_diffbind
