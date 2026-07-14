from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_homer_summary_preserves_hash_columns_and_extreme_significance():
    script = (ROOT / "analysis" / "motif_summary.R").read_text()
    assert 'comment = "#"' not in script
    assert '"Log.P.value" %in% colnames(motifs)' in script
    assert "logP = -Log.P.value / log(10)" in script


def test_diffbind_tracks_bam_indexes_as_inputs():
    rules = (ROOT / "workflow" / "analysis.smk").read_text()
    run_diffbind = rules.split("rule run_diffbind:", 1)[1].split("rule motif_summary:", 1)[0]
    assert "bams=factor_bams" in run_diffbind
    assert "bais=factor_bais" in run_diffbind
