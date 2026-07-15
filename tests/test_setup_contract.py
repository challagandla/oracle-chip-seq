from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_miniforge_release_asset_is_versioned_and_linux_scoped():
    setup = (ROOT / "setup.sh").read_text()
    assert 'asset="Miniforge3-${MINIFORGE_VERSION}-Linux-${ARCHITECTURE}.sh"' in setup
    assert '[[ "$SYSTEM" == "Linux" ]]' in setup
    assert "x86_64)" in setup
    assert "aarch64|arm64)" in setup
    assert "MacOSX" not in setup
    assert "ppc64le" not in setup


def test_every_entrypoint_discovers_custom_miniforge_prefix():
    custom_candidate = '"${MINIFORGE_HOME:+${MINIFORGE_HOME%/}/bin/conda}"'
    for script in ("setup.sh", "run.sh", "prepare_references.sh", "make_config.sh"):
        source = (ROOT / script).read_text()
        assert custom_candidate in source
        assert source.index(custom_candidate) < source.index('"$(type -P conda')
