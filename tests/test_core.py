import pytest


def test_engine_initialization_valid_email():
    pytest.importorskip("Bio")
    pytest.importorskip("numpy")
    pytest.importorskip("pandas")
    from enhancement_engine.core.engine import EnhancementEngine
    engine = EnhancementEngine("user@example.com")
    assert engine.email == "user@example.com"
    assert engine.config.researcher_email == "user@example.com"


def test_engine_initialization_invalid_email():
    pytest.importorskip("Bio")
    pytest.importorskip("numpy")
    pytest.importorskip("pandas")
    from enhancement_engine.core.engine import EnhancementEngine, EnhancementEngineError
    with pytest.raises(EnhancementEngineError):
        EnhancementEngine("invalid-email")
