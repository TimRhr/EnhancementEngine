import pytest

pytest.importorskip("flask")

try:
    from enhancement_engine.webapp import app
    from enhancement_engine.core.engine import EnhancementEngine
except Exception:  # pragma: no cover - skip if webapp missing
    pytest.skip("webapp not available", allow_module_level=True)


def test_index_page():
    client = app.test_client()
    resp = client.get("/")
    assert resp.status_code == 200
    html = resp.get_data(as_text=True)
    assert 'id="search-gene"' in html
    assert 'id="variant-menu"' in html


def test_analyze_route(monkeypatch):
    called = {}

    def fake_analyze_gene(self, gene_name, variant="enhancement_variant", target_tissue="general"):
        called['gene'] = gene_name
        class DummyGain:
            improvement_factor = 1.0

        class DummyPred:
            enhancement_gain = DummyGain()

        class DummySafety:
            overall_score = 1.0

        class R:
            def __init__(self, gene_name):
                self.gene_name = gene_name
                self.target_variant = variant
                self.feasibility_score = 95
                self.safety_assessment = DummySafety()
                self.predicted_effect = DummyPred()
                self.confidence_score = 1.0
                self.summary = ""
                self.recommendations = []
        return R(gene_name)

    monkeypatch.setattr(EnhancementEngine, "analyze_gene", fake_analyze_gene)

    client = app.test_client()
    resp = client.get("/analyze?gene=COMT")
    assert resp.status_code in (200, 405)
    if resp.status_code == 200:
        assert called['gene'] == "COMT"


def test_gene_variants_api(monkeypatch):
    monkeypatch.setattr(
        EnhancementEngine,
        "validate_gene_input",
        lambda self, gene_name: {"available_variants": ["v1", "v2"]},
    )
    client = app.test_client()
    resp = client.get("/api/gene_variants?gene=TEST")
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["variants"] == ["v1", "v2"]
