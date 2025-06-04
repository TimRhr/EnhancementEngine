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


def test_analyze_route(monkeypatch):
    called = {}

    def fake_analyze_gene(self, gene_name, variant="enhancement_variant", target_tissue="general"):
        called['gene'] = gene_name
        class R:
            def __init__(self, gene_name):
                self.gene_name = gene_name
        return R(gene_name)

    monkeypatch.setattr(EnhancementEngine, "analyze_gene", fake_analyze_gene)

    client = app.test_client()
    resp = client.get("/analyze?gene=COMT")
    assert resp.status_code == 200
    assert called['gene'] == "COMT"
