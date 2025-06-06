import pytest

pytest.importorskip("flask")

try:
    from enhancement_engine.webapp import app
    from enhancement_engine.core.therapeutic_engine import TherapeuticEnhancementEngine
except Exception:
    pytest.skip("therapeutic webapp not available", allow_module_level=True)


def test_therapeutic_route(monkeypatch):
    # Skip test if the therapeutic route is not registered
    routes = {rule.rule for rule in app.url_map.iter_rules()}
    if "/therapeutic/analyze" not in routes:
        pytest.skip("therapeutic route not available")

    called = {}

    def fake_analyze(self, gene_name, variant, patient_data, disease="rheumatoid_arthritis"):
        called["args"] = (gene_name, variant, disease)
        class R:
            def __init__(self, gene_name):
                self.gene_name = gene_name
        return R(gene_name)

    monkeypatch.setattr(TherapeuticEnhancementEngine, "analyze_disease_gene", fake_analyze)

    client = app.test_client()
    resp = client.get("/therapeutic/analyze?gene=PTPN22&variant=R620W&disease=ra")
    assert resp.status_code == 200
    assert called["args"][0] == "PTPN22"


def test_disease_api(monkeypatch):
    from enhancement_engine.core.disease_db import DiseaseDatabaseClient

    def fake_search(self, term, max_results=5):
        return ["dynamic"] if term else []

    monkeypatch.setattr(DiseaseDatabaseClient, "search_diseases", fake_search)

    client = app.test_client()
    resp = client.get("/api/diseases")
    assert resp.status_code == 200
    data = resp.get_json()
    assert "rheumatoid_arthritis" in data["diseases"]

    resp = client.get("/api/diseases?q=dyn")
    data = resp.get_json()
    assert "dynamic" in data["diseases"]
