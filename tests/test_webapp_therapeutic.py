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
        if term:
            return ["Dynamic", "dynamic", " dynamic "]
        return []

    monkeypatch.setattr(DiseaseDatabaseClient, "search_diseases", fake_search)

    client = app.test_client()
    resp = client.get("/api/diseases")
    assert resp.status_code == 200
    data = resp.get_json()
    assert "rheumatoid_arthritis" in data["diseases"]

    resp = client.get("/api/diseases?q=dyn")
    data = resp.get_json()
    assert "dynamic" in data["diseases"]
    lower = [d.lower() for d in data["diseases"]]
    assert len(lower) == len(set(lower))


def test_disease_info_api(monkeypatch):
    from enhancement_engine.webapp import run as run_module

    fake = {
        "cat": {
            "G1": {
                "disease_associations": {"dis": {}},
                "pathogenic_variants": {"v1": {}, "v2": {}}
            },
            "G2": {
                "disease_associations": {"other": {}},
                "pathogenic_variants": {"x": {}}
            }
        }
    }

    monkeypatch.setattr(run_module, "DISEASE_GENES", fake)

    client = app.test_client()
    resp = client.get("/api/disease_info?disease=dis")
    assert resp.status_code == 200
    data = resp.get_json()
    assert set(data["genes"]) == {"G1"}
    assert set(data["variants"]["G1"]) == {"v1", "v2"}


def test_disease_info_api_dynamic(monkeypatch):
    from enhancement_engine.webapp import run as run_module

    monkeypatch.setattr(run_module, "DISEASE_GENES", {})

    class DummyClient:
        def fetch_associated_genes(self, disease, max_results=5):
            return {"DG": ["v1", "v2"]}

    if getattr(run_module, "therapeutic_engine", None):
        monkeypatch.setattr(
            run_module.therapeutic_engine,
            "disease_db_client",
            DummyClient(),
            raising=False,
        )

    client = app.test_client()
    resp = client.get("/api/disease_info?disease=new")
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["genes"] == ["DG"]
    assert data["variants"]["DG"] == ["v1", "v2"]
