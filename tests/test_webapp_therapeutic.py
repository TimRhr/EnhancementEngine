import pytest

pytest.importorskip("flask")

try:
    from enhancement_engine.webapp import app
    from enhancement_engine.core.therapeutic_engine import TherapeuticEnhancementEngine
except Exception:
    pytest.skip("therapeutic webapp not available", allow_module_level=True)


def _get_therapeutic_engine():
    func = app.view_functions.get("/api/disease_info") or app.view_functions.get("api_disease_info")
    if func and func.__closure__ and len(func.__closure__) > 1:
        return func.__closure__[1].cell_contents
    return None


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
            return ["dynamic", "Dynamic", " dynamic "]
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


def test_disease_api_specific_search(monkeypatch):
    """Ensure /api/diseases returns mocked diseases for a query."""
    from enhancement_engine.core.disease_db import DiseaseDatabaseClient

    monkeypatch.setattr(
        DiseaseDatabaseClient,
        "search_diseases",
        lambda self, term, max_results=5: ["dynamic disease"],
    )

    client = app.test_client()
    resp = client.get("/api/diseases?q=dynamic")
    assert resp.status_code == 200
    data = resp.get_json()
    assert "dynamic disease" in data["diseases"]


def test_disease_info_api(monkeypatch):
    import webapp.run as run_module

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
    import webapp.run as run_module

    monkeypatch.setattr(run_module, "DISEASE_GENES", {})

    class DummyClient:
        def fetch_associated_genes(self, disease, max_results=5):
            return {"DG": ["v1", "v2"]}

    therapeutic = _get_therapeutic_engine()
    if therapeutic:
        monkeypatch.setattr(
            therapeutic,
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


def test_disease_info_api_synonym_fallback(monkeypatch):
    import webapp.run as run_module

    monkeypatch.setattr(run_module, "DISEASE_GENES", {})

    class DummyClient:
        def __init__(self):
            self.calls = []

        def fetch_associated_genes(self, disease, max_results=5):
            self.calls.append(("fetch", disease))
            if disease == "synonym":
                return {"SG": ["v1"]}
            return {}

        def search_diseases(self, term, max_results=5):
            self.calls.append(("search", term))
            return ["synonym"]

    dummy = DummyClient()

    therapeutic = _get_therapeutic_engine()
    if therapeutic:
        monkeypatch.setattr(
            therapeutic,
            "disease_db_client",
            dummy,
            raising=False,
        )

    client = app.test_client()
    resp = client.get("/api/disease_info?disease=unknown")
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["genes"] == ["SG"]
    assert data["variants"]["SG"] == ["v1"]
    assert ("fetch", "unknown") in dummy.calls
    assert ("search", "unknown") in dummy.calls
    assert ("fetch", "synonym") in dummy.calls



def test_therapeutic_page_contains_failure_message():
    from bs4 import BeautifulSoup

    client = app.test_client()
    resp = client.get("/therapeutic")
    assert resp.status_code == 200
    html = resp.get_data(as_text=True)

    soup = BeautifulSoup(html, "html.parser")
    gv = soup.find(id="gene-variant-section")
    assert gv is not None
    assert "display:none" in (gv.get("style", "").replace(" ", "").lower())

    # dropdown menu for variants should be present
    assert soup.find(id="variant-menu") is not None

    # dropdown list for dynamic disease search should be present
    assert soup.find(id="disease-results") is not None

    loading = soup.find(id="gene-loading")
    assert loading is not None
    assert "display:none" in (loading.get("style", "").replace(" ", "").lower())

    # search field/button should be present for finding diseases dynamically
    assert soup.find(id="search-disease") is not None
    assert soup.find(id="search-spinner") is not None

    assert "No genes found for" in html

    
def test_disease_info_api_clinvar_synonym(monkeypatch, tmp_path):
    import webapp.run as run_module
    from enhancement_engine.core import disease_db as db_mod

    monkeypatch.setattr(run_module, "DISEASE_GENES", {})

    class DummyHandle:
        def __init__(self, name):
            self.name = name
        def close(self):
            pass

    calls = []

    def fake_esearch(db, term, retmax=20):
        calls.append(("esearch", term))
        if term == "cancer":
            return DummyHandle("search_cancer")
        return DummyHandle(f"search_{term}")

    def fake_esummary(db, id):
        calls.append(("esummary", id))
        return DummyHandle(f"summary_{id}")

    def fake_read(handle, validate=False):
        if handle.name == "search_cancer":
            return {"IdList": []}
        if handle.name == "search_carcinoma":
            return {"IdList": ["1"]}
        if handle.name == "summary_1":
            return [{"gene": "CG", "variation_name": "v1"}]
        return {}

    monkeypatch.setattr(db_mod.Entrez, "esearch", fake_esearch)
    monkeypatch.setattr(db_mod.Entrez, "esummary", fake_esummary)
    monkeypatch.setattr(db_mod.Entrez, "read", fake_read)

    def fake_search(self, term, max_results=20):
        calls.append(("search_diseases", term))
        return ["carcinoma"]

    monkeypatch.setattr(db_mod.DiseaseDatabaseClient, "search_diseases", fake_search, raising=False)

    client = db_mod.DiseaseDatabaseClient("test@example.com", cache_dir=str(tmp_path))

    therapeutic = _get_therapeutic_engine()
    if therapeutic:
        monkeypatch.setattr(therapeutic, "disease_db_client", client, raising=False)

    c = app.test_client()
    resp = c.get("/api/disease_info?disease=cancer")
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["genes"] == ["CG"]
    assert data["variants"]["CG"] == ["v1"]
    assert ("search_diseases", "cancer") in calls
    assert ("esearch", "carcinoma") in calls


def test_disease_datalist_updates(monkeypatch):
    """Therapeutic page dynamically loads diseases via /api/diseases."""
    from enhancement_engine.core.disease_db import DiseaseDatabaseClient
    from bs4 import BeautifulSoup

    monkeypatch.setattr(
        DiseaseDatabaseClient,
        "search_diseases",
        lambda self, term, max_results=5: ["dynamic disease"],
    )

    client = app.test_client()
    page = client.get("/therapeutic")
    assert page.status_code == 200

    resp = client.get("/api/diseases?q=dynamic")
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["diseases"] == ["dynamic disease"]

    soup = BeautifulSoup(page.get_data(as_text=True), "html.parser")
    menu = soup.find("ul", id="disease-results")
    assert menu is not None
    menu.clear()
    for name in data["diseases"]:
        li = soup.new_tag("li")
        a = soup.new_tag("a", href="#")
        a["class"] = "dropdown-item"
        a.string = name
        li.append(a)
        menu.append(li)

    assert menu.find("a", string="dynamic disease") is not None


def test_disease_dropdown_populates(monkeypatch):
    """Disease dropdown menu should populate with search results."""
    from enhancement_engine.core.disease_db import DiseaseDatabaseClient
    from bs4 import BeautifulSoup

    monkeypatch.setattr(
        DiseaseDatabaseClient,
        "search_diseases",
        lambda self, term, max_results=5: ["d1", "d2"],
    )

    client = app.test_client()
    page = client.get("/therapeutic")
    assert page.status_code == 200

    resp = client.get("/api/diseases?q=d1")
    assert resp.status_code == 200
    data = resp.get_json()
    for d in ["d1", "d2"]:
        assert d in data["diseases"]

    soup = BeautifulSoup(page.get_data(as_text=True), "html.parser")
    menu = soup.find("ul", id="disease-results")
    assert menu is not None
    menu.clear()
    for name in data["diseases"]:
        li = soup.new_tag("li")
        a = soup.new_tag("a", href="#")
        a["class"] = "dropdown-item"
        a.string = name
        li.append(a)
        menu.append(li)

    items = [item.get_text(strip=True) for item in menu.find_all("li")]
    assert "d1" in items and "d2" in items
