import pytest
import shutil
from enhancement_engine.core.disease_db import DiseaseDatabaseClient


class DummyHandle:
    def __init__(self, name):
        self.name = name

    def close(self):
        pass


def test_fetch_disease_info_cached(monkeypatch, tmp_path):
    calls = {"esearch": 0, "esummary": 0}

    def fake_esearch(db, term, retmax=1):
        calls["esearch"] += 1
        return DummyHandle("search")

    def fake_esummary(db, id):
        calls["esummary"] += 1
        return DummyHandle("summary")

    def fake_read(handle, validate=False):
        if handle.name == "search":
            return {"IdList": ["1"]}
        return [{"OddsRatio": "2.0", "AlleleFrequency": "0.05", "Penetrance": "0.1"}]

    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esearch", fake_esearch
    )
    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esummary", fake_esummary
    )
    monkeypatch.setattr("enhancement_engine.core.disease_db.Entrez.read", fake_read)

    client = DiseaseDatabaseClient("test@example.com", cache_dir=str(tmp_path))
    info1 = client.fetch_disease_info("GENE", "VAR", "disease")
    assert info1 == {"odds_ratio": 2.0, "allele_frequency": 0.05, "penetrance": 0.1}
    info2 = client.fetch_disease_info("GENE", "VAR", "disease")
    assert info2 == info1
    assert calls["esearch"] == 1
    assert calls["esummary"] == 1


def test_risk_calculator_dynamic(monkeypatch):
    from enhancement_engine.core.disease_risk import DiseaseRiskCalculator

    class DummyClient:
        def fetch_disease_info(self, gene, variant, disease=None):
            return {"odds_ratio": 1.5, "allele_frequency": 0.1, "penetrance": 0.05}

    calc = DiseaseRiskCalculator(DummyClient())
    patient = {"ethnicity": "european", "age": 40, "sex": "female"}
    risk = calc.calculate_disease_risk("UNKNOWN", "VAR", patient, "rare")
    assert pytest.approx(risk.odds_ratio) == 1.5
    assert risk.population_frequency == 0.1
    assert risk.penetrance > 0


def test_disease_db_client_methods(monkeypatch, tmp_path):
    """Ensure key methods exist and are callable."""

    def fake_esearch(db, term, retmax=20):
        return DummyHandle("search")

    def fake_esummary(db, id):
        return DummyHandle("summary")

    def fake_read(handle, validate=False):
        return {"IdList": []} if handle.name == "search" else []

    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esearch",
        fake_esearch,
    )
    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esummary",
        fake_esummary,
    )
    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.read",
        fake_read,
    )

    client = DiseaseDatabaseClient("test@example.com", cache_dir=str(tmp_path))
    assert hasattr(client, "search_diseases")
    assert callable(client.search_diseases)
    assert hasattr(client, "fetch_associated_genes")
    assert callable(client.fetch_associated_genes)


def test_search_diseases_handles_dict_summary(monkeypatch, tmp_path):
    """search_diseases should handle dict summary structures."""

    def fake_esearch(db, term, retmax=20):
        return DummyHandle("search")

    def fake_esummary(db, id):
        return DummyHandle("summary")

    def fake_read(handle, validate=False):
        if handle.name == "search":
            return {"IdList": ["1"]}
        return {
            "DocumentSummarySet": {"DocumentSummary": [{"Description": "Test Disease"}]}
        }

    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esearch", fake_esearch
    )
    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esummary", fake_esummary
    )
    monkeypatch.setattr("enhancement_engine.core.disease_db.Entrez.read", fake_read)

    client = DiseaseDatabaseClient("test@example.com", cache_dir=str(tmp_path))
    result = client.search_diseases("test")
    assert result == ["test disease"]


def test_search_diseases_deduplicates(monkeypatch, tmp_path):
    """search_diseases should remove duplicate names."""

    def fake_esearch(db, term, retmax=20):
        return DummyHandle("search")

    def fake_esummary(db, id):
        return DummyHandle("summary")

    def fake_read(handle, validate=False):
        if handle.name == "search":
            return {"IdList": ["1", "2"]}
        return {
            "DocumentSummarySet": {"DocumentSummary": [{"Description": "Dup"}]}
        }

    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esearch", fake_esearch
    )
    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esummary", fake_esummary
    )
    monkeypatch.setattr("enhancement_engine.core.disease_db.Entrez.read", fake_read)

    client = DiseaseDatabaseClient("test@example.com", cache_dir=str(tmp_path))
    result = client.search_diseases("dup")
    assert result == ["dup"]


def test_fetch_associated_genes_creates_cache_dir(monkeypatch, tmp_path):
    """fetch_associated_genes should create cache dir and parse dict summary."""

    def fake_esearch(db, term, retmax=20):
        return DummyHandle("search")

    def fake_esummary(db, id):
        return DummyHandle("summary")

    def fake_read(handle, validate=False):
        if handle.name == "search":
            return {"IdList": ["1"]}
        return {
            "DocumentSummarySet": {
                "DocumentSummary": [{"gene": "GENE1", "variation_name": "VAR1"}]
            }
        }

    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esearch", fake_esearch
    )
    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.esummary", fake_esummary
    )
    monkeypatch.setattr("enhancement_engine.core.disease_db.Entrez.read", fake_read)

    client = DiseaseDatabaseClient("test@example.com", cache_dir=str(tmp_path))
    # remove disease cache directory to ensure method recreates it
    cache_dir = tmp_path / "disease"
    if cache_dir.exists():
        shutil.rmtree(cache_dir)

    result = client.fetch_associated_genes("disease")
    assert result == {"GENE1": ["VAR1"]}
    assert (tmp_path / "disease").exists()


def test_rate_limit(monkeypatch):
    """_rate_limit should pause when called rapidly."""

    times = iter([0.0, 0.0, 0.1, 0.1])
    sleep_calls = []

    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.time.time", lambda: next(times)
    )
    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.time.sleep",
        lambda s: sleep_calls.append(s),
    )

    client = DiseaseDatabaseClient("test@example.com")
    client._rate_limit()
    client._rate_limit()

    assert sleep_calls[0] == pytest.approx(client.request_delay)
    assert sleep_calls[1] == pytest.approx(client.request_delay - 0.1)


def test_init_sets_api_key(monkeypatch):
    monkeypatch.setattr(
        "enhancement_engine.core.disease_db.Entrez.api_key", None, raising=False
    )
    DiseaseDatabaseClient("test@example.com", api_key="key")
    from enhancement_engine.core.disease_db import Entrez

    assert Entrez.api_key == "key"
