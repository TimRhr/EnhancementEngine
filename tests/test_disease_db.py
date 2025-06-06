import pytest
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

    monkeypatch.setattr("enhancement_engine.core.disease_db.Entrez.esearch", fake_esearch)
    monkeypatch.setattr("enhancement_engine.core.disease_db.Entrez.esummary", fake_esummary)
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
