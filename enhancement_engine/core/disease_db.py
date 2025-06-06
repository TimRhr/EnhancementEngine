from __future__ import annotations
"""ClinVar/MedGen client for disease information."""

import logging
from typing import Dict, Optional, List, Set

from .database import CacheManager

try:
    from Bio import Entrez
except ImportError:  # pragma: no cover - biopython required
    Entrez = None  # type: ignore


class DiseaseDatabaseClient:
    """Query ClinVar/MedGen for disease risk data."""

    def __init__(self, email: str, cache_dir: str = "data/cache/disease_db"):
        if Entrez is None:
            raise ImportError("Biopython is required for DiseaseDatabaseClient")
        Entrez.email = email
        self.logger = logging.getLogger(__name__)
        self.cache = CacheManager(cache_dir)

    def fetch_disease_info(
        self, gene: str, variant: str, disease: Optional[str] = None
    ) -> Optional[Dict[str, float]]:
        """Fetch odds ratio, allele frequency and penetrance."""
        cache_key = f"{gene}_{variant}_{disease}"
        cached = self.cache.get(cache_key, "disease")
        if cached:
            return cached

        # ensure disease cache directory exists
        (self.cache.cache_dir / "disease").mkdir(exist_ok=True)

        query_parts = [gene]
        if variant:
            query_parts.append(variant)
        if disease:
            query_parts.append(disease)
        query = " AND ".join(query_parts)

        try:
            search_handle = Entrez.esearch(db="clinvar", term=query, retmax=1)
            search_result = Entrez.read(search_handle)
            search_handle.close()
            ids = search_result.get("IdList", [])
            if not ids:
                return None
            summary_handle = Entrez.esummary(db="clinvar", id=ids[0])
            summary = Entrez.read(summary_handle)
            summary_handle.close()
            doc = summary[0] if summary else {}
            data = {
                "odds_ratio": float(doc.get("OddsRatio", 1.0)),
                "allele_frequency": float(doc.get("AlleleFrequency", 0.0)),
                "penetrance": float(doc.get("Penetrance", 0.1)),
            }
            self.cache.set(cache_key, data, "disease")
            return data
        except Exception as exc:  # pragma: no cover - network errors
            self.logger.warning(f"Failed to query ClinVar: {exc}")
            return None

    def search_diseases(self, term: str, max_results: int = 20) -> List[str]:
        """Search MedGen for disease names matching the term."""
        cache_key = f"search_{term}"
        cached = self.cache.get(cache_key, "disease")
        if cached:
            return cached

        try:
            handle = Entrez.esearch(db="medgen", term=term, retmax=max_results)
            result = Entrez.read(handle)
            handle.close()
            ids = result.get("IdList", [])
            names = []
            for did in ids:
                sum_handle = Entrez.esummary(db="medgen", id=did)
                summary = Entrez.read(sum_handle)
                sum_handle.close()
                if summary:
                    name = (
                        summary[0].get("Description")
                        or summary[0].get("Title")
                        or summary[0].get("Name")
                    )
                    if name:
                        normalized = name.strip().lower()
                        if normalized:
                            names.append(normalized)
            self.cache.set(cache_key, names, "disease")
            return names
        except Exception as exc:  # pragma: no cover - network errors
            self.logger.warning(f"Failed to search diseases: {exc}")
            return []

    def fetch_associated_genes(
        self, disease: str, max_results: int = 20
    ) -> Dict[str, List[str]]:
        """Return genes and variants linked to a disease."""
        cache_key = f"assoc_{disease}"
        cached = self.cache.get(cache_key, "disease")
        if cached:
            return cached

        genes: Dict[str, Set[str]] = {}

        def query(term: str) -> List[str]:
            handle = Entrez.esearch(db="clinvar", term=term, retmax=max_results)
            result = Entrez.read(handle)
            handle.close()
            return result.get("IdList", [])

        def collect(ids: List[str]) -> None:
            for cid in ids:
                sum_handle = Entrez.esummary(db="clinvar", id=cid)
                summary = Entrez.read(sum_handle)
                sum_handle.close()
                if not summary:
                    continue
                doc = summary[0]
                gene = (
                    doc.get("gene")
                    or doc.get("gene_symbol")
                    or doc.get("geneName")
                    or ""
                )
                variant = doc.get("variation_name") or doc.get("title") or ""
                gene = str(gene).strip()
                variant = str(variant).strip()
                if not gene:
                    continue
                genes.setdefault(gene, set())
                if variant:
                    genes[gene].add(variant)

        try:
            ids = query(disease)
            if not ids:
                syns = self.search_diseases(disease, max_results=max_results)
                for s in syns:
                    s_ids = query(s)
                    result = {}
                    if s_ids:
                        collect(s_ids)
                        result = {g: sorted(v) for g, v in genes.items()}
                    self.cache.set(f"assoc_{s}", result, "disease")
                    if result:
                        break
            else:
                collect(ids)

            result = {g: sorted(v) for g, v in genes.items()}
            self.cache.set(cache_key, result, "disease")
            return result
        except Exception as exc:  # pragma: no cover - network errors
            self.logger.warning(f"Failed to fetch associated genes: {exc}")
            return {}
