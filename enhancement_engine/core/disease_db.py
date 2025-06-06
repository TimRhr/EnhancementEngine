from __future__ import annotations

"""ClinVar/MedGen client for disease information."""

import logging
import time
from typing import Dict, Optional, List, Set

from .database import CacheManager

try:
    from Bio import Entrez
except ImportError:  # pragma: no cover - biopython required
    Entrez = None  # type: ignore


class DiseaseDatabaseClient:
    """Query ClinVar/MedGen for disease risk data."""

    def __init__(
        self,
        email: str,
        cache_dir: str = "data/cache/disease_db",
        api_key: Optional[str] = None,
        request_delay: float = 0.34,
    ):
        if Entrez is None:
            raise ImportError("Biopython is required for DiseaseDatabaseClient")

        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        self.request_delay = request_delay
        self._last_request = 0.0

        self.logger = logging.getLogger(__name__)
        self.cache = CacheManager(cache_dir)

    def _read(self, handle):
        """Wrapper around ``Entrez.read`` disabling XML validation."""
        return Entrez.read(handle, validate=False)

    def _rate_limit(self) -> None:
        """Enforce rate limiting for NCBI requests."""
        current_time = time.time()
        elapsed = current_time - self._last_request

        if elapsed < self.request_delay:
            time.sleep(self.request_delay - elapsed)

        self._last_request = time.time()

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
            self._rate_limit()
            search_handle = Entrez.esearch(db="clinvar", term=query, retmax=1)
            search_result = self._read(search_handle)
            search_handle.close()
            ids = search_result.get("IdList", [])
            if not ids:
                return None
            self._rate_limit()
            summary_handle = Entrez.esummary(db="clinvar", id=ids[0])
            summary = self._read(summary_handle)
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
            self._rate_limit()
            handle = Entrez.esearch(db="medgen", term=term, retmax=max_results)
            result = self._read(handle)
            handle.close()
            ids = result.get("IdList", [])
            names = []
            for did in ids:
                self._rate_limit()
                sum_handle = Entrez.esummary(db="medgen", id=did)
                summary = self._read(sum_handle)
                sum_handle.close()

                if isinstance(summary, list):
                    doc = summary[0] if summary else None
                else:
                    doc = summary.get("DocumentSummarySet", {}).get(
                        "DocumentSummary", []
                    )
                    doc = doc[0] if isinstance(doc, list) and doc else None

                if doc:
                    name = doc.get("Description") or doc.get("Title") or doc.get("Name")
                    if name:
                        normalized = name.strip().lower()
                        if normalized:
                            names.append(normalized)

            unique: List[str] = []
            seen: Set[str] = set()
            for n in names:
                if n not in seen:
                    unique.append(n)
                    seen.add(n)
            names = unique

            self.cache.set(cache_key, names, "disease")
            return names
        except Exception as exc:
            # Verbesserte Fehlerbehandlung
            self.logger.warning(
                f"Failed to search diseases for term '{term}': {type(exc).__name__}: {exc}"
            )
            if "email" in str(exc).lower():
                self.logger.error(
                    "NCBI requires a valid email address. Please check your DEMO_EMAIL setting."
                )
            return []

    def fetch_associated_genes(
        self, disease: str, max_results: int = 20
    ) -> Dict[str, List[str]]:
        """Return genes and variants linked to a disease."""
        cache_key = f"assoc_{disease}"
        cached = self.cache.get(cache_key, "disease")
        if cached:
            return cached

        # ensure disease cache directory exists
        (self.cache.cache_dir / "disease").mkdir(exist_ok=True)

        genes: Dict[str, Set[str]] = {}

        def query(term: str) -> List[str]:
            try:
                self._rate_limit()
                handle = Entrez.esearch(db="clinvar", term=term, retmax=max_results)
                result = self._read(handle)
                handle.close()
                return result.get("IdList", [])
            except Exception as e:
                self.logger.warning(f"Query failed for term '{term}': {e}")
                return []

        def collect(ids: List[str]) -> None:
            for cid in ids:
                try:
                    self._rate_limit()
                    sum_handle = Entrez.esummary(db="clinvar", id=cid)
                    summary = self._read(sum_handle)
                    sum_handle.close()
                    if isinstance(summary, list):
                        doc = summary[0] if summary else None
                    else:
                        doc = summary.get("DocumentSummarySet", {}).get(
                            "DocumentSummary", []
                        )
                        doc = doc[0] if isinstance(doc, list) and doc else None
                    if not doc:
                        continue
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
                except Exception as e:
                    self.logger.warning(f"Failed to process ClinVar ID {cid}: {e}")
                    continue

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
        except Exception as exc:
            # Verbesserte Fehlerbehandlung
            self.logger.warning(
                f"Failed to fetch associated genes for disease '{disease}': {type(exc).__name__}: {exc}"
            )
            if "email" in str(exc).lower():
                self.logger.error(
                    "NCBI requires a valid email address. Please check your DEMO_EMAIL setting."
                )
            return {}
