"""
Database and data access layer for Enhancement Engine.

This module handles all external data sources including:
- NCBI databases (Gene, Nucleotide, dbSNP, ClinVar)
- UniProt protein database
- Local caching system
- Gene annotation data
"""

import os
import json
import pickle
import hashlib
import time
from typing import Dict, List, Optional, Union, Any, Tuple
from pathlib import Path
from datetime import datetime, timedelta
import logging

from ..utils import is_valid_email

try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError("BioPython is required. Install with: pip install biopython")

import requests
import pandas as pd

from ..models.data_classes import GeneInfo, VariantInfo, EnhancementCategory, CasType
from ..models.constants import (
    ENHANCEMENT_GENES,
    NCBI_DATABASES,
    EXTERNAL_RESOURCES,
    ERROR_MESSAGES,
)


class DatabaseError(Exception):
    """Custom exception for database-related errors."""

    pass


class CacheManager:
    """Manages local caching of database queries."""

    def __init__(self, cache_dir: str = "data/cache", max_age_days: int = 30):
        """
        Initialize cache manager.

        Args:
            cache_dir: Directory for cache files
            max_age_days: Maximum age of cached files in days
        """
        self.cache_dir = Path(cache_dir)
        self.max_age = timedelta(days=max_age_days)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        (self.cache_dir / "sequences").mkdir(exist_ok=True)
        (self.cache_dir / "genes").mkdir(exist_ok=True)
        (self.cache_dir / "variants").mkdir(exist_ok=True)
        (self.cache_dir / "disease").mkdir(exist_ok=True)

        self.logger = logging.getLogger(__name__)

    def _get_cache_path(self, key: str, category: str) -> Path:
        """Generate cache file path for a given key and category."""
        # Create hash of key for filename
        key_hash = hashlib.md5(key.encode()).hexdigest()
        return self.cache_dir / category / f"{key_hash}.pkl"

    def get(self, key: str, category: str = "general") -> Optional[Any]:
        """
        Retrieve item from cache.

        Args:
            key: Cache key
            category: Cache category (sequences, genes, variants, etc.)

        Returns:
            Cached object or None if not found/expired
        """
        cache_path = self._get_cache_path(key, category)

        if not cache_path.exists():
            return None

        # Check if cache is expired
        file_age = datetime.now() - datetime.fromtimestamp(cache_path.stat().st_mtime)
        if file_age > self.max_age:
            self.logger.debug(f"Cache expired for key: {key}")
            cache_path.unlink()  # Remove expired cache
            return None

        try:
            with open(cache_path, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            self.logger.warning(f"Failed to load cache for {key}: {e}")
            cache_path.unlink()  # Remove corrupted cache
            return None

    def set(self, key: str, value: Any, category: str = "general") -> None:
        """
        Store item in cache.

        Args:
            key: Cache key
            value: Object to cache
            category: Cache category
        """
        cache_path = self._get_cache_path(key, category)

        try:
            with open(cache_path, "wb") as f:
                pickle.dump(value, f)
            self.logger.debug(f"Cached {category}/{key}")
        except Exception as e:
            self.logger.warning(f"Failed to cache {key}: {e}")

    def clear_category(self, category: str) -> None:
        """Clear all cache files in a category."""
        cache_dir = self.cache_dir / category
        if cache_dir.exists():
            for cache_file in cache_dir.glob("*.pkl"):
                cache_file.unlink()
            self.logger.info(f"Cleared cache category: {category}")

    def clear_all(self) -> None:
        """Clear entire cache."""
        for cache_file in self.cache_dir.rglob("*.pkl"):
            cache_file.unlink()
        self.logger.info("Cleared all cache")


class NCBIClient:
    """Client for accessing NCBI databases."""

    def __init__(
        self, email: str, api_key: Optional[str] = None, request_delay: float = 0.34
    ):
        """
        Initialize NCBI client.

        Args:
            email: Email for NCBI (required)
            api_key: NCBI API key (optional, increases rate limit)
            request_delay: Delay between requests in seconds
        """
        if not is_valid_email(email):
            raise ValueError(ERROR_MESSAGES["invalid_email"])

        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        self.request_delay = request_delay
        self.last_request_time = 0
        self.logger = logging.getLogger(__name__)

    def _rate_limit(self) -> None:
        """Enforce rate limiting for NCBI requests."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time

        if time_since_last < self.request_delay:
            sleep_time = self.request_delay - time_since_last
            time.sleep(sleep_time)

        self.last_request_time = time.time()

    def search_gene(self, gene_name: str, organism: str = "homo sapiens") -> List[str]:
        """
        Search for gene IDs by name.

        Args:
            gene_name: Gene symbol or name
            organism: Target organism

        Returns:
            List of NCBI gene IDs
        """
        self._rate_limit()

        query = f"{gene_name}[Gene Name] AND {organism}[Organism]"

        try:
            handle = Entrez.esearch(db="gene", term=query, retmax=10)
            search_results = Entrez.read(handle)
            handle.close()

            gene_ids = search_results.get("IdList", [])
            self.logger.debug(f"Found {len(gene_ids)} genes for '{gene_name}'")
            return gene_ids

        except Exception as e:
            self.logger.error(f"Gene search failed for '{gene_name}': {e}")
            raise DatabaseError(f"Failed to search gene '{gene_name}': {e}")

    def get_gene_info(self, gene_id: str) -> Dict[str, Any]:
        """
        Get detailed gene information by ID.

        Args:
            gene_id: NCBI gene ID

        Returns:
            Gene information dictionary
        """
        self._rate_limit()

        try:
            handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
            gene_records = Entrez.read(handle)
            handle.close()

            if not gene_records:
                raise DatabaseError(f"No gene found with ID: {gene_id}")

            gene_record = gene_records[0]
            self.logger.debug(f"Retrieved gene info for ID: {gene_id}")
            return gene_record

        except Exception as e:
            self.logger.error(f"Failed to get gene info for ID {gene_id}: {e}")
            raise DatabaseError(f"Failed to retrieve gene {gene_id}: {e}")

    def get_sequence(self, accession: str) -> SeqRecord:
        """
        Get nucleotide sequence by accession number.

        Args:
            accession: RefSeq or GenBank accession

        Returns:
            BioPython SeqRecord object
        """
        self._rate_limit()

        try:
            handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="fasta", retmode="text"
            )
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()

            self.logger.debug(
                f"Retrieved sequence for {accession}: {len(seq_record)} bp"
            )
            return seq_record

        except Exception as e:
            self.logger.error(f"Failed to get sequence for {accession}: {e}")
            raise DatabaseError(f"Sequence not found: {accession}")

    def search_variants(self, gene_name: str, variant_type: str = "snp") -> List[Dict]:
        """
        Search for genetic variants in a gene.

        Args:
            gene_name: Gene symbol
            variant_type: Type of variant (snp, indel, etc.)

        Returns:
            List of variant information dictionaries
        """
        self._rate_limit()

        query = f"{gene_name}[Gene Name] AND {variant_type}[Variant Type]"

        try:
            handle = Entrez.esearch(db="snp", term=query, retmax=100)
            search_results = Entrez.read(handle)
            handle.close()

            variant_ids = search_results.get("IdList", [])

            # Get detailed info for each variant (limited to first 10)
            variants = []
            for var_id in variant_ids[:10]:
                var_info = self._get_variant_info(var_id)
                if var_info:
                    variants.append(var_info)

            self.logger.debug(f"Found {len(variants)} variants for {gene_name}")
            return variants

        except Exception as e:
            self.logger.error(f"Variant search failed for {gene_name}: {e}")
            return []

    def _get_variant_info(self, variant_id: str) -> Optional[Dict]:
        """Get detailed variant information."""
        try:
            self._rate_limit()
            handle = Entrez.efetch(db="snp", id=variant_id, retmode="xml")
            variant_data = Entrez.read(handle)
            handle.close()

            # Parse variant data (simplified)
            if variant_data:
                return {
                    "id": variant_id,
                    "rsid": f"rs{variant_id}",
                    "data": variant_data[0] if variant_data else None,
                }
            return None

        except Exception as e:
            self.logger.warning(f"Failed to get variant info for {variant_id}: {e}")
            return None


class UniProtClient:
    """Client for accessing UniProt protein database."""

    def __init__(self):
        """Initialize UniProt client."""
        self.base_url = EXTERNAL_RESOURCES["uniprot_base_url"]
        self.logger = logging.getLogger(__name__)

    def search_protein(self, gene_name: str, organism: str = "9606") -> List[Dict]:
        """
        Search for protein information by gene name.

        Args:
            gene_name: Gene symbol
            organism: Organism taxonomy ID (9606 = human)

        Returns:
            List of protein information dictionaries
        """
        query_url = f"{self.base_url}uniprotkb/search"
        params = {
            "query": f"gene:{gene_name} AND organism_id:{organism}",
            "format": "json",
            "size": 10,
        }

        try:
            response = requests.get(query_url, params=params, timeout=30)
            response.raise_for_status()

            data = response.json()
            proteins = data.get("results", [])

            self.logger.debug(f"Found {len(proteins)} proteins for {gene_name}")
            return proteins

        except Exception as e:
            self.logger.error(f"UniProt search failed for {gene_name}: {e}")
            return []


class GeneDatabase:
    """Main database interface for Enhancement Engine."""

    def __init__(
        self, email: str, cache_enabled: bool = True, cache_dir: str = "data/cache"
    ):
        """
        Initialize gene database.

        Args:
            email: Email for NCBI access
            cache_enabled: Enable local caching
            cache_dir: Cache directory path
        """
        self.ncbi = NCBIClient(email)
        self.uniprot = UniProtClient()

        self.cache_enabled = cache_enabled
        if cache_enabled:
            self.cache = CacheManager(cache_dir)

        self.logger = logging.getLogger(__name__)

        # Load enhancement genes from constants
        self._load_enhancement_genes()

    def _load_enhancement_genes(self) -> None:
        """Load enhancement gene database from constants."""
        self.enhancement_genes = {}

        for category_name, genes in ENHANCEMENT_GENES.items():
            category = EnhancementCategory(category_name)

            for gene_symbol, gene_data in genes.items():
                gene_info = GeneInfo(
                    name=gene_data["full_name"],
                    symbol=gene_symbol,
                    gene_id=gene_data["ncbi_id"],
                    chromosome=gene_data["chromosome"],
                    start_pos=0,  # Will be filled when sequence is loaded
                    end_pos=0,
                    description=gene_data["function"],
                    enhancement_category=category,
                    ncbi_id=gene_data["ncbi_id"],
                    refseq_id=gene_data["refseq_id"],
                )

                self.enhancement_genes[gene_symbol] = gene_info

        self.logger.info(f"Loaded {len(self.enhancement_genes)} enhancement genes")

    def search_gene(
        self, gene_name: str, organism: str = "homo sapiens"
    ) -> Optional[GeneInfo]:
        """
        Search for gene information.

        Args:
            gene_name: Gene symbol or name
            organism: Target organism

        Returns:
            GeneInfo object or None if not found
        """
        # Check if it's a known enhancement gene first
        if gene_name.upper() in self.enhancement_genes:
            return self.enhancement_genes[gene_name.upper()]

        # Check cache
        cache_key = f"{gene_name}_{organism}"
        if self.cache_enabled:
            cached_gene = self.cache.get(cache_key, "genes")
            if cached_gene:
                self.logger.debug(f"Retrieved {gene_name} from cache")
                return cached_gene

        try:
            # Search NCBI
            gene_ids = self.ncbi.search_gene(gene_name, organism)
            if not gene_ids:
                self.logger.warning(f"No genes found for '{gene_name}'")
                return None

            # Get detailed info for first result
            gene_info_raw = self.ncbi.get_gene_info(gene_ids[0])
            gene_info = self._parse_ncbi_gene(gene_info_raw)

            # Cache result
            if self.cache_enabled:
                self.cache.set(cache_key, gene_info, "genes")

            return gene_info

        except Exception as e:
            self.logger.error(f"Failed to search gene '{gene_name}': {e}")
            return None

    def _parse_ncbi_gene(self, gene_data: Dict) -> GeneInfo:
        """Parse NCBI gene data into GeneInfo object."""
        try:
            entrezgene = gene_data.get("Entrezgene", {})
            gene_ref = entrezgene.get("Entrezgene_gene", {}).get("Gene-ref", {})

            # Extract basic information
            symbol = gene_ref.get("Gene-ref_locus", "Unknown")
            description = gene_ref.get("Gene-ref_desc", "")

            # Extract location information
            locus_info = entrezgene.get("Entrezgene_locus", [])
            chromosome = "Unknown"
            start_pos = 0
            end_pos = 0

            if locus_info:
                location = locus_info[0].get("Gene-commentary_seqs", [])
                if location:
                    seq_loc = location[0].get("Seq-loc", {})
                    if "Seq-loc_int" in seq_loc:
                        interval = seq_loc["Seq-loc_int"]["Seq-interval"]
                        start_pos = interval.get("Seq-interval_from", 0)
                        end_pos = interval.get("Seq-interval_to", 0)

                        # Get chromosome from sequence ID
                        seq_id = interval.get("Seq-interval_id", {})
                        if "Seq-id_gi" in seq_id:
                            chromosome = str(seq_id["Seq-id_gi"])

            # Get aliases
            aliases = []
            synonyms = gene_ref.get("Gene-ref_syn", [])
            aliases.extend(synonyms)

            return GeneInfo(
                name=description,
                symbol=symbol,
                gene_id=str(
                    entrezgene.get("Entrezgene_track-info", {}).get(
                        "Gene-track_geneid", ""
                    )
                ),
                chromosome=chromosome,
                start_pos=start_pos,
                end_pos=end_pos,
                description=description,
                aliases=aliases,
            )

        except Exception as e:
            self.logger.error(f"Failed to parse NCBI gene data: {e}")
            raise DatabaseError(f"Failed to parse gene data: {e}")

    def get_sequence(self, accession: str) -> Optional[SeqRecord]:
        """
        Get nucleotide sequence by accession.

        Args:
            accession: RefSeq or GenBank accession

        Returns:
            BioPython SeqRecord or None if not found
        """
        # Check cache
        if self.cache_enabled:
            cached_seq = self.cache.get(accession, "sequences")
            if cached_seq:
                self.logger.debug(f"Retrieved sequence {accession} from cache")
                return cached_seq

        try:
            sequence = self.ncbi.get_sequence(accession)

            # Cache result
            if self.cache_enabled:
                self.cache.set(accession, sequence, "sequences")

            return sequence

        except Exception as e:
            self.logger.error(f"Failed to get sequence {accession}: {e}")
            return None

    def get_variants(self, gene_name: str) -> List[VariantInfo]:
        """
        Get known variants for a gene.

        Args:
            gene_name: Gene symbol

        Returns:
            List of VariantInfo objects
        """
        # Check cache
        cache_key = f"{gene_name}_variants"
        if self.cache_enabled:
            cached_variants = self.cache.get(cache_key, "variants")
            if cached_variants:
                self.logger.debug(f"Retrieved variants for {gene_name} from cache")
                return cached_variants

        try:
            variant_data = self.ncbi.search_variants(gene_name)
            variants = []

            for var_data in variant_data:
                variant = VariantInfo(
                    name=f"{gene_name}_variant",
                    rsid=var_data.get("rsid"),
                    variant_type="SNP",  # Simplified
                )
                variants.append(variant)

            # Cache result
            if self.cache_enabled:
                self.cache.set(cache_key, variants, "variants")

            return variants

        except Exception as e:
            self.logger.error(f"Failed to get variants for {gene_name}: {e}")
            return []

    def get_enhancement_genes(
        self, category: Optional[EnhancementCategory] = None
    ) -> Dict[str, GeneInfo]:
        """
        Get enhancement genes, optionally filtered by category.

        Args:
            category: Enhancement category to filter by

        Returns:
            Dictionary of gene symbol -> GeneInfo
        """
        if category is None:
            return self.enhancement_genes.copy()

        return {
            symbol: gene_info
            for symbol, gene_info in self.enhancement_genes.items()
            if gene_info.enhancement_category == category
        }

    def is_enhancement_gene(self, gene_name: str) -> bool:
        """Check if a gene is in the enhancement database."""
        return gene_name.upper() in self.enhancement_genes

    def get_gene_function(self, gene_name: str) -> Optional[str]:
        """Get functional description of a gene."""
        gene_info = self.search_gene(gene_name)
        return gene_info.description if gene_info else None

    def cache_data(self, key: str, data: Any, category: str = "general") -> None:
        """Manually cache data."""
        if self.cache_enabled:
            self.cache.set(key, data, category)

    def load_cached(self, key: str, category: str = "general") -> Optional[Any]:
        """Load data from cache."""
        if self.cache_enabled:
            return self.cache.get(key, category)
        return None

    def clear_cache(self, category: Optional[str] = None) -> None:
        """Clear cache for a category or all cache."""
        if self.cache_enabled:
            if category:
                self.cache.clear_category(category)
            else:
                self.cache.clear_all()

    def get_database_stats(self) -> Dict[str, Any]:
        """Get database statistics."""
        stats = {
            "enhancement_genes_total": len(self.enhancement_genes),
            "enhancement_categories": {},
            "cache_enabled": self.cache_enabled,
        }

        # Count genes by category
        for gene_info in self.enhancement_genes.values():
            category = gene_info.enhancement_category.value
            stats["enhancement_categories"][category] = (
                stats["enhancement_categories"].get(category, 0) + 1
            )

        # Cache statistics
        if self.cache_enabled:
            try:
                cache_files = list(self.cache.cache_dir.rglob("*.pkl"))
                stats["cache_files"] = len(cache_files)
                stats["cache_size_mb"] = sum(f.stat().st_size for f in cache_files) / (
                    1024 * 1024
                )
            except Exception:
                stats["cache_files"] = 0
                stats["cache_size_mb"] = 0

        return stats


# Convenience functions for common operations
def get_gene_info(gene_name: str, email: str) -> Optional[GeneInfo]:
    """Quick function to get gene information."""
    db = GeneDatabase(email)
    return db.search_gene(gene_name)


def get_enhancement_genes(
    category: Optional[str] = None, email: str = "temp@example.com"
) -> Dict[str, GeneInfo]:
    """Quick function to get enhancement genes."""
    db = GeneDatabase(email)
    cat_enum = EnhancementCategory(category) if category else None
    return db.get_enhancement_genes(cat_enum)


def search_sequences(
    accessions: List[str], email: str
) -> Dict[str, Optional[SeqRecord]]:
    """Batch sequence retrieval."""
    db = GeneDatabase(email)
    results = {}

    for accession in accessions:
        results[accession] = db.get_sequence(accession)
        time.sleep(0.34)  # Rate limiting

    return results
