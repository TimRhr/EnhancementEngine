"""
Constants and configuration parameters for Enhancement Engine.
"""

from typing import Dict, List, Any

# Enhancement-relevant genes database
ENHANCEMENT_GENES = {
    "cognitive": {
        "COMT": {
            "full_name": "Catechol-O-Methyltransferase",
            "chromosome": "22",
            "function": "Dopamine metabolism",
            "enhancement_variant": "Val158Met",
            "enhancement_effect": "Improved working memory and stress performance",
            "ncbi_id": "1312",
            "refseq_id": "NM_000754.3",
        },
        "BDNF": {
            "full_name": "Brain-Derived Neurotrophic Factor",
            "chromosome": "11",
            "function": "Neuroplasticity and neuronal survival",
            "enhancement_variant": "Val66Met",
            "enhancement_effect": "Enhanced learning and memory consolidation",
            "ncbi_id": "627",
            "refseq_id": "NM_170735.6",
        },
        "CACNA1C": {
            "full_name": "Calcium Voltage-Gated Channel Subunit Alpha1 C",
            "chromosome": "12",
            "function": "Calcium signaling in neurons",
            "enhancement_variant": "rs1006737",
            "enhancement_effect": "Improved memory formation",
            "ncbi_id": "775",
            "refseq_id": "NM_000719.8",
        },
    },
    "physical": {
        "ACTN3": {
            "full_name": "Actinin Alpha 3",
            "chromosome": "11",
            "function": "Fast-twitch muscle fiber composition",
            "enhancement_variant": "R577X",
            "enhancement_effect": "Enhanced power and sprint performance",
            "ncbi_id": "89",
            "refseq_id": "NM_001104.4",
        },
        "MSTN": {
            "full_name": "Myostatin",
            "chromosome": "2",
            "function": "Muscle growth inhibition",
            "enhancement_variant": "loss_of_function",
            "enhancement_effect": "Increased muscle mass and strength",
            "ncbi_id": "2660",
            "refseq_id": "NM_005259.4",
        },
        "EPO": {
            "full_name": "Erythropoietin",
            "chromosome": "7",
            "function": "Red blood cell production",
            "enhancement_variant": "overexpression",
            "enhancement_effect": "Enhanced oxygen delivery and endurance",
            "ncbi_id": "2056",
            "refseq_id": "NM_000799.3",
        },
    },
    "longevity": {
        "FOXO3": {
            "full_name": "Forkhead Box O3",
            "chromosome": "6",
            "function": "Cellular stress response and longevity",
            "enhancement_variant": "rs2802292",
            "enhancement_effect": "Extended lifespan and healthspan",
            "ncbi_id": "2309",
            "refseq_id": "NM_001455.3",
        },
        "APOE": {
            "full_name": "Apolipoprotein E",
            "chromosome": "19",
            "function": "Lipid metabolism and neuroprotection",
            "enhancement_variant": "E2_variant",
            "enhancement_effect": "Reduced Alzheimer's risk and cardiovascular protection",
            "ncbi_id": "348",
            "refseq_id": "NM_000041.3",
        },
    },
}

# CRISPR-Cas system parameters
PAM_PATTERNS = {
    "cas9": {
        "pattern": "NGG",
        "position": "3prime",
        "length": 3,
        "description": "SpCas9 PAM sequence",
    },
    "cas12a": {
        "pattern": "TTTV",
        "position": "5prime",
        "length": 4,
        "description": "Cas12a (Cpf1) PAM sequence",
    },
    "cas13": {
        "pattern": None,  # RNA-targeting, no PAM required
        "position": None,
        "length": 0,
        "description": "Cas13 for RNA targeting",
    },
}

CAS_TYPES = {
    "cas9": {
        "full_name": "CRISPR-associated protein 9",
        "target": "DNA",
        "cut_type": "blunt",
        "guide_length": 20,
        "pam": "NGG",
        "cut_position": -3,  # 3 bp upstream of PAM
    },
    "cas12a": {
        "full_name": "CRISPR-associated protein 12a",
        "target": "DNA",
        "cut_type": "staggered",
        "guide_length": 20,
        "pam": "TTTV",
        "cut_position": 18,  # Variable distance from PAM
    },
    "base_editor": {
        "full_name": "Base Editor (BE3, ABE)",
        "target": "DNA",
        "cut_type": "none",
        "guide_length": 20,
        "pam": "NGG",
        "edit_window": (4, 8),  # Editing window positions
    },
}

# NCBI database identifiers
NCBI_DATABASES = {
    "gene": "gene",
    "nucleotide": "nucleotide",
    "protein": "protein",
    "pubmed": "pubmed",
    "snp": "snp",
    "clinvar": "clinvar",
}

# Default analysis parameters
DEFAULT_PARAMETERS = {
    "guide_design": {
        "max_off_targets": 15,
        "min_efficiency_score": 0.3,
        "gc_content_range": (20, 80),
        "avoid_repeats": True,
        "max_mismatches": 4,
    },
    "safety_analysis": {
        "essential_genes_weight": 2.0,
        "off_target_threshold": 0.1,
        "immunogenicity_check": True,
        "genotoxicity_assessment": True,
    },
    "effect_simulation": {
        "confidence_interval": 0.95,
        "monte_carlo_iterations": 1000,
        "population_model": "mixed",
        "time_horizon_years": 10,
    },
}

# Scoring thresholds
SCORING_THRESHOLDS = {
    "guide_efficiency": {"excellent": 0.8, "good": 0.6, "fair": 0.4, "poor": 0.2},
    "safety_score": {
        "very_safe": 90,
        "safe": 70,
        "moderate_risk": 50,
        "high_risk": 30,
        "very_high_risk": 10,
    },
    "enhancement_effect": {
        "major": 2.0,  # >2x improvement
        "moderate": 1.5,  # 1.5-2x improvement
        "minor": 1.2,  # 1.2-1.5x improvement
        "negligible": 1.1,  # <1.2x improvement
    },
}

# File paths and URLs
EXTERNAL_RESOURCES = {
    "ncbi_base_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    "uniprot_base_url": "https://rest.uniprot.org/",
    "ensembl_base_url": "https://rest.ensembl.org/",
    "alphafold_base_url": "https://alphafold.ebi.ac.uk/api/",
}

# Error messages
ERROR_MESSAGES = {
    "invalid_gene": "Gene '{}' not found in database",
    "invalid_variant": "Variant '{}' not recognized for gene '{}'",
    "ncbi_connection": "Failed to connect to NCBI databases",
    "sequence_not_found": "Sequence not found for accession '{}'",
    "invalid_email": "Valid email required for NCBI access",
}

# Enhancement categories for reporting
ENHANCEMENT_CATEGORIES = [
    "cognitive",
    "physical",
    "longevity",
    "sensory",
    "metabolic",
    "immune",
]
