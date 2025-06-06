"""
Disease gene constants and therapeutic targets.

This module contains comprehensive databases of disease-associated genes,
their variants, and therapeutic intervention strategies.
"""

from typing import Dict, List, Any

# Disease-associated genes database
DISEASE_GENES = {
    "autoimmune": {
        "PTPN22": {
            "full_name": "Protein Tyrosine Phosphatase Non-Receptor Type 22",
            "chromosome": "1",
            "function": "T-cell activation regulation",
            "disease_associations": {
                "rheumatoid_arthritis": {
                    "odds_ratio": 2.5,
                    "population_frequency": 0.12,
                    "penetrance": 0.15
                },
                "type1_diabetes": {
                    "odds_ratio": 2.0,
                    "population_frequency": 0.12,
                    "penetrance": 0.08
                },
                "systemic_lupus_erythematosus": {
                    "odds_ratio": 1.8,
                    "population_frequency": 0.12,
                    "penetrance": 0.05
                }
            },
            "pathogenic_variants": {
                "R620W": {
                    "rs_id": "rs2476601",
                    "nucleotide_change": "C1858T",
                    "protein_change": "R620W",
                    "effect": "increased_phosphatase_activity",
                    "therapeutic_target": True,
                    "correction_strategy": "base_editing_C_to_T_reversion"
                }
            },
            "ncbi_id": "26191",
            "refseq_id": "NM_015967.7",
            "therapeutic_approaches": ["correction", "silencing"]
        },
        "HLA-DRB1": {
            "full_name": "Major Histocompatibility Complex Class II DR Beta 1",
            "chromosome": "6",
            "function": "Antigen presentation to T-cells",
            "disease_associations": {
                "rheumatoid_arthritis": {
                    "odds_ratio": 5.8,
                    "population_frequency": 0.25,
                    "penetrance": 0.12
                },
                "multiple_sclerosis": {
                    "odds_ratio": 3.2,
                    "population_frequency": 0.15,
                    "penetrance": 0.08
                }
            },
            "pathogenic_variants": {
                "DRB1*04:01": {
                    "rs_id": "rs6457617",
                    "hla_allele": "DRB1*04:01:01",
                    "shared_epitope": "QKRAA",
                    "positions": "70-74",
                    "effect": "altered_antigen_presentation",
                    "therapeutic_target": True,
                    "correction_strategy": "allele_replacement"
                },
                "DRB1*04:04": {
                    "hla_allele": "DRB1*04:04:01",
                    "shared_epitope": "QKRAA",
                    "positions": "70-74",
                    "effect": "altered_antigen_presentation",
                    "therapeutic_target": True,
                    "correction_strategy": "allele_replacement"
                },
                "DRB1*01:01": {
                    "hla_allele": "DRB1*01:01:01",
                    "shared_epitope": "QKRAA",
                    "positions": "70-74",
                    "effect": "altered_antigen_presentation",
                    "therapeutic_target": True,
                    "correction_strategy": "allele_replacement"
                }
            },
            "protective_variants": {
                "DRB1*13:01": {
                    "hla_allele": "DRB1*13:01:01",
                    "effect": "protective_antigen_presentation",
                    "replacement_target": True
                },
                "DRB1*07:01": {
                    "hla_allele": "DRB1*07:01:01",
                    "effect": "neutral_to_protective",
                    "replacement_target": True
                }
            },
            "ncbi_id": "3123",
            "refseq_id": "NM_002124.4",
            "therapeutic_approaches": ["replacement", "targeted_mutagenesis"]
        },
        "STAT4": {
            "full_name": "Signal Transducer And Activator Of Transcription 4",
            "chromosome": "2",
            "function": "Cytokine signaling and T-helper cell differentiation",
            "disease_associations": {
                "rheumatoid_arthritis": {
                    "odds_ratio": 1.6,
                    "population_frequency": 0.28,
                    "penetrance": 0.08
                },
                "systemic_lupus_erythematosus": {
                    "odds_ratio": 1.8,
                    "population_frequency": 0.28,
                    "penetrance": 0.06
                }
            },
            "pathogenic_variants": {
                "rs7574865": {
                    "rs_id": "rs7574865",
                    "nucleotide_change": "G>T",
                    "effect": "increased_stat4_expression",
                    "therapeutic_target": True,
                    "correction_strategy": "expression_silencing"
                }
            },
            "ncbi_id": "6775",
            "refseq_id": "NM_003151.4",
            "therapeutic_approaches": ["silencing", "correction"]
        },
        "TNFAIP3": {
            "full_name": "TNF Alpha Induced Protein 3",
            "chromosome": "6",
            "function": "NF-κB pathway regulation",
            "disease_associations": {
                "rheumatoid_arthritis": {
                    "odds_ratio": 1.9,
                    "population_frequency": 0.15,
                    "penetrance": 0.07
                }
            },
            "pathogenic_variants": {
                "rs5029937": {
                    "rs_id": "rs5029937",
                    "nucleotide_change": "G>T",
                    "effect": "reduced_expression",
                    "therapeutic_target": True,
                    "correction_strategy": "activation"
                }
            },
            "ncbi_id": "7128",
            "refseq_id": "NM_006290.4",
            "therapeutic_approaches": ["activation", "correction"]
        }
    },
    "cardiovascular": {
        "LDLR": {
            "full_name": "Low Density Lipoprotein Receptor",
            "chromosome": "19",
            "function": "Cholesterol uptake",
            "disease_associations": {
                "familial_hypercholesterolemia": {
                    "odds_ratio": 50.0,
                    "population_frequency": 0.003,
                    "penetrance": 0.95
                }
            },
            "pathogenic_variants": {
                "c.682G>A": {
                    "nucleotide_change": "c.682G>A",
                    "protein_change": "E228K",
                    "effect": "reduced_ldl_binding",
                    "therapeutic_target": True,
                    "correction_strategy": "base_editing"
                }
            },
            "ncbi_id": "3949",
            "refseq_id": "NM_000527.5",
            "therapeutic_approaches": ["correction", "replacement"]
        },
        "PCSK9": {
            "full_name": "Proprotein Convertase Subtilisin/Kexin Type 9",
            "chromosome": "1",
            "function": "LDLR degradation",
            "disease_associations": {
                "hypercholesterolemia": {
                    "odds_ratio": 3.5,
                    "population_frequency": 0.02,
                    "penetrance": 0.7
                }
            },
            "pathogenic_variants": {
                "D374Y": {
                    "protein_change": "D374Y",
                    "effect": "gain_of_function",
                    "therapeutic_target": True,
                    "correction_strategy": "knockout"
                }
            },
            "protective_variants": {
                "R46L": {
                    "protein_change": "R46L",
                    "effect": "loss_of_function",
                    "cardioprotective": True
                }
            },
            "ncbi_id": "255738",
            "refseq_id": "NM_174936.4",
            "therapeutic_approaches": ["knockout", "silencing"]
        }
    }
}

# Therapeutic delivery methods for different tissues
TISSUE_DELIVERY_METHODS = {
    "synovial_tissue": {
        "primary": "intra_articular",
        "alternatives": ["local_injection", "systemic_targeted"],
        "efficiency": 0.7,
        "duration_days": 60
    },
    "hematopoietic_system": {
        "primary": "ex_vivo_editing",
        "alternatives": ["intravenous", "bone_marrow_injection"],
        "efficiency": 0.9,
        "duration_days": 365
    },
    "liver": {
        "primary": "intravenous",
        "alternatives": ["portal_vein", "hepatic_artery"],
        "efficiency": 0.8,
        "duration_days": 180
    },
    "muscle": {
        "primary": "intramuscular",
        "alternatives": ["systemic", "regional_perfusion"],
        "efficiency": 0.6,
        "duration_days": 90
    }
}

# Base editing windows and constraints
BASE_EDITING_CONSTRAINTS = {
    "cytosine_base_editors": {
        "BE3": {
            "editing_window": (4, 8),
            "target_base": "C",
            "product_base": "T",
            "efficiency_range": (0.2, 0.8),
            "bystander_editing": True
        },
        "BE4max": {
            "editing_window": (4, 8),
            "target_base": "C",
            "product_base": "T",
            "efficiency_range": (0.4, 0.9),
            "bystander_editing": True
        },
        "AID_BE3": {
            "editing_window": (1, 20),
            "target_base": "C",
            "product_base": "T",
            "efficiency_range": (0.1, 0.6),
            "bystander_editing": True
        }
    },
    "adenine_base_editors": {
        "ABE7.10": {
            "editing_window": (4, 8),
            "target_base": "A",
            "product_base": "G",
            "efficiency_range": (0.3, 0.9),
            "bystander_editing": True
        },
        "ABE8e": {
            "editing_window": (4, 8),
            "target_base": "A",
            "product_base": "G",
            "efficiency_range": (0.4, 0.95),
            "bystander_editing": False
        }
    }
}

# Prime editing constraints
PRIME_EDITING_CONSTRAINTS = {
    "PE3": {
        "max_insertion": 10,
        "max_deletion": 80,
        "max_replacement": 20,
        "efficiency_range": (0.1, 0.5),
        "indel_frequency": 0.05
    },
    "PE3-SpRY": {
        "max_insertion": 15,
        "max_deletion": 100,
        "max_replacement": 30,
        "efficiency_range": (0.15, 0.6),
        "indel_frequency": 0.03
    }
}

# Clinical disease scoring systems
DISEASE_SCORING_SYSTEMS = {
    "rheumatoid_arthritis": {
        "DAS28": {
            "full_name": "Disease Activity Score 28",
            "range": (0, 10),
            "remission_threshold": 2.6,
            "low_activity": 3.2,
            "high_activity": 5.1,
            "components": ["tender_joints", "swollen_joints", "ESR", "patient_global"]
        },
        "CDAI": {
            "full_name": "Clinical Disease Activity Index",
            "range": (0, 76),
            "remission_threshold": 2.8,
            "low_activity": 10,
            "high_activity": 22,
            "components": ["tender_joints", "swollen_joints", "patient_global", "physician_global"]
        }
    },
    "systemic_lupus_erythematosus": {
        "SLEDAI": {
            "full_name": "SLE Disease Activity Index",
            "range": (0, 105),
            "active_disease": 6,
            "components": ["seizure", "psychosis", "visual_disturbance", "cranial_nerve"]
        }
    }
}

# Population-specific allele frequencies
POPULATION_FREQUENCIES = {
    "european": {
        "PTPN22_R620W": 0.12,
        "HLA-DRB1*04:01": 0.15,
        "HLA-DRB1*04:04": 0.05,
        "STAT4_rs7574865": 0.28
    },
    "east_asian": {
        "PTPN22_R620W": 0.001,  # Very rare
        "HLA-DRB1*04:01": 0.02,
        "HLA-DRB1*04:04": 0.01,
        "STAT4_rs7574865": 0.15
    },
    "african": {
        "PTPN22_R620W": 0.002,
        "HLA-DRB1*04:01": 0.03,
        "HLA-DRB1*04:04": 0.02,
        "STAT4_rs7574865": 0.20
    },
    "hispanic": {
        "PTPN22_R620W": 0.08,
        "HLA-DRB1*04:01": 0.12,
        "HLA-DRB1*04:04": 0.04,
        "STAT4_rs7574865": 0.25
    }
}

# Therapeutic safety parameters
THERAPEUTIC_SAFETY_THRESHOLDS = {
    "off_target_tolerance": {
        "essential_genes": 0,      # No off-targets in essential genes
        "tumor_suppressors": 0,    # No off-targets in tumor suppressors
        "general": 3,              # Max 3 off-targets with >2 mismatches
        "high_confidence": 1       # Max 1 high-confidence off-target
    },
    "editing_efficiency": {
        "minimum_therapeutic": 0.1,  # 10% correction minimum
        "optimal_range": (0.3, 0.8), # 30-80% optimal
        "maximum_safe": 0.95         # Don't exceed 95%
    },
    "immune_response": {
        "anti_cas_antibodies": 0.3,   # 30% threshold for concern
        "t_cell_response": 0.2,       # 20% threshold
        "inflammatory_markers": 2.0    # 2x baseline
    }
}

# Quality control specifications
QUALITY_CONTROL_SPECS = {
    "sequencing_requirements": {
        "coverage": 1000,           # Minimum read depth
        "variant_frequency": 0.01,  # Detect variants at 1% frequency
        "off_target_sites": 50,     # Check top 50 off-target sites
        "time_points": [7, 30, 90, 365]  # Days post-treatment
    },
    "functional_assays": {
        "protein_expression": True,
        "enzymatic_activity": True,
        "cellular_phenotype": True,
        "immune_function": True
    }
}

# Regulatory pathway requirements
REGULATORY_REQUIREMENTS = {
    "preclinical": {
        "in_vitro_validation": ["specificity", "efficiency", "toxicity"],
        "animal_studies": ["safety", "efficacy", "biodistribution", "persistence"],
        "manufacturing": ["GMP_production", "quality_control", "stability"]
    },
    "clinical": {
        "phase1": ["safety", "dose_escalation", "preliminary_efficacy"],
        "phase2": ["efficacy", "optimal_dose", "patient_selection"],
        "phase3": ["comparative_efficacy", "long_term_safety", "quality_of_life"]
    }
}

# Error messages for therapeutic applications
THERAPEUTIC_ERROR_MESSAGES = {
    "invalid_disease": "Disease '{}' not found in therapeutic database",
    "no_therapeutic_target": "No therapeutic targets identified for gene '{}'",
    "uncorrectable_variant": "Variant '{}' cannot be corrected with available methods",
    "delivery_not_possible": "No suitable delivery method for target tissue '{}'",
    "safety_threshold_exceeded": "Safety thresholds exceeded - intervention not recommended"
}