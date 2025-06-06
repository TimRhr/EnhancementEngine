"""
Disease gene constants and therapeutic targets.

This module contains comprehensive databases of disease-associated genes,
their variants, and therapeutic intervention strategies.
"""

from typing import Dict, List, Any, Tuple

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

# Reference HLA allele information used by the HLA analyzer
HLA_ALLELE_DATABASE = {
    "DRB1*01:01": {
        "sequence": "SEQ_01:01",
        "shared_epitope": "QKRAA",
        "risk_category": "high",
        "population_frequencies": {"European": 0.08, "Asian": 0.01, "African": 0.02},
    },
    "DRB1*04:01": {
        "sequence": "SEQ_04:01",
        "shared_epitope": "QKRAA",
        "risk_category": "very_high",
        "population_frequencies": {"European": 0.12, "Asian": 0.03, "African": 0.01},
    },
    "DRB1*04:04": {
        "sequence": "SEQ_04:04",
        "shared_epitope": "QKRAA",
        "risk_category": "high",
        "population_frequencies": {"European": 0.03, "Asian": 0.15, "African": 0.01},
    },
    "DRB1*04:05": {
        "sequence": "SEQ_04:05",
        "shared_epitope": "QRRAA",
        "risk_category": "moderate",
        "population_frequencies": {"European": 0.05, "Asian": 0.02, "African": 0.08},
    },
    "DRB1*03:01": {
        "sequence": "SEQ_03:01",
        "shared_epitope": None,
        "risk_category": "protective",
        "population_frequencies": {"European": 0.15, "Asian": 0.25, "African": 0.20},
    },
    "DRB1*13:01": {
        "sequence": "SEQ_13:01",
        "shared_epitope": None,
        "risk_category": "neutral",
        "population_frequencies": {"European": 0.12, "Asian": 0.08, "African": 0.15},
    },
    "DRB1*15:01": {
        "sequence": "SEQ_15:01",
        "shared_epitope": None,
        "risk_category": "protective",
        "population_frequencies": {"European": 0.18, "Asian": 0.05, "African": 0.25},
    },
}

# Shared epitope definitions
SHARED_EPITOPES = {
    "SE1": {
        "sequence": "QKRAA",
        "positions": (70, 74),
        "risk_weight": 1.0,
        "alleles": ["DRB1*01:01", "DRB1*04:01", "DRB1*04:04"],
    },
    "SE2": {
        "sequence": "QRRAA",
        "positions": (70, 74),
        "risk_weight": 0.6,
        "alleles": ["DRB1*04:05"],
    },
    "SE3": {
        "sequence": "RRRAA",
        "positions": (70, 74),
        "risk_weight": 0.4,
        "alleles": ["DRB1*04:08"],
    },
}

# Gene interaction network used by the combination therapy module
GENE_INTERACTIONS = {
    "PTPN22": {
        "HLA-DRB1": {"type": "epistatic", "strength": 0.7},
        "STAT4": {"type": "additive", "strength": 0.4},
        "TNFAIP3": {"type": "compensatory", "strength": 0.3},
    },
    "HLA-DRB1": {
        "PTPN22": {"type": "epistatic", "strength": 0.7},
        "STAT4": {"type": "modifying", "strength": 0.5},
        "PADI4": {"type": "synergistic", "strength": 0.6},
    },
    "STAT4": {
        "PTPN22": {"type": "additive", "strength": 0.4},
        "HLA-DRB1": {"type": "modifying", "strength": 0.5},
        "IRF4": {"type": "pathway", "strength": 0.8},
    },
}

# Synergy data for combination therapy
THERAPEUTIC_SYNERGIES = {
    "pairwise": {
        ("PTPN22", "HLA-DRB1"): {
            "synergy_type": "multiplicative",
            "synergy_factor": 1.8,
            "confidence": 0.8,
            "mechanism": "Combined T-cell regulation and antigen presentation",
        },
        ("HLA-DRB1", "STAT4"): {
            "synergy_type": "additive_plus",
            "synergy_factor": 1.3,
            "confidence": 0.6,
            "mechanism": "Antigen presentation and T-cell differentiation",
        },
        ("PTPN22", "STAT4"): {
            "synergy_type": "additive",
            "synergy_factor": 1.2,
            "confidence": 0.7,
            "mechanism": "Complementary T-cell regulation pathways",
        },
    },
    "multi_gene": {
        ("PTPN22", "HLA-DRB1", "STAT4"): {
            "synergy_factor": 2.5,
            "complexity_penalty": 0.3,
            "coordination_requirement": "high",
        }
    },
}

# High level therapeutic target summary used by delivery system
THERAPEUTIC_TARGETS = {
    "PTPN22": {"strategies": ["correction", "silencing"]},
    "HLA-DRB1": {"strategies": ["replacement", "targeted_mutagenesis"]},
    "STAT4": {"strategies": ["silencing", "correction"]},
    "TNFAIP3": {"strategies": ["activation", "correction"]},
    "LDLR": {"strategies": ["correction", "replacement"]},
    "PCSK9": {"strategies": ["knockout", "silencing"]},
}

# Default validation criteria for therapeutic strategies
VALIDATION_CRITERIA = {
    "base_editing": {
        "min_efficiency": 0.6,
        "max_off_targets": 3,
        "safety_threshold": 70.0,
        "efficacy_threshold": 0.5,
        "population_benefit_threshold": 0.3,
    },
    "gene_replacement": {
        "min_efficiency": 0.4,
        "max_off_targets": 5,
        "safety_threshold": 80.0,
        "efficacy_threshold": 0.7,
        "population_benefit_threshold": 0.5,
    },
    "gene_silencing": {
        "min_efficiency": 0.7,
        "max_off_targets": 2,
        "safety_threshold": 75.0,
        "efficacy_threshold": 0.4,
        "population_benefit_threshold": 0.2,
    },
}

# Population genetics data and disease epidemiology
POPULATION_GENETICS = {
    "demographics": {
        "total_population": 330000000,
        "age_distribution": {"0-17": 0.22, "18-34": 0.23, "35-54": 0.25, "55-74": 0.20, "75+": 0.10},
        "sex_distribution": {"male": 0.49, "female": 0.51},
        "ethnicity_distribution": {
            "caucasian": 0.60,
            "hispanic": 0.18,
            "african_american": 0.13,
            "asian": 0.06,
            "other": 0.03,
        },
    },
    "genetic_data": {
        "rheumatoid_arthritis": {
            "overall_prevalence": 0.01,
            "genetic_risk_distribution": {"low": 0.65, "moderate": 0.25, "high": 0.08, "very_high": 0.02},
            "heritability": 0.60,
            "environmental_factors": {
                "smoking": {"prevalence": 0.15, "risk_increase": 2.0},
                "infections": {"prevalence": 0.30, "risk_increase": 1.5},
            },
        }
    },
}

DISEASE_EPIDEMIOLOGY = {
    "rheumatoid_arthritis": {
        "incidence_rate_per_100k": 50,
        "age_incidence_curve": {
            "20-29": 0.1,
            "30-39": 0.15,
            "40-49": 0.25,
            "50-59": 0.30,
            "60-69": 0.20,
        },
        "progression_model": {
            "remission_rate": 0.20,
            "mild_to_moderate": 0.40,
            "moderate_to_severe": 0.30,
            "mortality_increase": 1.3,
        },
        "economic_burden": {
            "annual_direct_cost": 15000,
            "annual_indirect_cost": 12000,
            "qaly_loss_per_year": 0.3,
        },
    },
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