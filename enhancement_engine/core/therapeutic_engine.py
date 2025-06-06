"""
Therapeutic Enhancement Engine - Main interface for disease correction.

This module provides the primary interface for therapeutic genetic interventions,
integrating disease risk assessment, therapeutic design, and safety evaluation.
"""

import logging
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union, Any
from pathlib import Path

from .engine import EnhancementEngine
from .disease_risk import DiseaseRiskCalculator, MultiGeneRiskCalculator
from .disease_db import DiseaseDatabaseClient
from .therapeutic_crispr import TherapeuticCRISPRDesigner
from .therapeutic_safety import TherapeuticSafetyAnalyzer
from .database import GeneDatabase
from .sequence import SequenceAnalyzer

from ..models.data_classes import GeneInfo, GuideRNA, EnhancementReport
from ..models.therapeutic_data_classes import (
    TherapeuticReport,
    TherapeuticTarget,
    CorrectionStrategy,
    TherapeuticEfficacy,
    TherapeuticSafety,
    PatientStratification,
    DiseaseRisk,
    TherapeuticApproach,
    DeliveryMethod,
    DiseaseCategory,
    ClinicalEndpoint,
    ReversalStrategy,
)
from ..models.disease_constants import (
    DISEASE_GENES,
    THERAPEUTIC_SAFETY_THRESHOLDS,
    DISEASE_SCORING_SYSTEMS,
)
from ..utils import save_json, current_timestamp, is_valid_email


class TherapeuticEngineError(Exception):
    """Custom exception for Therapeutic Engine errors."""

    pass


class TherapeuticEnhancementEngine:
    """
    Main Therapeutic Enhancement Engine for disease correction.

    This class orchestrates all therapeutic analysis components to provide:
    - Disease risk assessment and genetic counseling
    - Therapeutic target identification and prioritization
    - CRISPR-based correction strategy design
    - Comprehensive safety evaluation
    - Clinical efficacy prediction
    - Patient stratification and personalized recommendations
    """

    def __init__(self, email: str, config: Optional[Dict[str, Any]] = None):
        """
        Initialize Therapeutic Enhancement Engine.

        Args:
            email: Email address for NCBI access (required)
            config: Optional configuration parameters
        """
        if not is_valid_email(email):
            raise TherapeuticEngineError("Valid email required for NCBI access")

        self.email = email
        self.config = config or {}

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # Initialize core components
        self._initialize_therapeutic_components()

        # Analysis cache for performance
        self._therapeutic_cache = {}

        self.logger.info(f"Therapeutic Enhancement Engine initialized for {email}")

    def _initialize_therapeutic_components(self):
        """Initialize all therapeutic analysis components."""
        try:
            # Core enhancement engine for basic functionality
            self.base_engine = EnhancementEngine(self.email)

            # Therapeutic-specific components
            self.disease_db_client = DiseaseDatabaseClient(self.email)
            self.disease_risk_calculator = DiseaseRiskCalculator(self.disease_db_client)
            self.multigene_risk_calculator = MultiGeneRiskCalculator()
            self.therapeutic_crispr_designer = TherapeuticCRISPRDesigner()
            self.therapeutic_safety_analyzer = TherapeuticSafetyAnalyzer()

            # Clinical assessment components
            self.patient_stratifier = PatientStratificationEngine()
            self.clinical_assessor = ClinicalAssessmentEngine()
            self.therapeutic_simulator = TherapeuticSimulationEngine()

            self.logger.info("All therapeutic components initialized successfully")

        except Exception as e:
            self.logger.error(f"Failed to initialize therapeutic components: {e}")
            raise TherapeuticEngineError(f"Component initialization failed: {e}")

    def analyze_disease_gene(
        self,
        gene_name: str,
        variant: str,
        patient_data: Dict[str, Any],
        disease: str = "rheumatoid_arthritis",
    ) -> TherapeuticReport:
        """
        Comprehensive analysis of disease gene for therapeutic intervention.

        Args:
            gene_name: Gene symbol (e.g., "PTPN22", "HLA-DRB1")
            variant: Disease-associated variant
            patient_data: Patient demographics and clinical data
            disease: Target disease for treatment

        Returns:
            TherapeuticReport with complete analysis
        """
        try:
            self.logger.info(
                f"Starting therapeutic analysis of {gene_name} {variant} for {disease}"
            )

            # Check cache first
            cache_key = (
                f"therapeutic_{gene_name}_{variant}_{disease}_{hash(str(patient_data))}"
            )
            if cache_key in self._therapeutic_cache:
                self.logger.debug(f"Retrieved {cache_key} from cache")
                return self._therapeutic_cache[cache_key]

            # Step 1: Disease risk assessment
            disease_risk = self._assess_disease_risk(
                gene_name, variant, patient_data, disease
            )

            # Step 2: Patient stratification
            patient_stratification = self._stratify_patient(patient_data, disease_risk)

            # Step 3: Identify therapeutic targets
            therapeutic_targets = self._identify_therapeutic_targets(
                gene_name, variant, disease, patient_stratification
            )

            # Step 4: Design correction strategies
            correction_strategies = self._design_correction_strategies(
                therapeutic_targets, gene_name
            )

            # Step 5: Safety assessment
            safety_assessment = self._assess_therapeutic_safety(
                therapeutic_targets, correction_strategies, patient_data
            )

            # Step 6: Efficacy prediction
            predicted_efficacy = self._predict_therapeutic_efficacy(
                correction_strategies, therapeutic_targets, patient_stratification
            )

            # Step 7: Clinical endpoints
            clinical_endpoints = self._define_clinical_endpoints(
                disease, patient_stratification, predicted_efficacy
            )

            # Step 8: Generate recommendations
            recommendations = self._generate_therapeutic_recommendations(
                disease_risk,
                safety_assessment,
                predicted_efficacy,
                patient_stratification,
            )

            # Step 9: Identify contraindications
            contraindications = self._identify_contraindications(
                safety_assessment, patient_data, disease_risk
            )

            # Step 10: Monitoring plan
            monitoring_plan = self._generate_monitoring_plan(
                safety_assessment, predicted_efficacy, therapeutic_targets
            )

            # Create comprehensive therapeutic report
            therapeutic_report = TherapeuticReport(
                patient_id=patient_data.get("patient_id", "anonymous"),
                target_disease=disease,
                therapeutic_targets=therapeutic_targets,
                correction_strategies=correction_strategies,
                predicted_efficacy=predicted_efficacy,
                safety_assessment=safety_assessment,
                patient_stratification=patient_stratification,
                clinical_endpoints=clinical_endpoints,
                recommendations=recommendations,
                contraindications=contraindications,
                monitoring_plan=monitoring_plan,
                follow_up_schedule=self._generate_follow_up_schedule(
                    predicted_efficacy
                ),
                confidence_score=self._calculate_therapeutic_confidence(
                    safety_assessment, predicted_efficacy, len(correction_strategies)
                ),
            )

            # Cache result
            self._therapeutic_cache[cache_key] = therapeutic_report

            self.logger.info(
                f"Therapeutic analysis completed for {gene_name}: "
                f"{therapeutic_report.treatment_recommendation}"
            )

            return therapeutic_report

        except Exception as e:
            self.logger.error(f"Therapeutic analysis failed for {gene_name}: {e}")
            raise TherapeuticEngineError(f"Therapeutic analysis failed: {e}")

    def _assess_disease_risk(
        self, gene_name: str, variant: str, patient_data: Dict[str, Any], disease: str
    ) -> DiseaseRisk:
        """Assess disease risk for genetic variant."""
        try:
            disease_risk = self.disease_risk_calculator.calculate_disease_risk(
                gene_name, variant, patient_data, disease
            )

            self.logger.debug(
                f"Disease risk calculated: OR={disease_risk.odds_ratio:.2f}"
            )
            return disease_risk

        except Exception as e:
            self.logger.error(f"Disease risk assessment failed: {e}")
            # Return default risk assessment
            return DiseaseRisk(
                disease_name=disease,
                odds_ratio=1.5,
                population_frequency=0.1,
                absolute_risk_increase=0.01,
                confidence_interval=(1.2, 1.8),
                penetrance=0.1,
            )

    def _stratify_patient(
        self, patient_data: Dict[str, Any], disease_risk: DiseaseRisk
    ) -> PatientStratification:
        """Stratify patient for treatment candidacy."""
        genotype_data = {
            gene: variant for gene, variant in patient_data.get("genotypes", {}).items()
        }

        # Determine disease severity
        disease_severity = patient_data.get("disease_severity", "moderate")

        # Calculate treatment priority based on risk and severity
        treatment_priority = self._calculate_treatment_priority(
            disease_risk, disease_severity
        )

        return PatientStratification(
            patient_genotype=genotype_data,
            disease_severity=disease_severity,
            prior_treatments=patient_data.get("prior_treatments", []),
            comorbidities=patient_data.get("comorbidities", []),
            ethnicity=patient_data.get("ethnicity", "european"),
            age=patient_data.get("age", 40),
            sex=patient_data.get("sex", "female"),
            treatment_priority=treatment_priority,
        )

    def _calculate_treatment_priority(
        self, disease_risk: DiseaseRisk, disease_severity: str
    ) -> float:
        """Calculate treatment priority score."""
        # Base priority from disease risk
        risk_priority = min(1.0, disease_risk.odds_ratio / 10.0)

        # Severity modifier
        severity_modifiers = {"mild": 0.3, "moderate": 0.6, "severe": 1.0}
        severity_priority = severity_modifiers.get(disease_severity, 0.6)

        # Combined priority
        return risk_priority * 0.4 + severity_priority * 0.6

    def _identify_therapeutic_targets(
        self,
        gene_name: str,
        variant: str,
        disease: str,
        patient_strat: PatientStratification,
    ) -> List[TherapeuticTarget]:
        """Identify and prioritize therapeutic targets."""
        targets = []

        # Get gene data from disease database
        gene_data = None
        for category in DISEASE_GENES.values():
            if gene_name in category:
                gene_data = category[gene_name]
                break

        if not gene_data:
            return targets

        # Determine therapeutic approach based on variant
        therapeutic_approach = self._determine_therapeutic_approach(
            gene_name, variant, gene_data
        )

        # Determine target tissue
        target_tissue = self._determine_target_tissue(disease, gene_name)

        # Determine delivery method
        delivery_method = self._determine_delivery_method(
            target_tissue, therapeutic_approach
        )

        # Get correction sequence if needed
        correction_sequence = self._get_correction_sequence(
            gene_name, variant, therapeutic_approach
        )

        # Create therapeutic target
        target = TherapeuticTarget(
            gene_symbol=gene_name,
            disease_association=disease,
            target_variant=variant,
            therapeutic_approach=therapeutic_approach,
            target_tissue=target_tissue,
            delivery_method=delivery_method,
            correction_sequence=correction_sequence,
            target_position=self._get_target_position(gene_name, variant),
            priority_score=patient_strat.treatment_priority,
        )

        targets.append(target)
        return targets

    def _determine_therapeutic_approach(
        self, gene_name: str, variant: str, gene_data: Dict[str, Any]
    ) -> TherapeuticApproach:
        """Determine optimal therapeutic approach."""
        # Check therapeutic approaches in gene data
        approaches = gene_data.get("therapeutic_approaches", [])

        # Gene-specific approach selection
        if gene_name == "PTPN22" and variant == "R620W":
            return TherapeuticApproach.CORRECTION  # Base editing to fix mutation
        elif gene_name.startswith("HLA-"):
            return TherapeuticApproach.REPLACEMENT  # Replace with protective allele
        elif gene_name == "PCSK9":
            return TherapeuticApproach.KNOCKOUT  # Knockout for cardioprotection
        elif "silencing" in approaches:
            return TherapeuticApproach.SILENCING
        elif "correction" in approaches:
            return TherapeuticApproach.CORRECTION
        else:
            return TherapeuticApproach.KNOCKOUT  # Default approach

    def _determine_target_tissue(self, disease: str, gene_name: str) -> str:
        """Determine target tissue for intervention."""
        tissue_mapping = {
            "rheumatoid_arthritis": "synovial_tissue",
            "systemic_lupus_erythematosus": "hematopoietic_system",
            "familial_hypercholesterolemia": "liver",
            "type1_diabetes": "pancreatic_islets",
        }

        # Special cases
        if gene_name.startswith("HLA-"):
            return "hematopoietic_system"  # HLA editing requires HSC modification

        return tissue_mapping.get(disease, "systemic")

    def _determine_delivery_method(
        self, target_tissue: str, therapeutic_approach: TherapeuticApproach
    ) -> DeliveryMethod:
        """Determine optimal delivery method."""
        delivery_mapping = {
            "synovial_tissue": DeliveryMethod.INTRA_ARTICULAR,
            "hematopoietic_system": DeliveryMethod.INTRAVENOUS,  # Ex vivo editing
            "liver": DeliveryMethod.INTRAVENOUS,
            "pancreatic_islets": DeliveryMethod.LOCAL_INJECTION,
            "muscle": DeliveryMethod.INTRAMUSCULAR,
        }

        return delivery_mapping.get(target_tissue, DeliveryMethod.INTRAVENOUS)

    def _get_correction_sequence(
        self, gene_name: str, variant: str, approach: TherapeuticApproach
    ) -> Optional[str]:
        """Get desired correction sequence."""
        if approach == TherapeuticApproach.CORRECTION:
            # For PTPN22 R620W: correct back to R620R
            if gene_name == "PTPN22" and variant == "R620W":
                return "R620R"  # Wild-type sequence
        elif approach == TherapeuticApproach.REPLACEMENT:
            # For HLA-DRB1: replace with protective allele
            if gene_name == "HLA-DRB1":
                return "DRB1*13:01"  # Protective allele

        return None

    def _get_target_position(self, gene_name: str, variant: str) -> int:
        """Get genomic position for targeting."""
        # Simplified position mapping
        position_mapping = {
            ("PTPN22", "R620W"): 1858,  # Position of C1858T mutation
            ("HLA-DRB1", "DRB1*04:01"): 200,  # Approximate exon 2 position
            ("STAT4", "rs7574865"): 500,  # Approximate position
        }

        return position_mapping.get((gene_name, variant), 1000)  # Default position

    def _design_correction_strategies(
        self, targets: List[TherapeuticTarget], gene_name: str
    ) -> List[CorrectionStrategy]:
        """Design correction strategies for therapeutic targets."""
        strategies = []

        for target in targets:
            try:
                # Get target sequence (simplified)
                target_sequence = self._get_target_sequence(gene_name)

                if target_sequence:
                    # Design strategies using therapeutic CRISPR designer
                    target_strategies = self.therapeutic_crispr_designer.design_therapeutic_intervention(
                        target, target_sequence
                    )
                    strategies.extend(target_strategies)

            except Exception as e:
                self.logger.warning(
                    f"Failed to design strategy for {target.gene_symbol}: {e}"
                )

        return strategies

    def _get_target_sequence(self, gene_name: str) -> Optional[str]:
        """Get target gene sequence."""
        try:
            # Use base engine to get gene sequence
            gene_info = self.base_engine._get_gene_information(gene_name)
            if gene_info and gene_info.refseq_id:
                sequence_record = self.base_engine.gene_database.get_sequence(
                    gene_info.refseq_id
                )
                if sequence_record:
                    return str(sequence_record.seq)
        except Exception as e:
            self.logger.warning(f"Failed to get sequence for {gene_name}: {e}")

        # Return mock sequence for demonstration
        return "ATCG" * 500  # 2000 bp mock sequence

    def _assess_therapeutic_safety(
        self,
        targets: List[TherapeuticTarget],
        strategies: List[CorrectionStrategy],
        patient_data: Dict[str, Any],
    ) -> TherapeuticSafety:
        """Assess therapeutic safety."""
        if not targets or not strategies:
            return TherapeuticSafety(
                off_target_risk=5.0,
                immunogenicity_risk=5.0,
                tissue_toxicity_risk=5.0,
                systemic_effects_risk=5.0,
                reversibility_score=0.3,
            )

        return self.therapeutic_safety_analyzer.analyze_therapeutic_safety(
            targets[0], strategies, patient_data
        )

    def _predict_therapeutic_efficacy(
        self,
        strategies: List[CorrectionStrategy],
        targets: List[TherapeuticTarget],
        patient_strat: PatientStratification,
    ) -> TherapeuticEfficacy:
        """Predict therapeutic efficacy."""
        if not strategies or not targets:
            return TherapeuticEfficacy(
                correction_efficiency=0.3,
                tissue_penetration=0.5,
                persistence_months=3.0,
                clinical_improvement=0.2,
                remission_probability=0.1,
            )

        return self.therapeutic_crispr_designer.predict_therapeutic_efficacy(
            strategies, targets[0].delivery_method, targets[0].target_tissue
        )

    def _define_clinical_endpoints(
        self,
        disease: str,
        patient_strat: PatientStratification,
        efficacy: TherapeuticEfficacy,
    ) -> List[ClinicalEndpoint]:
        """Define clinical endpoints for assessment."""
        endpoints = []

        if disease == "rheumatoid_arthritis":
            # DAS28 score endpoint
            baseline_das28 = 5.5 if patient_strat.disease_severity == "severe" else 4.0
            target_das28 = baseline_das28 * (1 - efficacy.clinical_improvement)

            endpoints.append(
                ClinicalEndpoint(
                    endpoint_name="DAS28",
                    baseline_value=baseline_das28,
                    target_value=target_das28,
                    measurement_unit="score",
                    time_to_assessment_weeks=12,
                    minimum_clinically_important_difference=1.2,
                )
            )

            # CDAI score endpoint
            baseline_cdai = 35.0 if patient_strat.disease_severity == "severe" else 25.0
            target_cdai = baseline_cdai * (1 - efficacy.clinical_improvement)

            endpoints.append(
                ClinicalEndpoint(
                    endpoint_name="CDAI",
                    baseline_value=baseline_cdai,
                    target_value=target_cdai,
                    measurement_unit="score",
                    time_to_assessment_weeks=12,
                    minimum_clinically_important_difference=10.0,
                )
            )

        return endpoints

    def _generate_therapeutic_recommendations(
        self,
        disease_risk: DiseaseRisk,
        safety: TherapeuticSafety,
        efficacy: TherapeuticEfficacy,
        patient_strat: PatientStratification,
    ) -> List[str]:
        """Generate therapeutic recommendations."""
        recommendations = []

        # Treatment recommendation based on benefit-risk
        overall_benefit = efficacy.overall_efficacy
        # Convert safety score (0-100) to risk scale (0-1)
        overall_risk = (100 - safety.overall_safety_score) / 100

        if overall_benefit > 0.7 and overall_risk < 0.3:
            recommendations.append("Strongly recommended for treatment")
        elif overall_benefit > 0.5 and overall_risk < 0.5:
            recommendations.append("Recommended with careful monitoring")
        elif overall_benefit > 0.3:
            recommendations.append(
                "Consider treatment if conventional therapies have failed"
            )
        else:
            recommendations.append("Not recommended - risks outweigh benefits")

        # Risk-specific recommendations
        if safety.off_target_risk > 3.0:
            recommendations.append("Extensive off-target monitoring required")

        if safety.immunogenicity_risk > 4.0:
            recommendations.append("Consider immunosuppressive premedication")

        if safety.tissue_toxicity_risk > 3.0:
            recommendations.append("Implement tissue-specific safety monitoring")

        # Efficacy-specific recommendations
        if efficacy.correction_efficiency < 0.3:
            recommendations.append(
                "Consider multiple treatment cycles for optimal effect"
            )

        if efficacy.persistence_months < 6:
            recommendations.append("Plan for repeat treatments to maintain benefit")

        # Patient-specific recommendations
        if patient_strat.treatment_priority > 0.8:
            recommendations.append(
                "High priority case - expedited treatment consideration"
            )

        if patient_strat.age > 65:
            recommendations.append("Enhanced monitoring due to advanced age")

        return recommendations

    def _identify_contraindications(
        self,
        safety: TherapeuticSafety,
        patient_data: Dict[str, Any],
        disease_risk: DiseaseRisk,
    ) -> List[str]:
        """Identify contraindications for treatment."""
        contraindications = []

        # Safety-based contraindications
        if safety.overall_safety_score < 30:
            contraindications.append("Prohibitively high safety risk")

        # Patient-based contraindications
        comorbidities = patient_data.get("comorbidities", [])

        if "pregnancy" in comorbidities:
            contraindications.append("Pregnancy - unknown fetal effects")

        if "active_malignancy" in comorbidities:
            contraindications.append("Active malignancy - genetic instability risk")

        if "severe_immunodeficiency" in comorbidities:
            contraindications.append("Severe immunodeficiency")

        # Age-based contraindications
        age = patient_data.get("age", 40)
        if age < 16:
            contraindications.append("Pediatric population - safety not established")

        return contraindications

    def _generate_monitoring_plan(
        self,
        safety: TherapeuticSafety,
        efficacy: TherapeuticEfficacy,
        targets: List[TherapeuticTarget],
    ) -> List[str]:
        """Generate monitoring plan."""
        monitoring = []

        # Standard monitoring
        monitoring.extend(
            [
                "Pre-treatment comprehensive medical evaluation",
                "Baseline disease activity assessment",
                "Genetic counseling and informed consent",
            ]
        )

        # Safety-based monitoring
        monitoring.extend(safety.monitoring_requirements)

        # Efficacy monitoring
        monitoring.extend(
            [
                "Weekly disease activity scores for first month",
                "Monthly clinical assessments for 6 months",
                "Quarterly long-term follow-up",
            ]
        )

        # Target-specific monitoring
        if targets and targets[0].target_tissue == "hematopoietic_system":
            monitoring.extend(
                [
                    "Complete blood count monitoring",
                    "Immune function assessment",
                    "Bone marrow evaluation if indicated",
                ]
            )

        return monitoring

    def _generate_follow_up_schedule(
        self, efficacy: TherapeuticEfficacy
    ) -> Dict[str, str]:
        """Generate follow-up schedule."""
        schedule = {
            "24_hours": "Safety assessment and vital signs",
            "1_week": "Initial safety and tolerability evaluation",
            "1_month": "Efficacy assessment and safety monitoring",
            "3_months": "Comprehensive efficacy and safety evaluation",
            "6_months": "Long-term efficacy and durability assessment",
            "1_year": "Annual comprehensive evaluation",
        }

        # Adjust based on efficacy
        if efficacy.persistence_months < 6:
            schedule["2_months"] = (
                "Early efficacy assessment for potential re-treatment"
            )

        return schedule

    def _calculate_therapeutic_confidence(
        self,
        safety: TherapeuticSafety,
        efficacy: TherapeuticEfficacy,
        num_strategies: int,
    ) -> float:
        """Calculate confidence in therapeutic analysis."""
        # Base confidence factors
        safety_confidence = safety.overall_safety_score / 100
        efficacy_confidence = efficacy.overall_efficacy
        strategy_confidence = min(
            1.0, num_strategies / 3
        )  # More strategies = higher confidence

        # Weighted average
        overall_confidence = (
            safety_confidence * 0.4
            + efficacy_confidence * 0.4
            + strategy_confidence * 0.2
        )

        return max(0.1, min(1.0, overall_confidence))

    def calculate_polygenic_risk(
        self,
        genotype_data: Dict[str, str],
        patient_data: Dict[str, Any],
        disease: str = "rheumatoid_arthritis",
    ) -> Dict[str, Any]:
        """
        Calculate polygenic risk score from multiple variants.

        Args:
            genotype_data: Gene -> genotype mapping
            patient_data: Patient demographics
            disease: Target disease

        Returns:
            Polygenic risk analysis
        """
        try:
            genetic_risk_score = (
                self.multigene_risk_calculator.calculate_polygenic_risk_score(
                    genotype_data, patient_data, disease
                )
            )

            # Generate therapeutic recommendations for polygenic risk
            therapeutic_recommendations = (
                self._generate_polygenic_therapeutic_recommendations(
                    genetic_risk_score, patient_data
                )
            )

            return {
                "genetic_risk_score": genetic_risk_score,
                "therapeutic_recommendations": therapeutic_recommendations,
                "priority_genes": self._identify_priority_therapeutic_targets(
                    genetic_risk_score.individual_contributions
                ),
                "intervention_strategy": self._design_polygenic_intervention_strategy(
                    genotype_data, genetic_risk_score
                ),
            }

        except Exception as e:
            self.logger.error(f"Polygenic risk calculation failed: {e}")
            raise TherapeuticEngineError(f"Polygenic risk calculation failed: {e}")

    def _generate_polygenic_therapeutic_recommendations(
        self, genetic_risk_score, patient_data: Dict[str, Any]
    ) -> List[str]:
        """Generate recommendations for polygenic risk management."""
        recommendations = []

        if genetic_risk_score.risk_category == "Very High Risk":
            recommendations.extend(
                [
                    "Consider multi-gene therapeutic intervention",
                    "Prioritize aggressive preventive measures",
                    "Frequent monitoring and early intervention protocols",
                ]
            )
        elif genetic_risk_score.risk_category == "High Risk":
            recommendations.extend(
                [
                    "Target highest-risk genetic variants for intervention",
                    "Enhanced monitoring and preventive care",
                    "Consider prophylactic treatment",
                ]
            )
        else:
            recommendations.extend(
                [
                    "Standard monitoring protocols",
                    "Lifestyle interventions and preventive care",
                    "Consider genetic counseling for family planning",
                ]
            )

        return recommendations

    def _identify_priority_therapeutic_targets(
        self, individual_contributions: Dict[str, float]
    ) -> List[str]:
        """Identify priority genes for therapeutic targeting."""
        # Sort genes by contribution to risk
        sorted_genes = sorted(
            individual_contributions.items(), key=lambda x: x[1], reverse=True
        )

        # Return top 3 contributors
        return [gene for gene, contribution in sorted_genes[:3] if contribution > 1.5]

    def _design_polygenic_intervention_strategy(
        self, genotype_data: Dict[str, str], genetic_risk_score
    ) -> Dict[str, Any]:
        """Design intervention strategy for polygenic risk."""
        strategy = {
            "approach": "sequential_targeting",
            "target_order": [],
            "expected_risk_reduction": 0.0,
            "intervention_complexity": "moderate",
        }

        # Prioritize targets by contribution
        priority_targets = self._identify_priority_therapeutic_targets(
            genetic_risk_score.individual_contributions
        )

        strategy["target_order"] = priority_targets

        # Estimate risk reduction
        total_reduction = 0
        for gene in priority_targets:
            contribution = genetic_risk_score.individual_contributions.get(gene, 1.0)
            if contribution > 1.5:
                total_reduction += (
                    contribution - 1.0
                ) * 0.7  # 70% correction efficiency

        strategy["expected_risk_reduction"] = min(0.8, total_reduction)

        # Assess complexity
        if len(priority_targets) > 2:
            strategy["intervention_complexity"] = "high"
        elif len(priority_targets) == 0:
            strategy["intervention_complexity"] = "low"

        return strategy

    def compare_therapeutic_options(
        self,
        gene_list: List[str],
        patient_data: Dict[str, Any],
        disease: str = "rheumatoid_arthritis",
    ) -> Dict[str, Any]:
        """
        Compare therapeutic options for multiple genes.

        Args:
            gene_list: List of genes to analyze
            patient_data: Patient data
            disease: Target disease

        Returns:
            Comparative analysis of therapeutic options
        """
        try:
            therapeutic_analyses = {}

            for gene in gene_list:
                # Get main pathogenic variant for gene
                variant = self._get_main_pathogenic_variant(gene, disease)
                if variant:
                    try:
                        analysis = self.analyze_disease_gene(
                            gene, variant, patient_data, disease
                        )
                        therapeutic_analyses[gene] = analysis
                    except Exception as e:
                        self.logger.warning(f"Failed to analyze {gene}: {e}")

            # Compare analyses
            comparison = self._perform_therapeutic_comparison(therapeutic_analyses)

            return comparison

        except Exception as e:
            self.logger.error(f"Therapeutic comparison failed: {e}")
            raise TherapeuticEngineError(f"Therapeutic comparison failed: {e}")

    def _get_main_pathogenic_variant(self, gene: str, disease: str) -> Optional[str]:
        """Get main pathogenic variant for gene-disease combination."""
        # Look up in disease database
        for category in DISEASE_GENES.values():
            if gene in category:
                gene_data = category[gene]
                pathogenic_variants = gene_data.get("pathogenic_variants", {})
                if pathogenic_variants:
                    return list(pathogenic_variants.keys())[0]
        return None

    def _perform_therapeutic_comparison(
        self, analyses: Dict[str, TherapeuticReport]
    ) -> Dict[str, Any]:
        """Perform comparative analysis of therapeutic options."""
        if not analyses:
            return {"error": "No successful analyses to compare"}

        # Extract key metrics
        comparison_data = {}
        for gene, report in analyses.items():
            comparison_data[gene] = {
                "safety_score": report.safety_assessment.overall_safety_score,
                "efficacy_score": report.predicted_efficacy.overall_efficacy * 100,
                "treatment_recommendation": report.treatment_recommendation,
                "confidence": report.confidence_score,
                "reversibility": report.safety_assessment.reversibility_score,
                "time_to_effect": report.predicted_efficacy.time_to_effect_days,
                "persistence": report.predicted_efficacy.persistence_months,
            }

        # Rank options
        rankings = {
            "by_safety": sorted(
                analyses.keys(),
                key=lambda x: comparison_data[x]["safety_score"],
                reverse=True,
            ),
            "by_efficacy": sorted(
                analyses.keys(),
                key=lambda x: comparison_data[x]["efficacy_score"],
                reverse=True,
            ),
            "by_overall_score": sorted(
                analyses.keys(),
                key=lambda x: (
                    comparison_data[x]["safety_score"] * 0.6
                    + comparison_data[x]["efficacy_score"] * 0.4
                ),
                reverse=True,
            ),
        }

        # Generate recommendations
        best_option = (
            rankings["by_overall_score"][0] if rankings["by_overall_score"] else None
        )

        return {
            "comparison_data": comparison_data,
            "rankings": rankings,
            "best_overall_option": best_option,
            "recommendation": self._generate_comparison_recommendation(
                comparison_data, rankings
            ),
            "summary_statistics": self._calculate_comparison_statistics(
                comparison_data
            ),
        }

    def _generate_comparison_recommendation(
        self, comparison_data: Dict, rankings: Dict
    ) -> str:
        """Generate recommendation from comparison."""
        if not rankings["by_overall_score"]:
            return "No suitable therapeutic options identified"

        best_gene = rankings["by_overall_score"][0]
        best_data = comparison_data[best_gene]

        if best_data["safety_score"] > 70 and best_data["efficacy_score"] > 60:
            return f"Recommend {best_gene} as primary therapeutic target - excellent safety and efficacy profile"
        elif best_data["safety_score"] > 50 and best_data["efficacy_score"] > 40:
            return f"Consider {best_gene} with enhanced monitoring - moderate benefit-risk ratio"
        else:
            return "All options show limited therapeutic potential - consider alternative approaches"

    def _calculate_comparison_statistics(self, comparison_data: Dict) -> Dict[str, Any]:
        """Calculate summary statistics for comparison."""
        if not comparison_data:
            return {}

        safety_scores = [data["safety_score"] for data in comparison_data.values()]
        efficacy_scores = [data["efficacy_score"] for data in comparison_data.values()]

        return {
            "mean_safety_score": sum(safety_scores) / len(safety_scores),
            "mean_efficacy_score": sum(efficacy_scores) / len(efficacy_scores),
            "high_safety_options": sum(1 for score in safety_scores if score > 70),
            "high_efficacy_options": sum(1 for score in efficacy_scores if score > 60),
            "recommended_options": sum(
                1
                for gene, data in comparison_data.items()
                if "recommend" in data["treatment_recommendation"].lower()
            ),
        }

    def export_therapeutic_report(
        self, report: TherapeuticReport, format: str = "json"
    ) -> str:
        """
        Export therapeutic report in specified format.

        Args:
            report: Therapeutic report to export
            format: Export format ("json", "clinical_summary")

        Returns:
            Exported content or filename
        """
        try:
            timestamp = current_timestamp()

            if format == "json":
                # Export full report as JSON
                filename = f"therapeutic_report_{report.patient_id}_{timestamp}.json"

                # Convert report to serializable format
                report_dict = {
                    "patient_id": report.patient_id,
                    "target_disease": report.target_disease,
                    "analysis_date": report.analysis_date.isoformat(),
                    "treatment_recommendation": report.treatment_recommendation,
                    "priority_level": report.priority_level,
                    "safety_assessment": {
                        "overall_score": report.safety_assessment.overall_safety_score,
                        "risk_category": report.safety_assessment.risk_category,
                        "reversibility": report.safety_assessment.reversibility_score,
                    },
                    "predicted_efficacy": {
                        "overall_efficacy": report.predicted_efficacy.overall_efficacy,
                        "efficacy_category": report.predicted_efficacy.efficacy_category,
                        "time_to_effect": report.predicted_efficacy.time_to_effect_days,
                        "persistence": report.predicted_efficacy.persistence_months,
                    },
                    "recommendations": report.recommendations,
                    "contraindications": report.contraindications,
                    "monitoring_plan": report.monitoring_plan,
                    "confidence_score": report.confidence_score,
                }

                save_json(report_dict, filename)
                return filename

            elif format == "clinical_summary":
                # Generate clinical summary
                summary = self._generate_clinical_summary(report)
                filename = f"clinical_summary_{report.patient_id}_{timestamp}.txt"

                with open(filename, "w") as f:
                    f.write(summary)

                return filename

            else:
                raise TherapeuticEngineError(f"Unsupported export format: {format}")

        except Exception as e:
            self.logger.error(f"Report export failed: {e}")
            raise TherapeuticEngineError(f"Report export failed: {e}")

    def _generate_clinical_summary(self, report: TherapeuticReport) -> str:
        """Generate clinical summary of therapeutic report."""
        summary = f"""
THERAPEUTIC GENETIC INTERVENTION REPORT

Patient ID: {report.patient_id}
Target Disease: {report.target_disease}
Analysis Date: {report.analysis_date.strftime('%Y-%m-%d')}

TREATMENT RECOMMENDATION: {report.treatment_recommendation}
Priority Level: {report.priority_level}

SAFETY ASSESSMENT:
- Overall Safety Score: {report.safety_assessment.overall_safety_score:.1f}/100
- Risk Category: {report.safety_assessment.risk_category}
- Reversibility Score: {report.safety_assessment.reversibility_score:.2f}

PREDICTED EFFICACY:
- Overall Efficacy: {report.predicted_efficacy.overall_efficacy:.2f}
- Efficacy Category: {report.predicted_efficacy.efficacy_category}
- Time to Effect: {report.predicted_efficacy.time_to_effect_days} days
- Duration of Effect: {report.predicted_efficacy.persistence_months:.1f} months

KEY RECOMMENDATIONS:
"""

        for i, rec in enumerate(report.recommendations[:5], 1):
            summary += f"{i}. {rec}\n"

        if report.contraindications:
            summary += "\nCONTRAINDICATIONS:\n"
            for i, contra in enumerate(report.contraindications, 1):
                summary += f"{i}. {contra}\n"

        summary += f"\nCONFIDENCE SCORE: {report.confidence_score:.2f}\n"
        summary += "\nDISCLAIMER: This analysis is for research purposes only. Clinical decisions should involve qualified healthcare professionals.\n"

        return summary

    def clear_cache(self):
        """Clear therapeutic analysis cache."""
        self._therapeutic_cache.clear()
        if hasattr(self.base_engine, "clear_cache"):
            self.base_engine.clear_cache()
        self.logger.info("Therapeutic analysis cache cleared")


# Additional helper classes that would be imported
# Additional helper classes that would be imported
class PatientStratificationEngine:
    """Basic engine for stratifying patients for therapy."""

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def stratify(
        self, patient_data: Dict[str, Any], disease_risk: DiseaseRisk
    ) -> PatientStratification:
        """Create a simple stratification record."""
        severity = patient_data.get("disease_severity", "moderate")
        priority = min(1.0, disease_risk.odds_ratio / 10.0)
        if severity == "severe":
            priority = max(priority, 0.8)
        return PatientStratification(
            patient_genotype=patient_data.get("genotypes", {}),
            disease_severity=severity,
            prior_treatments=patient_data.get("prior_treatments", []),
            comorbidities=patient_data.get("comorbidities", []),
            ethnicity=patient_data.get("ethnicity", "unknown"),
            age=patient_data.get("age", 0),
            sex=patient_data.get("sex", "unknown"),
            treatment_priority=priority,
        )


class ClinicalAssessmentEngine:
    """Assess basic clinical features."""

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def assess(self, patient_data: Dict[str, Any]) -> Dict[str, Any]:
        """Return a lightweight assessment summary."""
        assessment = {
            "disease_severity": patient_data.get("disease_severity", "unknown"),
            "comorbidities": patient_data.get("comorbidities", []),
        }
        self.logger.debug("Clinical assessment generated: %s", assessment)
        return assessment


class TherapeuticSimulationEngine:
    """Minimal therapeutic simulation engine."""

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def run(self, strategy: CorrectionStrategy) -> Dict[str, float]:
        """Simulate a strategy and return mock metrics."""
        complexity = getattr(strategy, "complexity_score", 0.5)
        return {
            "predicted_success_rate": max(0.0, 1.0 - complexity * 0.5),
            "predicted_off_target_rate": complexity * 0.01,
        }
