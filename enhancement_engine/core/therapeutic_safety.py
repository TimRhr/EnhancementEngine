"""
Therapeutic safety analysis for CRISPR-based disease correction.

This module provides specialized safety assessment for therapeutic genetic
interventions, including intervention-specific risks and monitoring requirements.
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any, Set
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum

from .safety import SafetyAnalyzer, OffTargetAnalyzer, GenotoxicityAssessor
from ..models.data_classes import GuideRNA, SafetyScore, SideEffect, RiskLevel
from ..models.therapeutic_data_classes import (
    TherapeuticSafety, TherapeuticTarget, CorrectionStrategy,
    TherapeuticApproach, DeliveryMethod, DiseaseCategory, ReversalStrategy
)
from ..models.disease_constants import THERAPEUTIC_SAFETY_THRESHOLDS


class TherapeuticSafetyError(Exception):
    """Custom exception for therapeutic safety analysis errors."""
    pass


@dataclass
class InterventionRisk:
    """Risk assessment for specific therapeutic intervention."""
    intervention_type: str
    target_tissue: str
    delivery_method: str
    acute_risks: Dict[str, float]      # Risk type -> probability
    chronic_risks: Dict[str, float]    # Risk type -> probability
    reversibility: float               # 0-1, higher = more reversible
    monitoring_frequency: str          # "daily", "weekly", "monthly"
    contraindications: List[str]
    special_populations: Dict[str, str]  # Population -> risk modifier


class TherapeuticOffTargetAnalyzer:
    """Specialized off-target analysis for therapeutic applications."""
    
    def __init__(self):
        """Initialize therapeutic off-target analyzer."""
        self.base_analyzer = OffTargetAnalyzer()
        self.logger = logging.getLogger(__name__)
    
    def analyze_therapeutic_off_targets(self, guide_rna: GuideRNA,
                                      target_disease: str,
                                      correction_strategy: CorrectionStrategy) -> Dict[str, Any]:
        """
        Analyze off-targets with therapeutic context.
        
        Args:
            guide_rna: Guide RNA to analyze
            target_disease: Disease being treated
            correction_strategy: Therapeutic correction strategy
            
        Returns:
            Therapeutic off-target analysis
        """
        try:
            # Standard off-target analysis
            base_analysis = self.base_analyzer.analyze_off_targets(guide_rna)
            
            # Add therapeutic-specific assessments
            therapeutic_risks = self._assess_therapeutic_off_target_risks(
                guide_rna.off_targets, target_disease, correction_strategy
            )
            
            # Evaluate intervention-specific consequences
            intervention_consequences = self._evaluate_intervention_consequences(
                guide_rna.off_targets, correction_strategy.editing_method
            )
            
            # Generate therapeutic recommendations
            therapeutic_recommendations = self._generate_therapeutic_recommendations(
                base_analysis, therapeutic_risks, intervention_consequences
            )
            
            # Combine analyses
            therapeutic_analysis = base_analysis.copy()
            therapeutic_analysis.update({
                "therapeutic_risks": therapeutic_risks,
                "intervention_consequences": intervention_consequences,
                "therapeutic_recommendations": therapeutic_recommendations,
                "clinical_monitoring": self._generate_monitoring_plan(therapeutic_risks),
                "patient_counseling": self._generate_counseling_points(therapeutic_risks)
            })
            
            return therapeutic_analysis
            
        except Exception as e:
            self.logger.error(f"Therapeutic off-target analysis failed: {e}")
            raise TherapeuticSafetyError(f"Failed to analyze therapeutic off-targets: {e}")
    
    def _assess_therapeutic_off_target_risks(self, off_targets: List, target_disease: str,
                                           correction_strategy: CorrectionStrategy) -> Dict[str, Any]:
        """Assess off-target risks in therapeutic context."""
        risks = {
            "disease_exacerbation": 0.0,
            "treatment_interference": 0.0,
            "new_disease_risk": 0.0,
            "immune_complications": 0.0
        }
        
        # Disease-specific risk assessment
        if target_disease == "rheumatoid_arthritis":
            risks.update(self._assess_autoimmune_off_target_risks(off_targets))
        elif target_disease == "familial_hypercholesterolemia":
            risks.update(self._assess_cardiovascular_off_target_risks(off_targets))
        
        # Intervention-specific risks
        if correction_strategy.editing_method == "base_editing":
            risks["bystander_editing"] = self._assess_bystander_editing_risk(off_targets)
        elif correction_strategy.editing_method == "HDR":
            risks["integration_errors"] = self._assess_integration_error_risk(off_targets)
        
        return risks

    
    def _assess_autoimmune_off_target_risks(self, off_targets: List) -> Dict[str, float]:
        """Assess off-target risks specific to autoimmune diseases."""
        autoimmune_genes = {
            "TNF", "IL1B", "IL6", "IFNG", "TLR4", "MYD88", "NFKB1", "STAT1", "STAT3"
        }
        
        autoimmune_risk = 0.0
        immune_dysfunction_risk = 0.0
        
        for ot in off_targets:
            if ot.gene_context and any(gene in ot.gene_context.upper() for gene in autoimmune_genes):
                # Higher risk if off-target in immune genes
                risk_contribution = (1.0 / (1.0 + ot.mismatches)) * 0.3
                autoimmune_risk += risk_contribution
                
                if ot.mismatches <= 1:
                    immune_dysfunction_risk += 0.5
        
        return {
            "autoimmune_exacerbation": min(1.0, autoimmune_risk),
            "immune_dysfunction": min(1.0, immune_dysfunction_risk)
        }
    
    def _assess_cardiovascular_off_target_risks(self, off_targets: List) -> Dict[str, float]:
        """Assess off-target risks for cardiovascular interventions."""
        cardio_genes = {
            "APOB", "LDLR", "PCSK9", "HMGCR", "NPC1L1", "ABCG8", "CYP7A1"
        }
        
        lipid_metabolism_risk = 0.0
        cardiovascular_risk = 0.0
        
        for ot in off_targets:
            if ot.gene_context and any(gene in ot.gene_context.upper() for gene in cardio_genes):
                risk_contribution = (1.0 / (1.0 + ot.mismatches)) * 0.4
                lipid_metabolism_risk += risk_contribution
                
                if ot.mismatches == 0:
                    cardiovascular_risk += 0.8
        
        return {
            "lipid_dysregulation": min(1.0, lipid_metabolism_risk),
            "cardiovascular_events": min(1.0, cardiovascular_risk)
        }
    
    def _assess_bystander_editing_risk(self, off_targets: List) -> float:
        """Assess risk of unintended bystander edits."""
        # Base editors can cause bystander edits in off-target sites
        bystander_risk = 0.0
        
        for ot in off_targets:
            if ot.mismatches <= 2:  # High-confidence off-targets
                # Risk depends on number of editable bases in window
                bystander_risk += 0.1 * (3 - ot.mismatches)
        
        return min(1.0, bystander_risk)
    
    def _assess_integration_error_risk(self, off_targets: List) -> float:
        """Assess risk of HDR template integration errors."""
        integration_risk = 0.0
        
        for ot in off_targets:
            if ot.mismatches <= 1:
                # Risk of template integrating at off-target site
                integration_risk += 0.05
        
        return min(1.0, integration_risk)
    
    def _evaluate_intervention_consequences(self, off_targets: List, 
                                          editing_method: str) -> Dict[str, Any]:
        """Evaluate consequences of off-target effects for specific intervention."""
        consequences = {
            "functional_disruption": [],
            "structural_damage": [],
            "regulatory_interference": [],
            "protein_misfunction": []
        }
        
        for ot in off_targets:
            if ot.gene_context:
                # Categorize potential consequences
                if ot.mismatches <= 1:
                    consequences["functional_disruption"].append({
                        "gene": ot.gene_context,
                        "probability": 0.8,
                        "severity": "high"
                    })
                
                if editing_method in ["knockout", "HDR"]:
                    consequences["structural_damage"].append({
                        "gene": ot.gene_context,
                        "probability": 0.6,
                        "severity": "moderate"
                    })
        
        return consequences
    
    def _generate_therapeutic_recommendations(self, base_analysis: Dict,
                                            therapeutic_risks: Dict,
                                            consequences: Dict) -> List[str]:
        """Generate therapeutic-specific safety recommendations."""
        recommendations = []
        
        # Base recommendations from standard analysis
        recommendations.extend(base_analysis.get("recommendations", []))
        
        # Therapeutic-specific recommendations
        if therapeutic_risks.get("disease_exacerbation", 0) > 0.3:
            recommendations.append("Monitor for disease exacerbation during treatment")
            recommendations.append("Have disease-specific rescue therapies available")
        
        if therapeutic_risks.get("immune_complications", 0) > 0.2:
            recommendations.append("Consider immunosuppressive premedication")
            recommendations.append("Monitor immune function markers closely")
        
        if therapeutic_risks.get("bystander_editing", 0) > 0.1:
            recommendations.append("Use high-fidelity base editors to minimize bystander edits")
            recommendations.append("Sequence target region thoroughly post-treatment")
        
        # Consequences-based recommendations
        if consequences.get("functional_disruption"):
            recommendations.append("Perform functional assays for affected pathways")
        
        if consequences.get("structural_damage"):
            recommendations.append("Monitor for structural protein abnormalities")
        
        return recommendations
    
    def _generate_monitoring_plan(self, therapeutic_risks: Dict) -> Dict[str, Any]:
        """Generate clinical monitoring plan based on risks."""
        monitoring_plan = {
            "baseline_assessments": [],
            "acute_monitoring": [],      # First 24-48 hours
            "short_term_monitoring": [], # First month
            "long_term_monitoring": [],  # >1 month
            "frequency": {}
        }
        
        # Risk-based monitoring
        if therapeutic_risks.get("immune_complications", 0) > 0.2:
            monitoring_plan["baseline_assessments"].extend([
                "Complete blood count", "Immune panel", "Inflammatory markers"
            ])
            monitoring_plan["acute_monitoring"].append("Monitor for infusion reactions")
            monitoring_plan["short_term_monitoring"].append("Weekly immune function tests")
            monitoring_plan["frequency"]["immune_monitoring"] = "weekly_x4_then_monthly"
        
        if therapeutic_risks.get("disease_exacerbation", 0) > 0.3:
            monitoring_plan["baseline_assessments"].append("Disease activity scores")
            monitoring_plan["short_term_monitoring"].append("Disease activity monitoring")
            monitoring_plan["frequency"]["disease_monitoring"] = "weekly_x8"
        
        # Genomic monitoring
        monitoring_plan["long_term_monitoring"].extend([
            "Targeted sequencing of edited sites",
            "Off-target site verification",
            "Chromosomal stability assessment"
        ])
        monitoring_plan["frequency"]["genomic_monitoring"] = "1_month_3_months_1_year"
        
        return monitoring_plan
    
    def _generate_counseling_points(self, therapeutic_risks: Dict) -> List[str]:
        """Generate patient counseling points."""
        counseling_points = []
        
        # Risk communication
        if therapeutic_risks.get("disease_exacerbation", 0) > 0.2:
            counseling_points.append(
                "There is a risk that your disease symptoms may temporarily worsen during treatment"
            )
        
        if therapeutic_risks.get("immune_complications", 0) > 0.1:
            counseling_points.append(
                "You may experience immune-related side effects; report fever, rash, or unusual symptoms immediately"
            )
        
        # General counseling
        counseling_points.extend([
            "This is an experimental treatment with potential unknown long-term effects",
            "Regular monitoring will be required for safety assessment",
            "Notify healthcare providers of any unusual symptoms",
            "Genetic counseling is recommended for family planning decisions"
        ])
        
        return counseling_points


class TherapeuticGenotoxicityAssessor:
    """Specialized genotoxicity assessment for therapeutic applications."""
    
    def __init__(self):
        """Initialize therapeutic genotoxicity assessor."""
        self.base_assessor = GenotoxicityAssessor()
        self.logger = logging.getLogger(__name__)
    
    def assess_therapeutic_genotoxicity(self, correction_strategy: CorrectionStrategy,
                                      target_tissue: str, delivery_method: DeliveryMethod,
                                      patient_factors: Dict[str, Any]) -> Dict[str, Any]:
        """
        Assess genotoxicity in therapeutic context.
        
        Args:
            correction_strategy: Correction strategy being used
            target_tissue: Target tissue for intervention
            delivery_method: Delivery method
            patient_factors: Patient-specific factors
            
        Returns:
            Therapeutic genotoxicity assessment
        """
        try:
            # Base genotoxicity assessment
            base_assessment = {}
            if correction_strategy.guide_rnas:
                base_assessment = self.base_assessor.assess_genotoxicity(
                    correction_strategy.guide_rnas[0],  # Primary guide
                    "therapeutic_target",
                    correction_strategy.editing_method
                )
            
            # Therapeutic-specific factors
            tissue_specific_risks = self._assess_tissue_specific_genotoxicity(
                target_tissue, correction_strategy.editing_method
            )
            
            delivery_risks = self._assess_delivery_genotoxicity(
                delivery_method, correction_strategy.editing_method
            )
            
            patient_specific_risks = self._assess_patient_specific_genotoxicity(
                patient_factors, correction_strategy
            )
            
            # Integration-specific risks
            integration_risks = self._assess_integration_genotoxicity(
                correction_strategy
            )
            
            # Combined assessment
            therapeutic_assessment = base_assessment.copy()
            therapeutic_assessment.update({
                "tissue_specific_risks": tissue_specific_risks,
                "delivery_risks": delivery_risks,
                "patient_specific_risks": patient_specific_risks,
                "integration_risks": integration_risks,
                "cumulative_risk_score": self._calculate_cumulative_genotoxicity(
                    base_assessment, tissue_specific_risks, delivery_risks, patient_specific_risks
                )
            })
            
            return therapeutic_assessment
            
        except Exception as e:
            self.logger.error(f"Therapeutic genotoxicity assessment failed: {e}")
            raise TherapeuticSafetyError(f"Failed to assess therapeutic genotoxicity: {e}")
    
    def _assess_tissue_specific_genotoxicity(self, target_tissue: str, 
                                           editing_method: str) -> Dict[str, float]:
        """Assess genotoxicity risks specific to target tissue."""
        tissue_risks = {
            "proliferation_disruption": 0.0,
            "differentiation_interference": 0.0,
            "stem_cell_damage": 0.0,
            "tissue_specific_toxicity": 0.0
        }
        
        # Tissue-specific risk factors
        if target_tissue == "hematopoietic_system":
            tissue_risks["stem_cell_damage"] = 0.3  # HSCs are sensitive
            tissue_risks["proliferation_disruption"] = 0.2
            
            if editing_method in ["knockout", "HDR"]:
                tissue_risks["stem_cell_damage"] += 0.2
        
        elif target_tissue == "synovial_tissue":
            tissue_risks["tissue_specific_toxicity"] = 0.1  # Lower risk
            tissue_risks["differentiation_interference"] = 0.15
        
        elif target_tissue == "liver":
            tissue_risks["proliferation_disruption"] = 0.25  # High proliferation
            tissue_risks["tissue_specific_toxicity"] = 0.2   # Metabolic stress
        
        return tissue_risks
    
    def _assess_delivery_genotoxicity(self, delivery_method: DeliveryMethod,
                                    editing_method: str) -> Dict[str, float]:
        """Assess genotoxicity from delivery method."""
        delivery_risks = {
            "vector_integration": 0.0,
            "inflammatory_genotoxicity": 0.0,
            "delivery_stress": 0.0
        }
        
        if delivery_method == DeliveryMethod.INTRAVENOUS:
            delivery_risks["inflammatory_genotoxicity"] = 0.2
            delivery_risks["delivery_stress"] = 0.15
        
        elif delivery_method == DeliveryMethod.INTRA_ARTICULAR:
            delivery_risks["inflammatory_genotoxicity"] = 0.1
            delivery_risks["delivery_stress"] = 0.05
        
        # Viral vectors increase integration risk
        if "viral" in str(delivery_method).lower():
            delivery_risks["vector_integration"] = 0.4
        
        return delivery_risks
    
    def _assess_patient_specific_genotoxicity(self, patient_factors: Dict[str, Any],
                                            correction_strategy: CorrectionStrategy) -> Dict[str, float]:
        """Assess patient-specific genotoxicity risks."""
        patient_risks = {
            "age_related_sensitivity": 0.0,
            "comorbidity_interactions": 0.0,
            "genetic_predisposition": 0.0
        }
        
        age = patient_factors.get("age", 40)
        if age > 65:
            patient_risks["age_related_sensitivity"] = 0.3
        elif age < 18:
            patient_risks["age_related_sensitivity"] = 0.2
        
        # Comorbidity interactions
        comorbidities = patient_factors.get("comorbidities", [])
        if "cancer_history" in comorbidities:
            patient_risks["comorbidity_interactions"] = 0.4
            patient_risks["genetic_predisposition"] = 0.2
        
        if "immunodeficiency" in comorbidities:
            patient_risks["comorbidity_interactions"] = 0.3
        
        return patient_risks
    
    def _assess_integration_genotoxicity(self, correction_strategy: CorrectionStrategy) -> Dict[str, float]:
        """Assess genotoxicity from DNA integration events."""
        integration_risks = {
            "random_integration": 0.0,
            "chromosomal_rearrangement": 0.0,
            "template_mutagenesis": 0.0
        }
        
        if correction_strategy.editing_method == "HDR":
            integration_risks["random_integration"] = 0.3
            integration_risks["chromosomal_rearrangement"] = 0.2
            
            if correction_strategy.template_sequence and len(correction_strategy.template_sequence) > 1000:
                integration_risks["template_mutagenesis"] = 0.15
        
        elif correction_strategy.editing_method == "prime_editing":
            integration_risks["template_mutagenesis"] = 0.05
        
        return integration_risks
    
    def _calculate_cumulative_genotoxicity(self, base_assessment: Dict,
                                         tissue_risks: Dict, delivery_risks: Dict,
                                         patient_risks: Dict) -> float:
        """Calculate cumulative genotoxicity risk score."""
        base_score = base_assessment.get("overall_score", 0.0)
        
        tissue_score = sum(tissue_risks.values()) / len(tissue_risks) * 2
        delivery_score = sum(delivery_risks.values()) / len(delivery_risks) * 2
        patient_score = sum(patient_risks.values()) / len(patient_risks) * 2
        
        # Weighted combination
        cumulative_score = (
            base_score * 0.4 +
            tissue_score * 0.25 +
            delivery_score * 0.2 +
            patient_score * 0.15
        )
        
        return min(10.0, cumulative_score)


class TherapeuticSafetyAnalyzer:
    """Main therapeutic safety analyzer integrating all therapeutic safety components."""
    
    def __init__(self):
        """Initialize therapeutic safety analyzer."""
        self.base_safety_analyzer = SafetyAnalyzer()
        self.therapeutic_off_target = TherapeuticOffTargetAnalyzer()
        self.therapeutic_genotoxicity = TherapeuticGenotoxicityAssessor()
        self.logger = logging.getLogger(__name__)
    
    def analyze_therapeutic_safety(self, therapeutic_target: TherapeuticTarget,
                                 correction_strategies: List[CorrectionStrategy],
                                 patient_data: Dict[str, Any]) -> TherapeuticSafety:
        """
        Comprehensive therapeutic safety analysis.
        
        Args:
            therapeutic_target: Therapeutic target definition
            correction_strategies: List of correction strategies
            patient_data: Patient-specific data
            
        Returns:
            TherapeuticSafety assessment
        """
        try:
            if not correction_strategies:
                raise TherapeuticSafetyError("No correction strategies provided")
            
            # Analyze best strategy
            best_strategy = correction_strategies[0]  # Assume first is best
            primary_guide = best_strategy.guide_rnas[0] if best_strategy.guide_rnas else None
            
            # Off-target analysis
            off_target_risk = 0.0
            if primary_guide:
                off_target_analysis = self.therapeutic_off_target.analyze_therapeutic_off_targets(
                    primary_guide, therapeutic_target.disease_association, best_strategy
                )
                off_target_risk = off_target_analysis.get("risk_score", 0.0)
            
            # Genotoxicity analysis
            genotoxicity_analysis = self.therapeutic_genotoxicity.assess_therapeutic_genotoxicity(
                best_strategy, therapeutic_target.target_tissue,
                therapeutic_target.delivery_method, patient_data
            )
            genotoxicity_risk = genotoxicity_analysis.get("cumulative_risk_score", 0.0)
            
            # Immunogenicity analysis (simplified)
            immunogenicity_risk = self._assess_therapeutic_immunogenicity(
                best_strategy, therapeutic_target, patient_data
            )
            
            # Tissue-specific toxicity
            tissue_toxicity_risk = self._assess_tissue_toxicity(
                therapeutic_target, best_strategy, patient_data
            )
            
            # Systemic effects
            systemic_effects_risk = self._assess_systemic_effects(
                therapeutic_target, best_strategy
            )
            
            # Reversibility assessment
            reversibility_score = self._assess_reversibility(best_strategy, therapeutic_target)
            
            # Generate safety recommendations
            recommendations = self._generate_therapeutic_safety_recommendations(
                off_target_risk, genotoxicity_risk, immunogenicity_risk,
                tissue_toxicity_risk, systemic_effects_risk
            )
            
            # Monitoring requirements
            monitoring_requirements = self._generate_monitoring_requirements(
                therapeutic_target, best_strategy, off_target_risk, genotoxicity_risk
            )
            
            # Contraindications
            contraindications = self._identify_contraindications(
                therapeutic_target, patient_data, genotoxicity_risk
            )
            
            therapeutic_safety = TherapeuticSafety(
                off_target_risk=off_target_risk,
                immunogenicity_risk=immunogenicity_risk,
                tissue_toxicity_risk=tissue_toxicity_risk,
                systemic_effects_risk=systemic_effects_risk,
                reversibility_score=reversibility_score,
                contraindications=contraindications,
                monitoring_requirements=monitoring_requirements
            )
            
            self.logger.info(f"Therapeutic safety analysis completed: {therapeutic_safety.risk_category}")
            return therapeutic_safety
            
        except Exception as e:
            self.logger.error(f"Therapeutic safety analysis failed: {e}")
            raise TherapeuticSafetyError(f"Failed to analyze therapeutic safety: {e}")
    
    def _assess_therapeutic_immunogenicity(self, correction_strategy: CorrectionStrategy,
                                         therapeutic_target: TherapeuticTarget,
                                         patient_data: Dict[str, Any]) -> float:
        """Assess immunogenicity risk for therapeutic intervention."""
        base_immunogenicity = 3.0  # Base risk score
        
        # Delivery method affects immunogenicity
        delivery_modifiers = {
            DeliveryMethod.INTRAVENOUS: 1.5,
            DeliveryMethod.INTRA_ARTICULAR: 0.8,
            DeliveryMethod.LOCAL_INJECTION: 0.6,
            DeliveryMethod.TOPICAL: 0.4
        }
        
        delivery_modifier = delivery_modifiers.get(therapeutic_target.delivery_method, 1.0)
        immunogenicity_risk = base_immunogenicity * delivery_modifier
        
        # Patient factors
        if "immunocompromised" in patient_data.get("comorbidities", []):
            immunogenicity_risk *= 0.7  # Lower immune response
        
        if patient_data.get("age", 40) > 65:
            immunogenicity_risk *= 1.2  # Higher risk in elderly
        
        return min(10.0, immunogenicity_risk)
    
    def _assess_tissue_toxicity(self, therapeutic_target: TherapeuticTarget,
                              correction_strategy: CorrectionStrategy,
                              patient_data: Dict[str, Any]) -> float:
        """Assess tissue-specific toxicity risk."""
        base_toxicity = 2.0
        
        # Tissue-specific risk factors
        tissue_modifiers = {
            "hematopoietic_system": 2.0,  # HSCs are sensitive
            "synovial_tissue": 1.0,       # Local treatment, lower risk
            "liver": 1.5,                 # Metabolically active
            "muscle": 0.8,                # Generally well-tolerated
            "brain": 3.0                  # High risk, critical tissue
        }
        
        tissue_modifier = tissue_modifiers.get(therapeutic_target.target_tissue, 1.0)
        toxicity_risk = base_toxicity * tissue_modifier
        
        # Editing method modifier
        method_modifiers = {
            "base_editing": 0.8,
            "prime_editing": 1.0,
            "knockout": 1.5,
            "HDR": 2.0
        }
        
        method_modifier = method_modifiers.get(correction_strategy.editing_method, 1.0)
        toxicity_risk *= method_modifier
        
        return min(10.0, toxicity_risk)
    
    def _assess_systemic_effects(self, therapeutic_target: TherapeuticTarget,
                               correction_strategy: CorrectionStrategy) -> float:
        """Assess systemic effects risk."""
        base_systemic_risk = 1.5
        
        # Delivery method affects systemic exposure
        if therapeutic_target.delivery_method == DeliveryMethod.INTRAVENOUS:
            base_systemic_risk *= 2.0
        elif therapeutic_target.delivery_method in [DeliveryMethod.INTRA_ARTICULAR, DeliveryMethod.LOCAL_INJECTION]:
            base_systemic_risk *= 0.5
        
        # Target tissue affects systemic risk
        if therapeutic_target.target_tissue == "hematopoietic_system":
            base_systemic_risk *= 1.8  # Systemic circulation
        
        return min(10.0, base_systemic_risk)
    
    def _assess_reversibility(self, correction_strategy: CorrectionStrategy,
                            therapeutic_target: TherapeuticTarget) -> float:
        """Assess reversibility of therapeutic intervention."""
        base_reversibility = {
            "base_editing": 0.2,     # Hard to reverse point mutations
            "prime_editing": 0.3,    # Slightly more reversible
            "knockout": 0.1,         # Very hard to reverse
            "HDR": 0.4,              # Potentially reversible
            "CRISPRi": 0.9,          # Highly reversible
            "activation": 0.8        # Mostly reversible
        }
        
        reversibility = base_reversibility.get(correction_strategy.editing_method, 0.3)
        
        # Target tissue affects reversibility
        if therapeutic_target.target_tissue == "hematopoietic_system":
            reversibility *= 1.2  # Can replace with new cells
        
        return min(1.0, reversibility)
    
    def _generate_therapeutic_safety_recommendations(self, off_target_risk: float,
                                                   genotoxicity_risk: float,
                                                   immunogenicity_risk: float,
                                                   tissue_toxicity_risk: float,
                                                   systemic_effects_risk: float) -> List[str]:
        """Generate safety recommendations for therapeutic intervention."""
        recommendations = []
        
        # Risk-specific recommendations
        if off_target_risk > 3.0:
            recommendations.extend([
                "Perform comprehensive off-target validation before treatment",
                "Consider using high-fidelity editing tools",
                "Monitor for unintended genetic changes"
            ])
        
        if genotoxicity_risk > 5.0:
            recommendations.extend([
                "Implement extensive genotoxicity monitoring",
                "Consider fractionated dosing to reduce DNA damage",
                "Monitor chromosomal stability long-term"
            ])
        
        if immunogenicity_risk > 4.0:
            recommendations.extend([
                "Consider immunosuppressive premedication",
                "Monitor for immune reactions during and after treatment",
                "Have emergency protocols for severe immune reactions"
            ])
        
        if tissue_toxicity_risk > 3.0:
            recommendations.extend([
                "Implement tissue-specific monitoring protocols",
                "Consider dose reduction or alternative approaches",
                "Monitor organ function carefully"
            ])
        
        if systemic_effects_risk > 3.0:
            recommendations.extend([
                "Monitor systemic biodistribution",
                "Assess for off-organ effects",
                "Consider local delivery methods"
            ])
        
        # General therapeutic recommendations
        recommendations.extend([
            "Obtain comprehensive informed consent",
            "Establish clear safety monitoring protocols",
            "Maintain emergency response capabilities",
            "Document all adverse events thoroughly"
        ])
        
        return recommendations
    
    def _generate_monitoring_requirements(self, therapeutic_target: TherapeuticTarget,
                                        correction_strategy: CorrectionStrategy,
                                        off_target_risk: float,
                                        genotoxicity_risk: float) -> List[str]:
        """Generate monitoring requirements for therapeutic safety."""
        monitoring = []
        
        # Standard monitoring
        monitoring.extend([
            "Pre-treatment safety assessments",
            "Real-time vital signs monitoring during administration",
            "Post-treatment observation period"
        ])
        
        # Risk-based monitoring
        if off_target_risk > 2.0:
            monitoring.extend([
                "Targeted sequencing of predicted off-target sites",
                "Whole-genome sequencing at 1, 6, and 12 months"
            ])
        
        if genotoxicity_risk > 3.0:
            monitoring.extend([
                "Chromosomal stability assessment (karyotype)",
                "Micronucleus assay for DNA damage",
                "Long-term cancer surveillance"
            ])
        
        # Target-specific monitoring
        if therapeutic_target.target_tissue == "hematopoietic_system":
            monitoring.extend([
                "Complete blood count monitoring",
                "Bone marrow assessment if indicated",
                "Immune function testing"
            ])
        
        if therapeutic_target.disease_association == "rheumatoid_arthritis":
            monitoring.extend([
                "Disease activity scores (DAS28, CDAI)",
                "Joint imaging studies",
                "Inflammatory marker monitoring"
            ])
        
        return monitoring
    
    def _identify_contraindications(self, therapeutic_target: TherapeuticTarget,
                                  patient_data: Dict[str, Any],
                                  genotoxicity_risk: float) -> List[str]:
        """Identify contraindications for therapeutic intervention."""
        contraindications = []
        
        # Absolute contraindications
        if genotoxicity_risk > 7.0:
            contraindications.append("Prohibitively high genotoxicity risk")
        
        # Patient-specific contraindications
        comorbidities = patient_data.get("comorbidities", [])
        
        if "pregnancy" in comorbidities:
            contraindications.append("Pregnancy - unknown fetal effects")
        
        if "active_malignancy" in comorbidities:
            contraindications.append("Active malignancy - risk of genetic instability")
        
        if "severe_immunodeficiency" in comorbidities:
            contraindications.append("Severe immunodeficiency - unpredictable immune response")
        
        # Age-related contraindications
        age = patient_data.get("age", 40)
        if age < 16:
            contraindications.append("Pediatric population - developmental considerations")
        
        # Disease-specific contraindications
        if therapeutic_target.disease_association == "rheumatoid_arthritis":
            if "severe_joint_destruction" in comorbidities:
                contraindications.append("Advanced joint destruction - limited benefit expected")
        
        return contraindications
    
    def design_reversal_strategy(self, original_strategy: CorrectionStrategy,
                               therapeutic_target: TherapeuticTarget) -> Optional[ReversalStrategy]:
        """Design strategy for reversing therapeutic intervention."""
        try:
            reversal_method = self._determine_reversal_method(original_strategy)
            
            if not reversal_method:
                return None
            
            # Design reversal guides
            reversal_guides = self._design_reversal_guides(original_strategy, reversal_method)
            
            # Estimate success probability
            success_probability = self._estimate_reversal_success(
                original_strategy, reversal_method
            )
            
            # Estimate time to reversal
            time_to_reversal = self._estimate_reversal_time(reversal_method)
            
            # Identify reversal risks
            reversal_risks = self._identify_reversal_risks(original_strategy, reversal_method)
            
            reversal_strategy = ReversalStrategy(
                reversal_method=reversal_method,
                reversal_guides=reversal_guides,
                success_probability=success_probability,
                time_to_reversal_days=time_to_reversal,
                reversal_risks=reversal_risks
            )
            
            return reversal_strategy
            
        except Exception as e:
            self.logger.error(f"Reversal strategy design failed: {e}")
            return None
    
    def _determine_reversal_method(self, original_strategy: CorrectionStrategy) -> Optional[str]:
        """Determine appropriate reversal method."""
        reversal_methods = {
            "base_editing": "counter_base_editing",
            "prime_editing": "counter_prime_editing",
            "knockout": "gene_replacement",
            "HDR": "targeted_excision",
            "CRISPRi": "discontinue_treatment",
            "activation": "counter_repression"
        }
        
        return reversal_methods.get(original_strategy.editing_method)
    
    def _design_reversal_guides(self, original_strategy: CorrectionStrategy,
                              reversal_method: str) -> List[GuideRNA]:
        """Design guide RNAs for reversal strategy."""
        # Simplified reversal guide design
        # In practice, would need sophisticated design algorithms
        
        if reversal_method == "counter_base_editing":
            # Design guides to reverse the base edit
            return original_strategy.guide_rnas[:1]  # Reuse original guide
        
        elif reversal_method == "gene_replacement":
            # Design guides for inserting functional gene copy
            return []  # Would design new guides
        
        elif reversal_method == "targeted_excision":
            # Design guides to excise integrated sequence
            return []  # Would design flanking guides
        
        return []
    
    def _estimate_reversal_success(self, original_strategy: CorrectionStrategy,
                                 reversal_method: str) -> float:
        """Estimate probability of successful reversal."""
        base_success_rates = {
            "counter_base_editing": 0.3,
            "counter_prime_editing": 0.4,
            "gene_replacement": 0.2,
            "targeted_excision": 0.3,
            "discontinue_treatment": 0.9,
            "counter_repression": 0.7
        }
        
        return base_success_rates.get(reversal_method, 0.2)
    
    def _estimate_reversal_time(self, reversal_method: str) -> int:
        """Estimate time required for reversal (days)."""
        reversal_times = {
            "counter_base_editing": 30,
            "counter_prime_editing": 45,
            "gene_replacement": 60,
            "targeted_excision": 30,
            "discontinue_treatment": 7,
            "counter_repression": 14
        }
        
        return reversal_times.get(reversal_method, 30)
    
    def _identify_reversal_risks(self, original_strategy: CorrectionStrategy,
                               reversal_method: str) -> List[str]:
        """Identify risks associated with reversal strategy."""
        risks = []
        
        if reversal_method in ["counter_base_editing", "counter_prime_editing"]:
            risks.extend([
                "Additional off-target effects from reversal editing",
                "Incomplete reversal leaving mosaic state",
                "Unintended bystander edits during reversal"
            ])
        
        elif reversal_method == "gene_replacement":
            risks.extend([
                "Integration complications",
                "Immune response to replacement gene",
                "Disruption of endogenous gene regulation"
            ])
        
        elif reversal_method == "targeted_excision":
            risks.extend([
                "Chromosomal instability from double-strand breaks",
                "Incomplete excision leaving fragments",
                "Disruption of flanking genes"
            ])
        
        # General reversal risks
        risks.extend([
            "Disease progression during reversal period",
            "Psychological impact of treatment failure",
            "Delay in alternative treatment approaches"
        ])
        
        return risks
