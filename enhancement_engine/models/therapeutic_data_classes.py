"""
Therapeutic data classes for CRISPR disease correction.

This module defines data structures for therapeutic genetic interventions,
disease risk assessment, and corrective strategies.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any, Union
from datetime import datetime
from enum import Enum

from .data_classes import GeneInfo, VariantInfo, GuideRNA, SideEffect, RiskLevel


class DiseaseCategory(Enum):
    """Categories of genetic diseases."""
    AUTOIMMUNE = "autoimmune"
    CARDIOVASCULAR = "cardiovascular"
    NEUROLOGICAL = "neurological"
    METABOLIC = "metabolic"
    CANCER_PREDISPOSITION = "cancer_predisposition"
    INFLAMMATORY = "inflammatory"
    IMMUNODEFICIENCY = "immunodeficiency"


class TherapeuticApproach(Enum):
    """Types of therapeutic genetic interventions."""
    CORRECTION = "correction"          # Fix disease-causing mutation
    KNOCKOUT = "knockout"              # Remove harmful gene function
    REPLACEMENT = "replacement"        # Replace with protective allele
    SILENCING = "silencing"           # Reduce expression (CRISPRi)
    ACTIVATION = "activation"         # Increase expression (CRISPRa)
    EPIGENETIC = "epigenetic"         # Modify methylation/chromatin


class DeliveryMethod(Enum):
    """Therapeutic delivery methods."""
    INTRAVENOUS = "intravenous"
    INTRA_ARTICULAR = "intra_articular"    # For joint diseases
    INTRATHECAL = "intrathecal"            # For CNS diseases
    LOCAL_INJECTION = "local_injection"
    TOPICAL = "topical"
    INHALED = "inhaled"
    ORAL = "oral"


@dataclass
class DiseaseRisk:
    """Disease risk assessment for genetic variants."""
    disease_name: str
    odds_ratio: float
    population_frequency: float
    absolute_risk_increase: float
    confidence_interval: Tuple[float, float]
    ethnicity_specific: Dict[str, float] = field(default_factory=dict)
    age_dependent: bool = False
    sex_specific: bool = False
    penetrance: float = 1.0  # 0-1, fraction who develop disease given genotype
    
    @property
    def risk_level(self) -> RiskLevel:
        """Categorize risk level based on odds ratio."""
        if self.odds_ratio >= 10.0:
            return RiskLevel.VERY_HIGH
        elif self.odds_ratio >= 5.0:
            return RiskLevel.HIGH
        elif self.odds_ratio >= 2.0:
            return RiskLevel.MODERATE
        elif self.odds_ratio >= 1.5:
            return RiskLevel.LOW
        else:
            return RiskLevel.VERY_LOW
    
    @property
    def clinical_significance(self) -> str:
        """Clinical significance of the risk."""
        if self.odds_ratio >= 5.0:
            return "High clinical significance - genetic counseling recommended"
        elif self.odds_ratio >= 2.0:
            return "Moderate significance - discuss with healthcare provider"
        else:
            return "Low significance - routine monitoring sufficient"


@dataclass
class TherapeuticTarget:
    """Target definition for therapeutic intervention."""
    gene_symbol: str
    disease_association: str
    target_variant: str
    therapeutic_approach: TherapeuticApproach
    target_tissue: str
    delivery_method: DeliveryMethod
    correction_sequence: Optional[str] = None  # Desired sequence after correction
    target_position: int = 0
    priority_score: float = 0.0  # 0-1, higher = more important target
    
    @property
    def is_correctable(self) -> bool:
        """Check if target is amenable to correction."""
        correctable_approaches = {
            TherapeuticApproach.CORRECTION,
            TherapeuticApproach.REPLACEMENT
        }
        return self.therapeutic_approach in correctable_approaches


@dataclass
class TherapeuticEfficacy:
    """Predicted efficacy of therapeutic intervention."""
    correction_efficiency: float  # 0-1, fraction of cells corrected
    tissue_penetration: float     # 0-1, fraction of target tissue reached
    persistence_months: float     # Duration of effect
    clinical_improvement: float   # Expected improvement in disease score
    remission_probability: float  # 0-1, chance of achieving remission
    time_to_effect_days: int = 30
    
    @property
    def overall_efficacy(self) -> float:
        """Calculate overall therapeutic efficacy score."""
        return (
            self.correction_efficiency * 0.3 +
            self.tissue_penetration * 0.2 +
            min(1.0, self.persistence_months / 12) * 0.2 +
            self.clinical_improvement * 0.3
        )
    
    @property
    def efficacy_category(self) -> str:
        """Categorize overall efficacy."""
        score = self.overall_efficacy
        if score >= 0.8:
            return "Excellent"
        elif score >= 0.6:
            return "Good"
        elif score >= 0.4:
            return "Moderate"
        elif score >= 0.2:
            return "Limited"
        else:
            return "Poor"


@dataclass
class CorrectionStrategy:
    """Strategy for correcting disease-causing mutations."""
    original_sequence: str
    corrected_sequence: str
    editing_method: str  # "base_editing", "prime_editing", "HDR"
    guide_rnas: List[GuideRNA]
    template_sequence: Optional[str] = None  # For HDR
    correction_window: Tuple[int, int] = (0, 0)  # Start, end positions
    
    @property
    def mutation_type(self) -> str:
        """Classify the type of correction needed."""
        if len(self.original_sequence) == len(self.corrected_sequence) == 1:
            return "point_mutation"
        elif len(self.original_sequence) > len(self.corrected_sequence):
            return "deletion"
        elif len(self.original_sequence) < len(self.corrected_sequence):
            return "insertion"
        else:
            return "substitution"
    
    @property
    def complexity_score(self) -> float:
        """Score correction complexity (0-1, higher = more complex)."""
        base_scores = {
            "point_mutation": 0.2,
            "substitution": 0.4,
            "deletion": 0.6,
            "insertion": 0.8
        }
        
        base_score = base_scores.get(self.mutation_type, 0.5)
        
        # Adjust for length
        length_factor = min(1.0, abs(len(self.corrected_sequence) - len(self.original_sequence)) / 10)
        
        return min(1.0, base_score + length_factor * 0.3)


@dataclass
class TherapeuticSafety:
    """Safety assessment for therapeutic interventions."""
    off_target_risk: float       # 0-10 scale
    immunogenicity_risk: float   # 0-10 scale
    tissue_toxicity_risk: float  # 0-10 scale
    systemic_effects_risk: float # 0-10 scale
    reversibility_score: float   # 0-1, higher = more reversible
    contraindications: List[str] = field(default_factory=list)
    monitoring_requirements: List[str] = field(default_factory=list)
    
    @property
    def overall_safety_score(self) -> float:
        """Calculate overall safety score (0-100)."""
        risk_sum = (
            self.off_target_risk +
            self.immunogenicity_risk +
            self.tissue_toxicity_risk +
            self.systemic_effects_risk
        )
        
        # Convert to 0-100 scale (lower risk = higher score)
        safety_score = max(0, 100 - (risk_sum / 4) * 10)
        
        # Adjust for reversibility
        safety_score += self.reversibility_score * 10
        
        return min(100.0, safety_score)
    
    @property
    def risk_category(self) -> str:
        """Categorize overall risk."""
        score = self.overall_safety_score
        if score >= 90:
            return "Very Low Risk"
        elif score >= 70:
            return "Low Risk"
        elif score >= 50:
            return "Moderate Risk"
        elif score >= 30:
            return "High Risk"
        else:
            return "Very High Risk"


@dataclass
class PatientStratification:
    """Patient-specific treatment stratification."""
    patient_genotype: Dict[str, str]  # Gene -> genotype
    disease_severity: str             # "mild", "moderate", "severe"
    prior_treatments: List[str]
    comorbidities: List[str]
    ethnicity: str
    age: int
    sex: str
    treatment_priority: float = 0.0  # 0-1, urgency of intervention
    
    @property
    def treatment_candidacy(self) -> str:
        """Assess candidacy for genetic therapy."""
        if self.disease_severity == "severe" and self.treatment_priority > 0.7:
            return "Excellent candidate"
        elif self.disease_severity in ["moderate", "severe"] and self.treatment_priority > 0.5:
            return "Good candidate"
        elif self.treatment_priority > 0.3:
            return "Consider for treatment"
        else:
            return "Conservative management preferred"


@dataclass
class ClinicalEndpoint:
    """Clinical endpoints for therapeutic assessment."""
    endpoint_name: str
    baseline_value: float
    target_value: float
    measurement_unit: str
    time_to_assessment_weeks: int
    minimum_clinically_important_difference: float
    
    @property
    def expected_improvement(self) -> float:
        """Calculate expected improvement."""
        return abs(self.target_value - self.baseline_value)
    
    @property
    def is_clinically_significant(self) -> bool:
        """Check if expected improvement is clinically significant."""
        return self.expected_improvement >= self.minimum_clinically_important_difference


@dataclass
class TherapeuticReport:
    """Comprehensive therapeutic intervention report."""
    patient_id: str
    target_disease: str
    therapeutic_targets: List[TherapeuticTarget]
    correction_strategies: List[CorrectionStrategy]
    predicted_efficacy: TherapeuticEfficacy
    safety_assessment: TherapeuticSafety
    patient_stratification: PatientStratification
    clinical_endpoints: List[ClinicalEndpoint]
    recommendations: List[str] = field(default_factory=list)
    contraindications: List[str] = field(default_factory=list)
    monitoring_plan: List[str] = field(default_factory=list)
    follow_up_schedule: Dict[str, str] = field(default_factory=dict)
    analysis_date: datetime = field(default_factory=datetime.now)
    confidence_score: float = 0.0
    
    @property
    def treatment_recommendation(self) -> str:
        """Overall treatment recommendation."""
        efficacy_score = self.predicted_efficacy.overall_efficacy
        safety_score = self.safety_assessment.overall_safety_score / 100
        
        combined_score = (efficacy_score * 0.6 + safety_score * 0.4)
        
        if combined_score >= 0.8:
            return "Strongly recommended"
        elif combined_score >= 0.6:
            return "Recommended with monitoring"
        elif combined_score >= 0.4:
            return "Consider carefully - benefits may outweigh risks"
        else:
            return "Not recommended - risks outweigh benefits"
    
    @property
    def priority_level(self) -> str:
        """Treatment priority level."""
        urgency = self.patient_stratification.treatment_priority
        if urgency >= 0.8:
            return "Urgent"
        elif urgency >= 0.6:
            return "High"
        elif urgency >= 0.4:
            return "Medium"
        else:
            return "Low"


@dataclass
class ReversalStrategy:
    """Strategy for reversing therapeutic interventions."""
    reversal_method: str          # "counter_editing", "replacement", "silencing"
    reversal_guides: List[GuideRNA]
    success_probability: float    # 0-1
    time_to_reversal_days: int
    reversal_risks: List[str] = field(default_factory=list)
    
    @property
    def is_feasible(self) -> bool:
        """Check if reversal is feasible."""
        return self.success_probability >= 0.3 and self.time_to_reversal_days <= 90


# Type aliases for therapeutic applications
TherapeuticTargets = List[TherapeuticTarget]
CorrectionStrategies = List[CorrectionStrategy]
DiseaseRisks = List[DiseaseRisk]