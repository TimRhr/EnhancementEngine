"""
Data classes for Enhancement Engine.

This module defines all data structures used throughout the Enhancement Engine
for representing genes, CRISPR guides, analysis results, and reports.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any, Union
from datetime import datetime
from enum import Enum


class EnhancementCategory(Enum):
    """Categories of genetic enhancement."""
    COGNITIVE = "cognitive"
    PHYSICAL = "physical"
    LONGEVITY = "longevity"
    SENSORY = "sensory"
    METABOLIC = "metabolic"
    IMMUNE = "immune"


class RiskLevel(Enum):
    """Risk assessment levels."""
    VERY_LOW = "very_low"
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    VERY_HIGH = "very_high"


class CasType(Enum):
    """Supported Cas protein types."""
    CAS9 = "cas9"
    CAS12A = "cas12a"
    CAS13 = "cas13"
    BASE_EDITOR = "base_editor"
    PRIME_EDITOR = "prime_editor"


@dataclass
class GeneInfo:
    """Information about a gene."""
    name: str
    symbol: str
    gene_id: str
    chromosome: str
    start_pos: int
    end_pos: int
    description: str
    aliases: List[str] = field(default_factory=list)
    organism: str = "homo sapiens"
    ncbi_id: Optional[str] = None
    refseq_id: Optional[str] = None
    enhancement_category: Optional[EnhancementCategory] = None
    
    @property
    def length(self) -> int:
        """Gene length in base pairs."""
        return self.end_pos - self.start_pos + 1


@dataclass
class VariantInfo:
    """Information about a genetic variant."""
    name: str
    rsid: Optional[str] = None
    chromosome: str = ""
    position: int = 0
    ref_allele: str = ""
    alt_allele: str = ""
    variant_type: str = "SNP"  # SNP, INDEL, CNV, etc.
    clinical_significance: Optional[str] = None
    allele_frequency: Optional[float] = None
    
    def __str__(self) -> str:
        if self.rsid:
            return f"{self.name} ({self.rsid})"
        return self.name


@dataclass
class PAMSite:
    """PAM site information."""
    sequence: str
    position: int
    strand: str  # "+" or "-"
    cas_type: CasType
    distance_to_target: int = 0
    
    @property
    def is_valid(self) -> bool:
        """Check if PAM site is valid for the Cas type."""
        pam_patterns = {
            CasType.CAS9: "NGG",
            CasType.CAS12A: "TTTV",
        }
        pattern = pam_patterns.get(self.cas_type)
        if not pattern:
            return True  # Unknown patterns are considered valid
        
        # Simple pattern matching (N=any, V=A,C,G)
        if len(self.sequence) != len(pattern):
            return False
            
        for i, (seq_char, pat_char) in enumerate(zip(self.sequence, pattern)):
            if pat_char == 'N':
                continue
            elif pat_char == 'V' and seq_char in 'ACG':
                continue
            elif seq_char != pat_char:
                return False
        return True


@dataclass
class OffTarget:
    """Off-target binding site information."""
    sequence: str
    chromosome: str
    position: int
    strand: str
    mismatches: int
    bulges: int = 0
    score: float = 0.0
    gene_context: Optional[str] = None
    essential_gene: bool = False
    
    @property
    def risk_level(self) -> RiskLevel:
        """Assess risk level based on mismatches and context."""
        if self.essential_gene and self.mismatches <= 2:
            return RiskLevel.VERY_HIGH
        elif self.mismatches == 0:
            return RiskLevel.VERY_HIGH
        elif self.mismatches == 1:
            return RiskLevel.HIGH
        elif self.mismatches == 2:
            return RiskLevel.MODERATE
        elif self.mismatches >= 3:
            return RiskLevel.LOW
        return RiskLevel.MODERATE


@dataclass
class EfficiencyScore:
    """CRISPR guide efficiency metrics."""
    on_target_score: float
    doench_score: Optional[float] = None
    azimuth_score: Optional[float] = None
    chopchop_score: Optional[float] = None
    combined_score: Optional[float] = None
    
    @property
    def overall_efficiency(self) -> float:
        """Calculate overall efficiency score."""
        if self.combined_score is not None:
            return self.combined_score
        
        scores = [s for s in [self.on_target_score, self.doench_score, 
                             self.azimuth_score, self.chopchop_score] if s is not None]
        
        if not scores:
            return 0.0
        return sum(scores) / len(scores)


@dataclass
class GuideRNA:
    """CRISPR guide RNA information."""
    sequence: str
    pam_site: PAMSite
    target_position: int
    efficiency_score: EfficiencyScore
    gc_content: float
    off_targets: List[OffTarget] = field(default_factory=list)
    specificity_score: float = 0.0
    cas_type: CasType = CasType.CAS9
    strand: str = "+"
    
    @property
    def length(self) -> int:
        """Guide RNA length."""
        return len(self.sequence)
    
    @property
    def total_off_targets(self) -> int:
        """Total number of predicted off-targets."""
        return len(self.off_targets)
    
    @property
    def high_risk_off_targets(self) -> List[OffTarget]:
        """Off-targets with high risk (≤2 mismatches)."""
        return [ot for ot in self.off_targets if ot.mismatches <= 2]
    
    @property
    def overall_safety_score(self) -> float:
        """Calculate overall safety score (0-100)."""
        if not self.off_targets:
            return 100.0
        
        # Penalize based on number and severity of off-targets
        high_risk_count = len(self.high_risk_off_targets)
        total_count = len(self.off_targets)
        
        # Base score starts at 100
        safety_score = 100.0
        
        # Subtract points for high-risk off-targets
        safety_score -= high_risk_count * 20
        
        # Subtract points for total off-targets
        safety_score -= (total_count - high_risk_count) * 5
        
        # Essential gene penalty
        essential_hits = sum(1 for ot in self.off_targets if ot.essential_gene)
        safety_score -= essential_hits * 30
        
        return max(0.0, min(100.0, safety_score))


@dataclass
class CodingRegion:
    """Coding sequence region information."""
    start: int
    end: int
    exons: List[Tuple[int, int]] = field(default_factory=list)
    introns: List[Tuple[int, int]] = field(default_factory=list)
    
    @property
    def length(self) -> int:
        """Total coding region length."""
        return self.end - self.start + 1
    
    @property
    def exon_count(self) -> int:
        """Number of exons."""
        return len(self.exons)


@dataclass
class VariantPosition:
    """Position information for a genetic variant."""
    gene_position: int
    codon_position: int
    codon_number: int
    amino_acid_position: int
    reference_codon: str
    alternative_codon: str
    reference_aa: str
    alternative_aa: str
    
    @property
    def is_synonymous(self) -> bool:
        """Check if variant is synonymous (no amino acid change)."""
        return self.reference_aa == self.alternative_aa
    
    @property
    def variant_type(self) -> str:
        """Classify variant type."""
        if self.is_synonymous:
            return "synonymous"
        elif self.alternative_aa == "*":
            return "nonsense"
        elif self.reference_aa == "*":
            return "readthrough"
        else:
            return "missense"


@dataclass
class ProteinEffect:
    """Predicted protein-level effects of a variant."""
    stability_change: float  # ΔΔG in kcal/mol
    activity_change: float   # Relative activity change
    structure_disrupted: bool = False
    domain_affected: Optional[str] = None
    confidence_score: float = 0.0
    
    @property
    def effect_magnitude(self) -> str:
        """Classify effect magnitude."""
        abs_stability = abs(self.stability_change)
        if abs_stability > 2.0:
            return "large"
        elif abs_stability > 1.0:
            return "moderate"
        elif abs_stability > 0.5:
            return "small"
        else:
            return "minimal"


@dataclass
class EnhancementGain:
    """Quantified enhancement benefits."""
    category: EnhancementCategory
    primary_metric: str
    baseline_value: float
    enhanced_value: float
    improvement_factor: float
    confidence_interval: Tuple[float, float] = (0.0, 0.0)
    population_percentile: Optional[float] = None
    
    @property
    def absolute_gain(self) -> float:
        """Absolute improvement value."""
        return self.enhanced_value - self.baseline_value
    
    @property
    def relative_gain(self) -> float:
        """Relative improvement as percentage."""
        if self.baseline_value == 0:
            return 0.0
        return (self.absolute_gain / self.baseline_value) * 100
    
    @property
    def effect_size_category(self) -> str:
        """Categorize effect size."""
        if self.improvement_factor >= 2.0:
            return "major"
        elif self.improvement_factor >= 1.5:
            return "moderate"
        elif self.improvement_factor >= 1.2:
            return "minor"
        else:
            return "negligible"


@dataclass
class SideEffect:
    """Potential side effect information."""
    description: str
    probability: float  # 0-1
    severity: str  # "mild", "moderate", "severe"
    reversible: bool = True
    onset_time: Optional[str] = None  # "immediate", "days", "weeks", "months", "years"
    evidence_level: str = "theoretical"  # "theoretical", "animal_studies", "human_case", "clinical_trial"
    
    @property
    def risk_score(self) -> float:
        """Calculate overall risk score (0-10)."""
        severity_weights = {"mild": 1, "moderate": 3, "severe": 5}
        reversibility_factor = 1.0 if self.reversible else 2.0
        
        base_score = self.probability * severity_weights.get(self.severity, 3)
        return min(10.0, base_score * reversibility_factor)


@dataclass
class SafetyScore:
    """Comprehensive safety assessment."""
    overall_score: float  # 0-100
    off_target_score: float
    genotoxicity_score: float
    immunogenicity_score: float
    essential_gene_score: float
    confidence_level: float = 0.8
    
    @property
    def risk_level(self) -> RiskLevel:
        """Overall risk level assessment."""
        if self.overall_score >= 90:
            return RiskLevel.VERY_LOW
        elif self.overall_score >= 70:
            return RiskLevel.LOW
        elif self.overall_score >= 50:
            return RiskLevel.MODERATE
        elif self.overall_score >= 30:
            return RiskLevel.HIGH
        else:
            return RiskLevel.VERY_HIGH
    
    @property
    def recommendation(self) -> str:
        """Safety recommendation."""
        risk = self.risk_level
        if risk == RiskLevel.VERY_LOW:
            return "Proceed with standard monitoring"
        elif risk == RiskLevel.LOW:
            return "Proceed with enhanced monitoring"
        elif risk == RiskLevel.MODERATE:
            return "Proceed with caution and intensive monitoring"
        elif risk == RiskLevel.HIGH:
            return "Consider alternative approaches"
        else:
            return "Not recommended - significant safety concerns"


@dataclass
class VariantEffect:
    """Complete effect prediction for a genetic variant."""
    variant: VariantInfo
    protein_effect: ProteinEffect
    enhancement_gain: EnhancementGain
    side_effects: List[SideEffect] = field(default_factory=list)
    population_frequency: Optional[float] = None
    
    @property
    def benefit_risk_ratio(self) -> float:
        """Calculate benefit-to-risk ratio."""
        total_risk = sum(se.risk_score for se in self.side_effects)
        if total_risk == 0:
            return float('inf')
        return self.enhancement_gain.improvement_factor / (total_risk / 10)
    
    @property
    def overall_recommendation(self) -> str:
        """Overall recommendation based on benefit-risk analysis."""
        ratio = self.benefit_risk_ratio
        if ratio > 5:
            return "Highly recommended"
        elif ratio > 2:
            return "Recommended"
        elif ratio > 1:
            return "Consider carefully"
        elif ratio > 0.5:
            return "Not recommended"
        else:
            return "Strongly discouraged"


@dataclass
class ExperimentalValidation:
    """Experimental validation data."""
    study_type: str  # "in_vitro", "animal_model", "clinical_trial"
    organism: str
    sample_size: int
    success_rate: float
    effect_size: float
    statistical_significance: float
    publication_date: Optional[datetime] = None
    reference: Optional[str] = None


@dataclass
class EnhancementReport:
    """Comprehensive enhancement analysis report."""
    gene_name: str
    target_variant: str
    analysis_date: datetime
    best_guide: GuideRNA
    safety_assessment: SafetyScore
    predicted_effect: VariantEffect
    alternative_guides: List[GuideRNA] = field(default_factory=list)
    experimental_evidence: List[ExperimentalValidation] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)
    confidence_score: float = 0.0
    warnings: List[str] = field(default_factory=list)
    
    @property
    def feasibility_score(self) -> float:
        """Overall feasibility score (0-100)."""
        guide_score = self.best_guide.efficiency_score.overall_efficiency * 100
        safety_score = self.safety_assessment.overall_score
        effect_score = min(100, self.predicted_effect.enhancement_gain.improvement_factor * 50)
        
        # Weighted average
        return (guide_score * 0.3 + safety_score * 0.4 + effect_score * 0.3)
    
    @property
    def summary(self) -> str:
        """Generate summary text."""
        gene = self.gene_name
        variant = self.target_variant
        feasibility = self.feasibility_score
        safety = self.safety_assessment.risk_level.value
        
        return (f"Enhancement of {gene} via {variant}: "
                f"Feasibility {feasibility:.1f}/100, "
                f"Safety risk: {safety}, "
                f"Expected improvement: {self.predicted_effect.enhancement_gain.improvement_factor:.1f}x")


@dataclass
class BatchReport:
    """Results from batch analysis of multiple genes."""
    analysis_date: datetime
    total_genes: int
    successful_analyses: int
    results: Dict[str, EnhancementReport] = field(default_factory=dict)
    failed_analyses: Dict[str, str] = field(default_factory=dict)  # gene -> error message
    
    @property
    def success_rate(self) -> float:
        """Success rate of batch analysis."""
        if self.total_genes == 0:
            return 0.0
        return self.successful_analyses / self.total_genes
    
    @property
    def top_candidates(self) -> List[Tuple[str, float]]:
        """Top enhancement candidates by feasibility score."""
        candidates = [(gene, report.feasibility_score) 
                     for gene, report in self.results.items()]
        return sorted(candidates, key=lambda x: x[1], reverse=True)
    
    def get_by_category(self, category: EnhancementCategory) -> Dict[str, EnhancementReport]:
        """Filter results by enhancement category."""
        return {gene: report for gene, report in self.results.items() 
                if report.predicted_effect.enhancement_gain.category == category}


@dataclass
class ProjectConfig:
    """Configuration for an enhancement project."""
    project_name: str
    researcher_email: str
    target_organism: str = "homo sapiens"
    cas_type: CasType = CasType.CAS9
    safety_threshold: float = 70.0
    efficiency_threshold: float = 0.5
    max_off_targets: int = 5
    enable_caching: bool = True
    cache_directory: str = "data/cache"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'project_name': self.project_name,
            'researcher_email': self.researcher_email,
            'target_organism': self.target_organism,
            'cas_type': self.cas_type.value,
            'safety_threshold': self.safety_threshold,
            'efficiency_threshold': self.efficiency_threshold,
            'max_off_targets': self.max_off_targets,
            'enable_caching': self.enable_caching,
            'cache_directory': self.cache_directory
        }


# Type aliases for convenience
GeneName = str
VariantName = str
Sequence = str
Position = int
Score = float

# Collection types
GeneDatabase = Dict[GeneName, GeneInfo]
VariantDatabase = Dict[VariantName, VariantInfo]
GuideCollection = List[GuideRNA]
OffTargetCollection = List[OffTarget]