"""
Safety analysis module for Enhancement Engine.

This module provides comprehensive safety assessment for genetic enhancement
including:
- Off-target risk analysis
- Essential gene impact assessment
- Genotoxicity evaluation
- Immunogenicity prediction
- Overall risk scoring and recommendations
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any, Set
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError("BioPython is required. Install with: pip install biopython")

from ..models.data_classes import (
    GuideRNA, OffTarget, SafetyScore, RiskLevel, SideEffect,
    GeneInfo, VariantInfo, CasType, EnhancementCategory
)
from ..models.constants import (
    DEFAULT_PARAMETERS, SCORING_THRESHOLDS, ERROR_MESSAGES
)


class SafetyAnalysisError(Exception):
    """Custom exception for safety analysis errors."""
    pass


class EssentialGeneDatabase:
    """Database of essential genes for safety assessment."""
    
    
    def __init__(self):
        """Initialize essential gene database."""
        self.logger = logging.getLogger(__name__)
        
        # Load essential genes (simplified list - in practice would be comprehensive)
        self.essential_genes = {
            # DNA repair and cell cycle
            'TP53', 'BRCA1', 'BRCA2', 'ATM', 'ATR', 'CHEK1', 'CHEK2',
            'RAD51', 'MDM2', 'CDKN1A', 'RB1', 'E2F1',
            
            # Transcription and RNA processing
            'POLR2A', 'POLR2B', 'TAF1', 'TBP', 'GTF2F1', 'HNRNPA1',
            'SF3B1', 'SRSF1', 'U2AF1',
            
            # Translation and protein folding
            'RPL5', 'RPL11', 'RPS6', 'RPS14', 'EIF4E', 'EIF4G1',
            'HSP90AA1', 'HSPA8', 'HSPD1',
            
            # Metabolism
            'GAPDH', 'ENO1', 'LDHA', 'PFKL', 'ALDOA', 'TPI1',
            'PGAM1', 'PKM', 'HK1', 'HK2',
            
            # Cell structure and movement
            'ACTB', 'ACTG1', 'TUBB', 'TUBA1A', 'MYH9', 'VIM',
            
            # Signal transduction
            'AKT1', 'PIK3CA', 'MTOR', 'MAPK1', 'MAPK3', 'JUN',
            'FOS', 'MYC', 'EGFR', 'ERBB2'
        }
        
        # Essential gene categories with risk weights
        self.essential_categories = {
            'DNA_repair': {'weight': 3.0, 'genes': {'TP53', 'BRCA1', 'BRCA2', 'ATM', 'ATR'}},
            'cell_cycle': {'weight': 2.5, 'genes': {'CDKN1A', 'RB1', 'E2F1', 'CHEK1', 'CHEK2'}},
            'transcription': {'weight': 2.0, 'genes': {'POLR2A', 'POLR2B', 'TAF1', 'TBP'}},
            'translation': {'weight': 2.0, 'genes': {'RPL5', 'RPL11', 'RPS6', 'EIF4E'}},
            'metabolism': {'weight': 1.5, 'genes': {'GAPDH', 'ENO1', 'LDHA', 'PFKL'}},
            'structure': {'weight': 1.0, 'genes': {'ACTB', 'TUBB', 'TUBA1A'}},
            'signaling': {'weight': 1.5, 'genes': {'AKT1', 'PIK3CA', 'MTOR', 'MYC'}}
        }
    
    def analyze_safety(self, guide_rna: GuideRNA, target_gene: str,
                      modification_type: str = "knockout",
                      delivery_method: str = "LNP") -> SafetyScore:
        """
        Comprehensive safety analysis of a genetic modification.
        
        Args:
            guide_rna: GuideRNA object to analyze
            target_gene: Target gene symbol
            modification_type: Type of modification
            delivery_method: Delivery method
            
        Returns:
            SafetyScore object with comprehensive assessment
        """
        try:
            # Perform individual assessments
            off_target_analysis = self.off_target_analyzer.analyze_off_targets(guide_rna)
            genotoxicity_analysis = self.genotoxicity_assessor.assess_genotoxicity(
                guide_rna, target_gene, modification_type
            )
            immunogenicity_analysis = self.immunogenicity_predictor.predict_immunogenicity(
                guide_rna, guide_rna.cas_type, delivery_method
            )
            
            # Convert individual scores to 0-100 scale
            off_target_score = max(0, 100 - (off_target_analysis['risk_score'] * 10))
            genotoxicity_score = max(0, 100 - (genotoxicity_analysis['overall_score'] * 10))
            immunogenicity_score = max(0, 100 - (immunogenicity_analysis['overall_score'] * 10))
            
            # Essential gene penalty
            essential_gene_score = self._calculate_essential_gene_score(
                off_target_analysis['essential_analysis']
            )
            
            # Calculate overall safety score (weighted average)
            overall_score = self._calculate_overall_safety_score(
                off_target_score, genotoxicity_score, immunogenicity_score, essential_gene_score
            )
            
            # Determine confidence level
            confidence = self._calculate_confidence_level(
                off_target_analysis, genotoxicity_analysis, immunogenicity_analysis
            )
            
            safety_score = SafetyScore(
                overall_score=overall_score,
                off_target_score=off_target_score,
                genotoxicity_score=genotoxicity_score,
                immunogenicity_score=immunogenicity_score,
                essential_gene_score=essential_gene_score,
                confidence_level=confidence
            )
            
            self.logger.info(f"Safety analysis completed for {target_gene}: {overall_score:.1f}/100")
            return safety_score
            
        except Exception as e:
            self.logger.error(f"Safety analysis failed: {e}")
            raise SafetyAnalysisError(f"Failed to analyze safety: {e}")
    
    def _calculate_essential_gene_score(self, essential_analysis: Dict[str, Any]) -> float:
        """Calculate essential gene safety score."""
        if essential_analysis['total_essential_hits'] == 0:
            return 100.0
        
        # Start with base score
        score = 100.0
        
        # Penalty for essential gene hits
        score -= essential_analysis['total_essential_hits'] * 15
        
        # Extra penalty for critical hits
        score -= len(essential_analysis['critical_hits']) * 25
        
        # Penalty based on highest risk score
        score -= essential_analysis['highest_risk_score'] * 10
        
        return max(0.0, score)
    
    def _calculate_overall_safety_score(self, off_target: float, genotoxicity: float,
                                      immunogenicity: float, essential_gene: float) -> float:
        """Calculate weighted overall safety score."""
        weights = {
            'off_target': 0.35,
            'genotoxicity': 0.25,
            'immunogenicity': 0.20,
            'essential_gene': 0.20
        }
        
        overall = (
            off_target * weights['off_target'] +
            genotoxicity * weights['genotoxicity'] +
            immunogenicity * weights['immunogenicity'] +
            essential_gene * weights['essential_gene']
        )
        
        return max(0.0, min(100.0, overall))
    
    def _calculate_confidence_level(self, off_target_analysis: Dict,
                                  genotoxicity_analysis: Dict,
                                  immunogenicity_analysis: Dict) -> float:
        """Calculate confidence level of safety assessment."""
        base_confidence = 0.7
        
        # Higher confidence with more data
        if len(off_target_analysis.get('high_confidence_off_targets', [])) > 5:
            base_confidence += 0.1
        
        # Lower confidence for novel modifications
        if genotoxicity_analysis.get('mutagenesis', {}).get('high_risk_gene', False):
            base_confidence -= 0.2
        
        # Adjust for delivery method uncertainty
        delivery_score = immunogenicity_analysis.get('delivery_immunogenicity', {}).get('score', 0)
        if delivery_score > 5:
            base_confidence -= 0.1
        
        return max(0.1, min(1.0, base_confidence))
    
    def generate_comprehensive_safety_report(self, guide_rna: GuideRNA, target_gene: str,
                                           modification_type: str = "knockout",
                                           delivery_method: str = "LNP") -> Dict[str, Any]:
        """
        Generate comprehensive safety report.
        
        Args:
            guide_rna: GuideRNA to analyze
            target_gene: Target gene symbol
            modification_type: Type of modification
            delivery_method: Delivery method
            
        Returns:
            Comprehensive safety report dictionary
        """
        # Perform all analyses
        safety_score = self.analyze_safety(guide_rna, target_gene, modification_type, delivery_method)
        
        off_target_analysis = self.off_target_analyzer.analyze_off_targets(guide_rna)
        genotoxicity_analysis = self.genotoxicity_assessor.assess_genotoxicity(
            guide_rna, target_gene, modification_type
        )
        immunogenicity_analysis = self.immunogenicity_predictor.predict_immunogenicity(
            guide_rna, guide_rna.cas_type, delivery_method
        )
        
        # Compile all recommendations
        all_recommendations = []
        all_recommendations.extend(off_target_analysis.get('recommendations', []))
        all_recommendations.extend(genotoxicity_analysis.get('recommendations', []))
        all_recommendations.extend(immunogenicity_analysis.get('recommendations', []))
        
        # Add overall recommendations
        all_recommendations.extend(self._generate_overall_recommendations(safety_score))
        
        # Generate side effects prediction
        predicted_side_effects = self._predict_side_effects(
            off_target_analysis, genotoxicity_analysis, immunogenicity_analysis
        )
        
        return {
            'safety_score': safety_score,
            'detailed_analyses': {
                'off_target': off_target_analysis,
                'genotoxicity': genotoxicity_analysis,
                'immunogenicity': immunogenicity_analysis
            },
            'recommendations': all_recommendations,
            'predicted_side_effects': predicted_side_effects,
            'monitoring_requirements': self._generate_monitoring_requirements(safety_score),
            'risk_mitigation_strategies': self._generate_risk_mitigation_strategies(
                off_target_analysis, genotoxicity_analysis, immunogenicity_analysis
            )
        }
    
    def _generate_overall_recommendations(self, safety_score: SafetyScore) -> List[str]:
        """Generate overall safety recommendations."""
        recommendations = []
        
        if safety_score.overall_score >= 90:
            recommendations.append("Excellent safety profile - proceed with standard protocols")
        elif safety_score.overall_score >= 70:
            recommendations.append("Good safety profile - proceed with enhanced monitoring")
        elif safety_score.overall_score >= 50:
            recommendations.append("Moderate safety concerns - consider risk mitigation strategies")
        elif safety_score.overall_score >= 30:
            recommendations.append("Significant safety concerns - extensive evaluation required")
        else:
            recommendations.append("High safety risk - consider alternative approaches")
        
        # Specific score-based recommendations
        if safety_score.off_target_score < 50:
            recommendations.append("Off-target risk is primary concern - validate experimentally")
        
        if safety_score.genotoxicity_score < 50:
            recommendations.append("Genotoxicity risk elevated - monitor chromosomal stability")
        
        if safety_score.immunogenicity_score < 50:
            recommendations.append("Immunogenicity risk high - consider immunosuppression")
        
        if safety_score.essential_gene_score < 50:
            recommendations.append("Essential gene impact detected - critical safety evaluation needed")
        
        return recommendations
    
    def _predict_side_effects(self, off_target_analysis: Dict,
                            genotoxicity_analysis: Dict,
                            immunogenicity_analysis: Dict) -> List[SideEffect]:
        """Predict potential side effects based on safety analysis."""
        side_effects = []
        
        # Off-target related side effects
        if off_target_analysis['risk_score'] > 2.0:
            side_effects.append(SideEffect(
                description="Unintended genetic modifications in off-target sites",
                probability=min(0.8, off_target_analysis['risk_score'] / 5.0),
                severity="moderate",
                reversible=False,
                onset_time="immediate",
                evidence_level="theoretical"
            ))
        
        # Essential gene side effects
        if off_target_analysis['essential_analysis']['total_essential_hits'] > 0:
            side_effects.append(SideEffect(
                description="Disruption of essential cellular functions",
                probability=0.3 * off_target_analysis['essential_analysis']['total_essential_hits'],
                severity="severe",
                reversible=False,
                onset_time="days",
                evidence_level="theoretical"
            ))
        
        # Genotoxicity side effects
        if genotoxicity_analysis['overall_score'] > 5.0:
            side_effects.append(SideEffect(
                description="Chromosomal instability and increased mutation rate",
                probability=min(0.6, genotoxicity_analysis['overall_score'] / 10.0),
                severity="severe",
                reversible=False,
                onset_time="weeks",
                evidence_level="theoretical"
            ))
        
        # Immunogenicity side effects
        if immunogenicity_analysis['overall_score'] > 4.0:
            side_effects.append(SideEffect(
                description="Immune response against CRISPR components",
                probability=min(0.7, immunogenicity_analysis['overall_score'] / 8.0),
                severity="moderate",
                reversible=True,
                onset_time="days",
                evidence_level="clinical_trial"
            ))
        
        return side_effects
    
    def _generate_monitoring_requirements(self, safety_score: SafetyScore) -> List[str]:
        """Generate monitoring requirements based on safety score."""
        monitoring = []
        
        # Base monitoring
        monitoring.append("Regular clinical assessment and vital signs monitoring")
        
        # Off-target monitoring
        if safety_score.off_target_score < 70:
            monitoring.append("Periodic whole-genome sequencing to detect off-target effects")
            monitoring.append("Targeted sequencing of predicted off-target sites")
        
        # Genotoxicity monitoring
        if safety_score.genotoxicity_score < 70:
            monitoring.append("Chromosomal stability assessment via karyotyping")
            monitoring.append("Micronucleus assay for DNA damage detection")
            monitoring.append("Long-term cancer surveillance")
        
        # Immunogenicity monitoring
        if safety_score.immunogenicity_score < 70:
            monitoring.append("Anti-Cas protein antibody monitoring")
            monitoring.append("T-cell response assessment")
            monitoring.append("Inflammatory marker surveillance")
        
        # Essential gene monitoring
        if safety_score.essential_gene_score < 70:
            monitoring.append("Functional assessment of essential cellular pathways")
            monitoring.append("Multi-organ system monitoring")
        
        # Risk-based monitoring frequency
        if safety_score.overall_score < 50:
            monitoring.append("Weekly monitoring for first month, then monthly for 6 months")
        elif safety_score.overall_score < 70:
            monitoring.append("Bi-weekly monitoring for first month, then monthly for 3 months")
        else:
            monitoring.append("Standard monitoring schedule with safety checkpoints")
        
        return monitoring
    
    def _generate_risk_mitigation_strategies(self, off_target_analysis: Dict,
                                           genotoxicity_analysis: Dict,
                                           immunogenicity_analysis: Dict) -> List[str]:
        """Generate risk mitigation strategies."""
        strategies = []
        
        # Off-target mitigation
        if off_target_analysis['risk_score'] > 2.0:
            strategies.append("Use high-fidelity Cas variants (e.g., SpCas9-HF1)")
            strategies.append("Optimize guide RNA design for improved specificity")
            strategies.append("Use truncated guide RNAs (17-18 nucleotides)")
            strategies.append("Implement dual-guide approaches for increased specificity")
        
        # Genotoxicity mitigation
        if genotoxicity_analysis['overall_score'] > 5.0:
            strategies.append("Use base editors or prime editors instead of nucleases")
            strategies.append("Minimize Cas protein expression time")
            strategies.append("Co-deliver DNA repair pathway enhancers")
        
        # Immunogenicity mitigation
        if immunogenicity_analysis['overall_score'] > 4.0:
            strategies.append("Use immunosuppressive protocols during treatment")
            strategies.append("Consider humanized Cas proteins")
            strategies.append("Use alternative delivery methods with lower immunogenicity")
            strategies.append("Pre-screen for pre-existing immunity to Cas proteins")
        
        # General mitigation
        strategies.append("Establish clear stopping criteria for adverse events")
        strategies.append("Develop reversal strategies where possible")
        strategies.append("Maintain comprehensive adverse event reporting")
        
        return strategies
    
    def compare_guide_safety(self, guides: List[GuideRNA], target_gene: str) -> List[Tuple[GuideRNA, SafetyScore]]:
        """
        Compare safety of multiple guide RNAs.
        
        Args:
            guides: List of GuideRNA objects to compare
            target_gene: Target gene symbol
            
        Returns:
            List of (GuideRNA, SafetyScore) tuples sorted by safety
        """
        guide_safety_pairs = []
        
        for guide in guides:
            try:
                safety_score = self.analyze_safety(guide, target_gene)
                guide_safety_pairs.append((guide, safety_score))
            except Exception as e:
                self.logger.warning(f"Failed to analyze safety for guide {guide.sequence}: {e}")
        
        # Sort by overall safety score (descending)
        guide_safety_pairs.sort(key=lambda x: x[1].overall_score, reverse=True)
        
        return guide_safety_pairs
    
    def assess_enhancement_safety(self, target_gene: str, enhancement_category: EnhancementCategory,
                                modification_type: str = "base_edit") -> Dict[str, Any]:
        """
        Assess safety considerations specific to enhancement applications.
        
        Args:
            target_gene: Gene to be enhanced
            enhancement_category: Type of enhancement
            modification_type: Type of modification
            
        Returns:
            Enhancement-specific safety assessment
        """
        # Enhancement-specific risk factors
        enhancement_risks = {
            EnhancementCategory.COGNITIVE: {
                'brain_specific_risks': True,
                'developmental_concerns': True,
                'behavioral_effects': True
            },
            EnhancementCategory.PHYSICAL: {
                'muscle_metabolism_effects': True,
                'cardiovascular_impact': True,
                'bone_density_changes': True
            },
            EnhancementCategory.LONGEVITY: {
                'cancer_risk_modulation': True,
                'metabolic_disruption': True,
                'immune_system_effects': True
            }
        }
        
        category_risks = enhancement_risks.get(enhancement_category, {})
        
        # Generate enhancement-specific recommendations
        recommendations = []
        
        if enhancement_category == EnhancementCategory.COGNITIVE:
            recommendations.extend([
                "Monitor neurological function and cognitive assessments",
                "Assess for psychiatric side effects",
                "Consider developmental stage-specific risks"
            ])
        elif enhancement_category == EnhancementCategory.PHYSICAL:
            recommendations.extend([
                "Monitor cardiovascular function during enhanced performance",
                "Assess musculoskeletal system for overuse injuries",
                "Monitor metabolic parameters"
            ])
        elif enhancement_category == EnhancementCategory.LONGEVITY:
            recommendations.extend([
                "Long-term cancer surveillance protocols",
                "Monitor for accelerated aging in other systems",
                "Assess immune system function regularly"
            ])
        
        return {
            'enhancement_category': enhancement_category.value,
            'category_specific_risks': category_risks,
            'recommendations': recommendations,
            'monitoring_duration': self._get_enhancement_monitoring_duration(enhancement_category),
            'reversibility_assessment': self._assess_enhancement_reversibility(
                target_gene, modification_type
            )
        }
    
    def _get_enhancement_monitoring_duration(self, category: EnhancementCategory) -> str:
        """Get recommended monitoring duration for enhancement category."""
        durations = {
            EnhancementCategory.COGNITIVE: "Lifelong (neurological changes)",
            EnhancementCategory.PHYSICAL: "5-10 years (performance monitoring)",
            EnhancementCategory.LONGEVITY: "Lifelong (longevity assessment)",
            EnhancementCategory.SENSORY: "10-20 years (sensory function)",
            EnhancementCategory.METABOLIC: "Lifelong (metabolic stability)",
            EnhancementCategory.IMMUNE: "Lifelong (immune competence)"
        }
        
        return durations.get(category, "Long-term monitoring recommended")
    
    def _assess_enhancement_reversibility(self, target_gene: str, modification_type: str) -> Dict[str, Any]:
        """Assess reversibility of enhancement modifications."""
        reversibility_scores = {
            'base_edit': 0.2,     # Single nucleotide changes - difficult to reverse
            'knockout': 0.1,      # Gene loss - very difficult to reverse
            'knockin': 0.3,       # Integration - potentially reversible
            'activation': 0.8,    # Epigenetic - more reversible
            'interference': 0.9   # RNA-based - highly reversible
        }
        
        base_reversibility = reversibility_scores.get(modification_type, 0.5)
        
        # Gene-specific factors
        gene_factors = 1.0
        
        # Some genes are harder to reverse due to developmental roles
        developmental_genes = {'FOXO3', 'MYC', 'TP53'}
        if target_gene.upper() in developmental_genes:
            gene_factors *= 0.5
        
        overall_reversibility = base_reversibility * gene_factors
        
        return {
            'reversibility_score': overall_reversibility,
            'modification_factor': base_reversibility,
            'gene_factor': gene_factors,
            'reversibility_level': 'High' if overall_reversibility > 0.7 else
                                 'Moderate' if overall_reversibility > 0.4 else 'Low',
            'reversal_strategies': self._get_reversal_strategies(modification_type)
        }
    
    def _get_reversal_strategies(self, modification_type: str) -> List[str]:
        """Get potential reversal strategies for modification type."""
        strategies = {
            'base_edit': [
                "Counter base editing to revert nucleotide change",
                "Gene therapy with wild-type gene copy"
            ],
            'knockout': [
                "Gene replacement therapy",
                "Complementation with functional gene copy"
            ],
            'knockin': [
                "Targeted excision of inserted sequence",
                "Recombinase-mediated removal"
            ],
            'activation': [
                "Remove activating elements",
                "Counter-repressions"
            ],
            'interference': [
                "Discontinue interfering RNA",
                "Natural RNA degradation"
            ]
        }
        
        return strategies.get(modification_type, ["Consult genetic counselor for reversal options"])
    
    def is_essential(self, gene_symbol: str) -> bool:
        """Check if a gene is essential."""
        return gene_symbol.upper() in self.essential_genes
    
    def get_essentiality_score(self, gene_symbol: str) -> float:
        """Get essentiality risk score for a gene."""
        gene_upper = gene_symbol.upper()
        
        if not self.is_essential(gene_upper):
            return 0.0
        
        # Find category and return weight
        for category, info in self.essential_categories.items():
            if gene_upper in info['genes']:
                return info['weight']
        
        return 1.0  # Default for essential genes not in specific categories
    
    def get_gene_function_category(self, gene_symbol: str) -> Optional[str]:
        """Get functional category of an essential gene."""
        gene_upper = gene_symbol.upper()
        
        for category, info in self.essential_categories.items():
            if gene_upper in info['genes']:
                return category
        
        return None


class OffTargetAnalyzer:
    """Analyzes off-target effects and their risks."""
    
    def __init__(self):
        """Initialize off-target analyzer."""
        self.essential_db = EssentialGeneDatabase()
        self.logger = logging.getLogger(__name__)
    
    def analyze_off_targets(self, guide_rna: GuideRNA) -> Dict[str, Any]:
        """
        Comprehensive analysis of off-target effects.
        
        Args:
            guide_rna: GuideRNA object with off-targets
            
        Returns:
            Dictionary with off-target analysis results
        """
        off_targets = guide_rna.off_targets
        
        if not off_targets:
            return {
                'total_off_targets': 0,
                'high_risk_count': 0,
                'essential_gene_hits': 0,
                'risk_score': 0.0,
                'risk_level': RiskLevel.VERY_LOW,
                'recommendations': ['No off-targets detected - excellent specificity']
            }
        
        # Categorize off-targets by risk
        risk_categories = self._categorize_off_targets(off_targets)
        
        # Analyze essential gene impacts
        essential_analysis = self._analyze_essential_gene_impacts(off_targets)
        
        # Calculate overall risk score
        risk_score = self._calculate_off_target_risk_score(off_targets, essential_analysis)
        
        # Generate recommendations
        recommendations = self._generate_off_target_recommendations(
            risk_categories, essential_analysis, risk_score
        )
        
        return {
            'total_off_targets': len(off_targets),
            'risk_categories': risk_categories,
            'essential_analysis': essential_analysis,
            'risk_score': risk_score,
            'risk_level': self._score_to_risk_level(risk_score),
            'recommendations': recommendations,
            'high_confidence_off_targets': [ot for ot in off_targets if ot.score > 0.5],
            'chromosome_distribution': self._analyze_chromosome_distribution(off_targets)
        }
    
    def _categorize_off_targets(self, off_targets: List[OffTarget]) -> Dict[str, int]:
        """Categorize off-targets by risk level."""
        categories = {
            'very_high': 0,  # 0 mismatches
            'high': 0,       # 1 mismatch
            'moderate': 0,   # 2 mismatches
            'low': 0,        # 3 mismatches
            'very_low': 0    # 4+ mismatches
        }
        
        for ot in off_targets:
            if ot.mismatches == 0:
                categories['very_high'] += 1
            elif ot.mismatches == 1:
                categories['high'] += 1
            elif ot.mismatches == 2:
                categories['moderate'] += 1
            elif ot.mismatches == 3:
                categories['low'] += 1
            else:
                categories['very_low'] += 1
        
        return categories
    
    def _analyze_essential_gene_impacts(self, off_targets: List[OffTarget]) -> Dict[str, Any]:
        """Analyze impacts on essential genes."""
        essential_hits = [ot for ot in off_targets if ot.essential_gene]
        
        if not essential_hits:
            return {
                'total_essential_hits': 0,
                'by_category': {},
                'highest_risk_score': 0.0,
                'critical_hits': []
            }
        
        # Group by functional category
        by_category = defaultdict(list)
        for hit in essential_hits:
            if hit.gene_context:
                category = self.essential_db.get_gene_function_category(hit.gene_context)
                if category:
                    by_category[category].append(hit)
        
        # Find critical hits (low mismatches + essential)
        critical_hits = [
            hit for hit in essential_hits 
            if hit.mismatches <= 1 and self.essential_db.get_essentiality_score(
                hit.gene_context or "") >= 2.0
        ]
        
        # Calculate highest risk score
        risk_scores = [
            (1.0 / (1.0 + hit.mismatches)) * self.essential_db.get_essentiality_score(
                hit.gene_context or "")
            for hit in essential_hits if hit.gene_context
        ]
        
        return {
            'total_essential_hits': len(essential_hits),
            'by_category': dict(by_category),
            'highest_risk_score': max(risk_scores) if risk_scores else 0.0,
            'critical_hits': critical_hits
        }
    
    def _calculate_off_target_risk_score(self, off_targets: List[OffTarget], 
                                       essential_analysis: Dict[str, Any]) -> float:
        """Calculate overall off-target risk score (0-10)."""
        if not off_targets:
            return 0.0
        
        total_risk = 0.0
        
        for ot in off_targets:
            # Base risk from mismatches
            if ot.mismatches == 0:
                base_risk = 5.0
            elif ot.mismatches == 1:
                base_risk = 2.0
            elif ot.mismatches == 2:
                base_risk = 0.5
            else:
                base_risk = 0.1
            
            # Multiply by off-target score
            base_risk *= ot.score
            
            # Essential gene penalty
            if ot.essential_gene and ot.gene_context:
                essentiality_score = self.essential_db.get_essentiality_score(ot.gene_context)
                base_risk *= (1.0 + essentiality_score)
            
            total_risk += base_risk
        
        # Normalize and cap at 10
        return min(10.0, total_risk)
    
    def _score_to_risk_level(self, score: float) -> RiskLevel:
        """Convert risk score to risk level."""
        if score < 0.5:
            return RiskLevel.VERY_LOW
        elif score < 1.5:
            return RiskLevel.LOW
        elif score < 3.0:
            return RiskLevel.MODERATE
        elif score < 6.0:
            return RiskLevel.HIGH
        else:
            return RiskLevel.VERY_HIGH
    
    def _generate_off_target_recommendations(self, risk_categories: Dict[str, int],
                                           essential_analysis: Dict[str, Any],
                                           risk_score: float) -> List[str]:
        """Generate safety recommendations based on off-target analysis."""
        recommendations = []
        
        if risk_score < 0.5:
            recommendations.append("Excellent specificity - proceed with confidence")
        elif risk_score < 1.5:
            recommendations.append("Good specificity - standard monitoring recommended")
        elif risk_score < 3.0:
            recommendations.append("Moderate off-target risk - enhanced monitoring required")
            recommendations.append("Consider experimental validation of predicted off-targets")
        else:
            recommendations.append("High off-target risk - consider alternative guides")
            recommendations.append("Experimental off-target validation strongly recommended")
        
        # Essential gene specific recommendations
        if essential_analysis['total_essential_hits'] > 0:
            recommendations.append(
                f"Warning: {essential_analysis['total_essential_hits']} potential "
                f"essential gene hit(s) detected"
            )
            
            if essential_analysis['critical_hits']:
                recommendations.append(
                    "Critical: High-confidence off-targets in essential genes detected"
                )
                recommendations.append("Consider using a different guide RNA")
        
        # High-risk category recommendations
        if risk_categories['very_high'] > 0:
            recommendations.append(
                f"Warning: {risk_categories['very_high']} perfect match off-target(s) found"
            )
        
        if risk_categories['high'] > 3:
            recommendations.append("Consider guides with fewer single-mismatch off-targets")
        
        return recommendations
    
    def _analyze_chromosome_distribution(self, off_targets: List[OffTarget]) -> Dict[str, int]:
        """Analyze distribution of off-targets across chromosomes."""
        distribution = defaultdict(int)
        
        for ot in off_targets:
            distribution[ot.chromosome] += 1
        
        return dict(distribution)


class GenotoxicityAssessor:
    """Assesses genotoxicity risks of genetic modifications."""
    
    def __init__(self):
        """Initialize genotoxicity assessor."""
        self.logger = logging.getLogger(__name__)
    
    def assess_genotoxicity(self, guide_rna: GuideRNA, target_gene: str,
                          modification_type: str = "knockout") -> Dict[str, Any]:
        """
        Assess genotoxicity risks of a genetic modification.
        
        Args:
            guide_rna: GuideRNA object
            target_gene: Target gene symbol
            modification_type: Type of modification (knockout, knockin, base_edit)
            
        Returns:
            Genotoxicity assessment results
        """
        # Assess different genotoxicity factors
        dna_damage_risk = self._assess_dna_damage_risk(guide_rna)
        chromosomal_instability_risk = self._assess_chromosomal_instability(guide_rna)
        mutagenesis_risk = self._assess_mutagenesis_risk(modification_type, target_gene)
        
        # Calculate overall genotoxicity score
        overall_score = self._calculate_genotoxicity_score(
            dna_damage_risk, chromosomal_instability_risk, mutagenesis_risk
        )
        
        # Generate recommendations
        recommendations = self._generate_genotoxicity_recommendations(
            overall_score, dna_damage_risk, chromosomal_instability_risk, mutagenesis_risk
        )
        
        return {
            'overall_score': overall_score,
            'dna_damage_risk': dna_damage_risk,
            'chromosomal_instability_risk': chromosomal_instability_risk,
            'mutagenesis_risk': mutagenesis_risk,
            'risk_level': self._score_to_risk_level(overall_score),
            'recommendations': recommendations
        }
    
    def _assess_dna_damage_risk(self, guide_rna: GuideRNA) -> Dict[str, Any]:
        """Assess risk of DNA damage from CRISPR cutting."""
        # Higher risk with more off-targets (multiple cuts)
        off_target_count = len(guide_rna.off_targets)
        
        # Base risk from double-strand breaks
        base_risk = 1.0
        
        # Risk increases with off-targets
        off_target_risk = min(3.0, off_target_count * 0.5)
        
        # Cas type specific risks
        cas_risk_modifiers = {
            CasType.CAS9: 1.0,      # Standard DSB
            CasType.CAS12A: 0.8,    # Potentially less toxic
            CasType.BASE_EDITOR: 0.3,  # No DSB
            CasType.PRIME_EDITOR: 0.4  # Minimal DSB
        }
        
        cas_modifier = cas_risk_modifiers.get(guide_rna.cas_type, 1.0)
        
        total_risk = (base_risk + off_target_risk) * cas_modifier
        
        return {
            'score': min(5.0, total_risk),
            'base_risk': base_risk,
            'off_target_contribution': off_target_risk,
            'cas_modifier': cas_modifier
        }
    
    def _assess_chromosomal_instability(self, guide_rna: GuideRNA) -> Dict[str, Any]:
        """Assess risk of chromosomal instability."""
        # Risk from multiple cuts on same chromosome
        chromosome_hits = defaultdict(int)
        
        for ot in guide_rna.off_targets:
            chromosome_hits[ot.chromosome] += 1
        
        # Calculate instability risk
        instability_risk = 0.0
        
        for chrom, hits in chromosome_hits.items():
            if hits > 1:
                # Multiple hits on same chromosome increase instability risk
                instability_risk += (hits - 1) * 0.5
        
        # Risk from cuts near heterochromatin or repetitive regions
        repetitive_risk = sum(0.2 for ot in guide_rna.off_targets 
                            if "repeat" in (ot.gene_context or "").lower())
        
        total_risk = instability_risk + repetitive_risk
        
        return {
            'score': min(5.0, total_risk),
            'multi_hit_chromosomes': {k: v for k, v in chromosome_hits.items() if v > 1},
            'repetitive_hits': repetitive_risk / 0.2 if repetitive_risk > 0 else 0
        }
    
    def _assess_mutagenesis_risk(self, modification_type: str, target_gene: str) -> Dict[str, Any]:
        """Assess mutagenesis risk from the specific modification."""
        # Base risks by modification type
        modification_risks = {
            'knockout': 2.0,      # Complete gene loss
            'knockin': 1.5,       # Foreign DNA integration
            'base_edit': 0.5,     # Single nucleotide change
            'prime_edit': 0.8,    # Targeted insertion/deletion
            'activation': 0.3,    # No DNA change
            'interference': 0.2   # No DNA change
        }
        
        base_risk = modification_risks.get(modification_type, 1.0)
        
        # Gene-specific risk modifiers
        gene_risk_modifier = 1.0
        
        # Higher risk for certain gene categories
        high_risk_genes = {
            'TP53', 'BRCA1', 'BRCA2', 'APC', 'MLH1', 'MSH2', 'MSH6', 'PMS2'  # Tumor suppressors
        }
        
        if target_gene.upper() in high_risk_genes:
            gene_risk_modifier = 2.0
        
        total_risk = base_risk * gene_risk_modifier
        
        return {
            'score': min(5.0, total_risk),
            'modification_base_risk': base_risk,
            'gene_risk_modifier': gene_risk_modifier,
            'high_risk_gene': target_gene.upper() in high_risk_genes
        }
    
    def _calculate_genotoxicity_score(self, dna_damage: Dict, chromosomal: Dict, 
                                    mutagenesis: Dict) -> float:
        """Calculate overall genotoxicity score."""
        # Weighted combination of risk factors
        weights = {
            'dna_damage': 0.4,
            'chromosomal': 0.3,
            'mutagenesis': 0.3
        }
        
        total_score = (
            dna_damage['score'] * weights['dna_damage'] +
            chromosomal['score'] * weights['chromosomal'] +
            mutagenesis['score'] * weights['mutagenesis']
        )
        
        return min(10.0, total_score)
    
    def _score_to_risk_level(self, score: float) -> RiskLevel:
        """Convert genotoxicity score to risk level."""
        if score < 1.0:
            return RiskLevel.VERY_LOW
        elif score < 2.5:
            return RiskLevel.LOW
        elif score < 5.0:
            return RiskLevel.MODERATE
        elif score < 7.5:
            return RiskLevel.HIGH
        else:
            return RiskLevel.VERY_HIGH
    
    def _generate_genotoxicity_recommendations(self, overall_score: float,
                                             dna_damage: Dict, chromosomal: Dict,
                                             mutagenesis: Dict) -> List[str]:
        """Generate genotoxicity-specific recommendations."""
        recommendations = []
        
        if overall_score < 1.0:
            recommendations.append("Low genotoxicity risk - standard monitoring sufficient")
        elif overall_score < 2.5:
            recommendations.append("Moderate genotoxicity risk - regular safety monitoring")
        elif overall_score < 5.0:
            recommendations.append("Elevated genotoxicity risk - enhanced monitoring required")
            recommendations.append("Consider periodic genomic stability assessments")
        else:
            recommendations.append("High genotoxicity risk - extensive safety evaluation needed")
            recommendations.append("Long-term genomic stability monitoring essential")
        
        # Specific recommendations based on risk factors
        if dna_damage['score'] > 3.0:
            recommendations.append("High DNA damage risk due to multiple off-targets")
            recommendations.append("Consider using high-fidelity Cas variants")
        
        if chromosomal['score'] > 2.0:
            recommendations.append("Chromosomal instability risk detected")
            recommendations.append("Monitor for large chromosomal rearrangements")
        
        if mutagenesis['high_risk_gene']:
            recommendations.append("Targeting tumor suppressor gene - cancer risk assessment needed")
        
        return recommendations


class ImmunogenicityPredictor:
    """Predicts immunogenic responses to CRISPR components."""
    
    def __init__(self):
        """Initialize immunogenicity predictor."""
        self.logger = logging.getLogger(__name__)
    
    def predict_immunogenicity(self, guide_rna: GuideRNA, cas_type: CasType,
                             delivery_method: str = "LNP") -> Dict[str, Any]:
        """
        Predict immunogenic responses.
        
        Args:
            guide_rna: GuideRNA object
            cas_type: Cas protein type
            delivery_method: Delivery method (LNP, AAV, etc.)
            
        Returns:
            Immunogenicity prediction results
        """
        # Assess different immunogenicity factors
        cas_immunogenicity = self._assess_cas_immunogenicity(cas_type)
        guide_immunogenicity = self._assess_guide_immunogenicity(guide_rna)
        delivery_immunogenicity = self._assess_delivery_immunogenicity(delivery_method)
        
        # Calculate overall immunogenicity score
        overall_score = self._calculate_immunogenicity_score(
            cas_immunogenicity, guide_immunogenicity, delivery_immunogenicity
        )
        
        # Generate recommendations
        recommendations = self._generate_immunogenicity_recommendations(
            overall_score, cas_immunogenicity, delivery_immunogenicity
        )
        
        return {
            'overall_score': overall_score,
            'cas_immunogenicity': cas_immunogenicity,
            'guide_immunogenicity': guide_immunogenicity,
            'delivery_immunogenicity': delivery_immunogenicity,
            'risk_level': self._score_to_risk_level(overall_score),
            'recommendations': recommendations
        }
    
    def _assess_cas_immunogenicity(self, cas_type: CasType) -> Dict[str, Any]:
        """Assess immunogenicity of Cas proteins."""
        # Base immunogenicity scores for different Cas types
        cas_scores = {
            CasType.CAS9: 3.0,      # Well-studied, moderate immunogenicity
            CasType.CAS12A: 2.5,    # Less immunogenic than Cas9
            CasType.BASE_EDITOR: 3.5,  # Larger protein complex
            CasType.PRIME_EDITOR: 4.0   # Very large protein complex
        }
        
        base_score = cas_scores.get(cas_type, 3.0)
        
        # Factors affecting immunogenicity
        factors = {
            'protein_size': 1.0,    # Larger proteins more immunogenic
            'bacterial_origin': 0.5,  # Bacterial proteins are foreign
            'expression_level': 0.8   # Higher expression = more immune response
        }
        
        # Adjust based on Cas type
        if cas_type in [CasType.BASE_EDITOR, CasType.PRIME_EDITOR]:
            factors['protein_size'] = 1.5  # Larger fusion proteins
        
        total_score = base_score + sum(factors.values())
        
        return {
            'score': min(10.0, total_score),
            'base_immunogenicity': base_score,
            'modifying_factors': factors
        }
    
    def _assess_guide_immunogenicity(self, guide_rna: GuideRNA) -> Dict[str, Any]:
        """Assess immunogenicity of guide RNA."""
        # Guide RNAs generally have low immunogenicity
        base_score = 0.5
        
        # Factors that might increase immunogenicity
        immunogenic_motifs = 0
        guide_seq = guide_rna.sequence.upper()
        
        # Look for potentially immunogenic sequences
        if 'GGGG' in guide_seq:  # G-quadruplex forming
            immunogenic_motifs += 1
        
        if guide_seq.count('CG') > 3:  # CpG motifs
            immunogenic_motifs += 1
        
        # High GC content can be immunogenic
        if guide_rna.gc_content > 70:
            immunogenic_motifs += 1
        
        total_score = base_score + (immunogenic_motifs * 0.3)
        
        return {
            'score': min(5.0, total_score),
            'immunogenic_motifs': immunogenic_motifs,
            'gc_content_risk': guide_rna.gc_content > 70
        }
    
    def _assess_delivery_immunogenicity(self, delivery_method: str) -> Dict[str, Any]:
        """Assess immunogenicity of delivery method."""
        # Immunogenicity scores for different delivery methods
        delivery_scores = {
            'LNP': 2.0,           # Lipid nanoparticles - moderate
            'AAV': 4.0,           # Adeno-associated virus - higher
            'lentivirus': 5.0,    # Lentivirus - high
            'electroporation': 1.0,  # Direct delivery - low
            'microinjection': 0.5,   # Direct injection - very low
            'transfection': 1.5   # Chemical transfection - low-moderate
        }
        
        base_score = delivery_scores.get(delivery_method.lower(), 2.0)
        
        # Additional factors
        factors = {}
        
        if 'virus' in delivery_method.lower():
            factors['viral_components'] = 1.0
        
        if delivery_method.upper() == 'AAV':
            factors['capsid_immunogenicity'] = 0.5
        
        total_score = base_score + sum(factors.values())
        
        return {
            'score': min(10.0, total_score),
            'base_score': base_score,
            'additional_factors': factors
        }
    
    def _calculate_immunogenicity_score(self, cas_immun: Dict, guide_immun: Dict,
                                      delivery_immun: Dict) -> float:
        """Calculate overall immunogenicity score."""
        # Weighted combination
        weights = {
            'cas': 0.5,      # Cas protein is major component
            'guide': 0.1,    # Guide RNA has low immunogenicity
            'delivery': 0.4  # Delivery method is important
        }
        
        total_score = (
            cas_immun['score'] * weights['cas'] +
            guide_immun['score'] * weights['guide'] +
            delivery_immun['score'] * weights['delivery']
        )
        
        return min(10.0, total_score)
    
    def _score_to_risk_level(self, score: float) -> RiskLevel:
        """Convert immunogenicity score to risk level."""
        if score < 2.0:
            return RiskLevel.VERY_LOW
        elif score < 4.0:
            return RiskLevel.LOW
        elif score < 6.0:
            return RiskLevel.MODERATE
        elif score < 8.0:
            return RiskLevel.HIGH
        else:
            return RiskLevel.VERY_HIGH
    
    def _generate_immunogenicity_recommendations(self, overall_score: float,
                                               cas_immun: Dict, delivery_immun: Dict) -> List[str]:
        """Generate immunogenicity-specific recommendations."""
        recommendations = []
        
        if overall_score < 2.0:
            recommendations.append("Low immunogenicity risk - standard monitoring")
        elif overall_score < 4.0:
            recommendations.append("Moderate immunogenicity - monitor for immune responses")
        elif overall_score < 6.0:
            recommendations.append("Elevated immunogenicity - immunosuppression may be needed")
        else:
            recommendations.append("High immunogenicity risk - consider alternative approaches")
        
        # Specific recommendations
        if cas_immun['score'] > 5.0:
            recommendations.append("Consider using smaller Cas variants or base editors")
        
        if delivery_immun['score'] > 5.0:
            recommendations.append("Consider alternative delivery methods")
            if 'virus' in str(delivery_immun.get('additional_factors', {})):
                recommendations.append("Viral delivery detected - monitor for anti-viral responses")
        
        return recommendations


class SafetyAnalyzer:
    """Main safety analysis class integrating all safety assessments."""
    
    def __init__(self):
        """Initialize safety analyzer."""
        self.off_target_analyzer = OffTargetAnalyzer()
        self.genotoxicity_assessor = GenotoxicityAssessor()
        self.immunogenicity_predictor = ImmunogenicityPredictor()
        
        self.logger = logging.getLogger(__name__)
    
    def analyze_safety(self, guide_rna: GuideRNA, target_gene: str,
                      modification_type: str = "knockout",
                      delivery_method: str = "LNP") -> SafetyScore:
        """
        Comprehensive safety analysis of a genetic modification.
        
        Args:
            guide_rna: GuideRNA object to analyze
            target_gene: Target gene symbol
            modification_type: Type of modification
            delivery_method: Delivery method
            
        Returns:
            SafetyScore object with comprehensive assessment
        """
        try:
            # Perform individual assessments
            off_target_analysis = self.off_target_analyzer.analyze_off_targets(guide_rna)
            genotoxicity_analysis = self.genotoxicity_assessor.assess_genotoxicity(
                guide_rna, target_gene, modification_type
            )
            immunogenicity_analysis = self.immunogenicity_predictor.predict_immunogenicity(
                guide_rna, guide_rna.cas_type, delivery_method
            )
            
            # Convert individual scores to 0-100 scale
            off_target_score = max(0, 100 - (off_target_analysis['risk_score'] * 10))
            genotoxicity_score = max(0, 100 - (genotoxicity_analysis['overall_score'] * 10))
            immunogenicity_score = max(0, 100 - (immunogenicity_analysis['overall_score'] * 10))
            
            # Essential gene penalty
            essential_gene_score = self._calculate_essential_gene_score(
                off_target_analysis['essential_analysis']
            )
            
            # Calculate overall safety score (weighted average)
            overall_score = self._calculate_overall_safety_score(
                off_target_score, genotoxicity_score, immunogenicity_score, essential_gene_score
            )
            
            # Determine confidence level
            confidence = self._calculate_confidence_level(
                off_target_analysis, genotoxicity_analysis, immunogenicity_analysis
            )
            
            safety_score = SafetyScore(
                overall_score=overall_score,
                off_target_score=off_target_score,
                genotoxicity_score=genotoxicity_score,
                immunogenicity_score=immunogenicity_score,
                essential_gene_score=essential_gene_score,
                confidence_level=confidence
            )
            
            self.logger.info(f"Safety analysis completed for {target_gene}: {overall_score:.1f}/100")
            return safety_score
            
        except Exception as e:
            self.logger.error(f"Safety analysis failed: {e}")
            raise SafetyAnalysisError(f"Failed to analyze safety: {e}")
    
    def _calculate_essential_gene_score(self, essential_analysis: Dict[str, Any]) -> float:
        """Calculate essential gene safety score."""
        if essential_analysis['total_essential_hits'] == 0:
            return 100.0
        
        # Start with base score
        score = 100.0
        
        # Penalty for essential gene hits
        score -= essential_analysis['total_essential_hits'] * 15
        
        # Extra penalty for critical hits
        score -= len(essential_analysis['critical_hits']) * 25
        
        # Penalty based on highest risk score
        score -= essential_analysis['highest_risk_score'] * 10
        
        return max(0.0, score)
    
    def _calculate_overall_safety_score(self, off_target: float, genotoxicity: float,
                                      immunogenicity: float, essential_gene: float) -> float:
        """Calculate weighted overall safety score."""
        weights = {
            'off_target': 0.35,
            'genotoxicity': 0.25,
            'immunogenicity': 0.20,
            'essential_gene': 0.20
        }
        
        overall = (
            off_target * weights['off_target'] +
            genotoxicity * weights['genotoxicity'] +
            immunogenicity * weights['immunogenicity'] +
            essential_gene * weights['essential_gene']
        )
        
        return max(0.0, min(100.0, overall))
    
    def _calculate_confidence_level(self, off_target_analysis: Dict,
                                  genotoxicity_analysis: Dict,
                                  immunogenicity_analysis: Dict) -> float:
        """Calculate confidence level of safety assessment."""
        base_confidence = 0.7
        
        # Higher confidence with more data
        if len(off_target_analysis.get('high_confidence_off_targets', [])) > 5:
            base_confidence += 0.1
        
        # Lower confidence for novel modifications
        if genotoxicity_analysis.get('mutagenesis', {}).get('high_risk_gene', False):
            base_confidence -= 0.2
        
        # Adjust for delivery method uncertainty
        delivery_score = immunogenicity_analysis.get('delivery_immunogenicity', {}).get('score', 0)
        if delivery_score > 5:
            base_confidence -= 0.1
        
        return max(0.1, min(1.0, base_confidence))
    
    def generate_comprehensive_safety_report(self, guide_rna: GuideRNA, target_gene: str,
                                           modification_type: str = "knockout",
                                           delivery_method: str = "LNP") -> Dict[str, Any]:
        """
        Generate comprehensive safety report.
        
        Args:
            guide_rna: GuideRNA to analyze
            target_gene: Target gene symbol
            modification_type: Type of modification
            delivery_method: Delivery method
            
        Returns:
            Comprehensive safety report dictionary
        """
        # Perform all analyses
        safety_score = self.analyze_safety(guide_rna, target_gene, modification_type, delivery_method)
        
        off_target_analysis = self.off_target_analyzer.analyze_off_targets(guide_rna)
        genotoxicity_analysis = self.genotoxicity_assessor.assess_genotoxicity(
            guide_rna, target_gene, modification_type
        )
        immunogenicity_analysis = self.immunogenicity_predictor.predict_immunogenicity(
            guide_rna, guide_rna.cas_type, delivery_method
        )
        
        # Compile all recommendations
        all_recommendations = []
        all_recommendations.extend(off_target_analysis.get('recommendations', []))
        all_recommendations.extend(genotoxicity_analysis.get('recommendations', []))
        all_recommendations.extend(immunogenicity_analysis.get('recommendations', []))
        
        # Add overall recommendations
        all_recommendations.extend(self._generate_overall_recommendations(safety_score))
        
        # Generate side effects prediction
        predicted_side_effects = self._predict_side_effects(
            off_target_analysis, genotoxicity_analysis, immunogenicity_analysis
        )
        
        return {
            'safety_score': safety_score,
            'detailed_analyses': {
                'off_target': off_target_analysis,
                'genotoxicity': genotoxicity_analysis,
                'immunogenicity': immunogenicity_analysis
            },
            'recommendations': all_recommendations,
            'predicted_side_effects': predicted_side_effects,
            'monitoring_requirements': self._generate_monitoring_requirements(safety_score),
            'risk_mitigation_strategies': self._generate_risk_mitigation_strategies(
                off_target_analysis, genotoxicity_analysis, immunogenicity_analysis
            )
        }
    
    def _generate_overall_recommendations(self, safety_score: SafetyScore) -> List[str]:
        """Generate overall safety recommendations."""
        recommendations = []
        
        if safety_score.overall_score >= 90:
            recommendations.append("Excellent safety profile - proceed with standard protocols")
        elif safety_score.overall_score >= 70:
            recommendations.append("Good safety profile - proceed with enhanced monitoring")
        elif safety_score.overall_score >= 50:
            recommendations.append("Moderate safety concerns - consider risk mitigation strategies")
        elif safety_score.overall_score >= 30:
            recommendations.append("Significant safety concerns - extensive evaluation required")
        else:
            recommendations.append("High safety risk - consider alternative approaches")
        
        # Specific score-based recommendations
        if safety_score.off_target_score < 50:
            recommendations.append("Off-target risk is primary concern - validate experimentally")
        
        if safety_score.genotoxicity_score < 50:
            recommendations.append("Genotoxicity risk elevated - monitor chromosomal stability")
        
        if safety_score.immunogenicity_score < 50:
            recommendations.append("Immunogenicity risk high - consider immunosuppression")
        
        if safety_score.essential_gene_score < 50:
            recommendations.append("Essential gene impact detected - critical safety evaluation needed")
        
        return recommendations
    
    def _predict_side_effects(self, off_target_analysis: Dict,
                            genotoxicity_analysis: Dict,
                            immunogenicity_analysis: Dict) -> List[SideEffect]:
        """Predict potential side effects based on safety analysis."""
        side_effects = []
        
        # Off-target related side effects
        if off_target_analysis['risk_score'] > 2.0:
            side_effects.append(SideEffect(
                description="Unintended genetic modifications in off-target sites",
                probability=min(0.8, off_target_analysis['risk_score'] / 5.0),
                severity="moderate",
                reversible=False,
                onset_time="immediate",
                evidence_level="theoretical"
            ))
        
        # Essential gene side effects
        if off_target_analysis['essential_analysis']['total_essential_hits'] > 0:
            side_effects.append(SideEffect(
                description="Disruption of essential cellular functions",
                probability=0.3 * off_target_analysis['essential_analysis']['total_essential_hits'],
                severity="severe",
                reversible=False,
                onset_time="days",
                evidence_level="theoretical"
            ))
        
        # Genotoxicity side effects
        if genotoxicity_analysis['overall_score'] > 5.0:
            side_effects.append(SideEffect(
                description="Chromosomal instability and increased mutation rate",
                probability=min(0.6, genotoxicity_analysis['overall_score'] / 10.0),
                severity="severe",
                reversible=False,
                onset_time="weeks",
                evidence_level="theoretical"
            ))
        
        # Immunogenicity side effects
        if immunogenicity_analysis['overall_score'] > 4.0:
            side_effects.append(SideEffect(
                description="Immune response against CRISPR components",
                probability=min(0.7, immunogenicity_analysis['overall_score'] / 8.0),
                severity="moderate",
                reversible=True,
                onset_time="days",
                evidence_level="clinical_trial"
            ))
        
        return side_effects
    
    def _generate_monitoring_requirements(self, safety_score: SafetyScore) -> List[str]:
        """Generate monitoring requirements based on safety score."""
        monitoring = []
        
        # Base monitoring
        monitoring.append("Regular clinical assessment and vital signs monitoring")
        
        # Off-target monitoring
        if safety_score.off_target_score < 70:
            monitoring.append("Periodic whole-genome sequencing to detect off-target effects")
            monitoring.append("Targeted sequencing of predicted off-target sites")
        
        # Genotoxicity monitoring
        if safety_score.genotoxicity_score < 70:
            monitoring.append("Chromosomal stability assessment via karyotyping")
            monitoring.append("Micronucleus assay for DNA damage detection")
            monitoring.append("Long-term cancer surveillance")
        
        # Immunogenicity monitoring
        if safety_score.immunogenicity_score < 70:
            monitoring.append("Anti-Cas protein antibody monitoring")
            monitoring.append("T-cell response assessment")
            monitoring.append("Inflammatory marker surveillance")
        
        # Essential gene monitoring
        if safety_score.essential_gene_score < 70:
            monitoring.append("Functional assessment of essential cellular pathways")
            monitoring.append("Multi-organ system monitoring")
        
        # Risk-based monitoring frequency
        if safety_score.overall_score < 50:
            monitoring.append("Weekly monitoring for first month, then monthly for 6 months")
        elif safety_score.overall_score < 70:
            monitoring.append("Bi-weekly monitoring for first month, then monthly for 3 months")
        else:
            monitoring.append("Standard monitoring schedule with safety checkpoints")
        
        return monitoring
    
    def _generate_risk_mitigation_strategies(self, off_target_analysis: Dict,
                                           genotoxicity_analysis: Dict,
                                           immunogenicity_analysis: Dict) -> List[str]:
        """Generate risk mitigation strategies."""
        strategies = []
        
        # Off-target mitigation
        if off_target_analysis['risk_score'] > 2.0:
            strategies.append("Use high-fidelity Cas variants (e.g., SpCas9-HF1)")
            strategies.append("Optimize guide RNA design for improved specificity")
            strategies.append("Use truncated guide RNAs (17-18 nucleotides)")
            strategies.append("Implement dual-guide approaches for increased specificity")
        
        # Genotoxicity mitigation
        if genotoxicity_analysis['overall_score'] > 5.0:
            strategies.append("Use base editors or prime editors instead of nucleases")
            strategies.append("Minimize Cas protein expression time")
            strategies.append("Co-deliver DNA repair pathway enhancers")
        
        # Immunogenicity mitigation
        if immunogenicity_analysis['overall_score'] > 4.0:
            strategies.append("Use immunosuppressive protocols during treatment")
            strategies.append("Consider humanized Cas proteins")
            strategies.append("Use alternative delivery methods with lower immunogenicity")
            strategies.append("Pre-screen for pre-existing immunity to Cas proteins")
        
        # General mitigation
        strategies.append("Establish clear stopping criteria for adverse events")
        strategies.append("Develop reversal strategies where possible")
        strategies.append("Maintain comprehensive adverse event reporting")
        
        return strategies
    
    def compare_guide_safety(self, guides: List[GuideRNA], target_gene: str) -> List[Tuple[GuideRNA, SafetyScore]]:
        """
        Compare safety of multiple guide RNAs.
        
        Args:
            guides: List of GuideRNA objects to compare
            target_gene: Target gene symbol
            
        Returns:
            List of (GuideRNA, SafetyScore) tuples sorted by safety
        """
        guide_safety_pairs = []
        
        for guide in guides:
            try:
                safety_score = self.analyze_safety(guide, target_gene)
                guide_safety_pairs.append((guide, safety_score))
            except Exception as e:
                self.logger.warning(f"Failed to analyze safety for guide {guide.sequence}: {e}")
        
        # Sort by overall safety score (descending)
        guide_safety_pairs.sort(key=lambda x: x[1].overall_score, reverse=True)
        
        return guide_safety_pairs
    
    def assess_enhancement_safety(self, target_gene: str, enhancement_category: EnhancementCategory,
                                modification_type: str = "base_edit") -> Dict[str, Any]:
        """
        Assess safety considerations specific to enhancement applications.
        
        Args:
            target_gene: Gene to be enhanced
            enhancement_category: Type of enhancement
            modification_type: Type of modification
            
        Returns:
            Enhancement-specific safety assessment
        """
        # Enhancement-specific risk factors
        enhancement_risks = {
            EnhancementCategory.COGNITIVE: {
                'brain_specific_risks': True,
                'developmental_concerns': True,
                'behavioral_effects': True
            },
            EnhancementCategory.PHYSICAL: {
                'muscle_metabolism_effects': True,
                'cardiovascular_impact': True,
                'bone_density_changes': True
            },
            EnhancementCategory.LONGEVITY: {
                'cancer_risk_modulation': True,
                'metabolic_disruption': True,
                'immune_system_effects': True
            }
        }
        
        category_risks = enhancement_risks.get(enhancement_category, {})
        
        # Generate enhancement-specific recommendations
        recommendations = []
        
        if enhancement_category == EnhancementCategory.COGNITIVE:
            recommendations.extend([
                "Monitor neurological function and cognitive assessments",
                "Assess for psychiatric side effects",
                "Consider developmental stage-specific risks"
            ])
        elif enhancement_category == EnhancementCategory.PHYSICAL:
            recommendations.extend([
                "Monitor cardiovascular function during enhanced performance",
                "Assess musculoskeletal system for overuse injuries",
                "Monitor metabolic parameters"
            ])
        elif enhancement_category == EnhancementCategory.LONGEVITY:
            recommendations.extend([
                "Long-term cancer surveillance protocols",
                "Monitor for accelerated aging in other systems",
                "Assess immune system function regularly"
            ])
        
        return {
            'enhancement_category': enhancement_category.value,
            'category_specific_risks': category_risks,
            'recommendations': recommendations,
            'monitoring_duration': self._get_enhancement_monitoring_duration(enhancement_category),
            'reversibility_assessment': self._assess_enhancement_reversibility(
                target_gene, modification_type
            )
        }
    
    def _get_enhancement_monitoring_duration(self, category: EnhancementCategory) -> str:
        """Get recommended monitoring duration for enhancement category."""
        durations = {
            EnhancementCategory.COGNITIVE: "Lifelong (neurological changes)",
            EnhancementCategory.PHYSICAL: "5-10 years (performance monitoring)",
            EnhancementCategory.LONGEVITY: "Lifelong (longevity assessment)",
            EnhancementCategory.SENSORY: "10-20 years (sensory function)",
            EnhancementCategory.METABOLIC: "Lifelong (metabolic stability)",
            EnhancementCategory.IMMUNE: "Lifelong (immune competence)"
        }
        
        return durations.get(category, "Long-term monitoring recommended")
    
    def _assess_enhancement_reversibility(self, target_gene: str, modification_type: str) -> Dict[str, Any]:
        """Assess reversibility of enhancement modifications."""
        reversibility_scores = {
            'base_edit': 0.2,     # Single nucleotide changes - difficult to reverse
            'knockout': 0.1,      # Gene loss - very difficult to reverse
            'knockin': 0.3,       # Integration - potentially reversible
            'activation': 0.8,    # Epigenetic - more reversible
            'interference': 0.9   # RNA-based - highly reversible
        }
        
        base_reversibility = reversibility_scores.get(modification_type, 0.5)
        
        # Gene-specific factors
        gene_factors = 1.0
        
        # Some genes are harder to reverse due to developmental roles
        developmental_genes = {'FOXO3', 'MYC', 'TP53'}
        if target_gene.upper() in developmental_genes:
            gene_factors *= 0.5
        
        overall_reversibility = base_reversibility * gene_factors
        
        return {
            'reversibility_score': overall_reversibility,
            'modification_factor': base_reversibility,
            'gene_factor': gene_factors,
            'reversibility_level': 'High' if overall_reversibility > 0.7 else
                                 'Moderate' if overall_reversibility > 0.4 else 'Low',
            'reversal_strategies': self._get_reversal_strategies(modification_type)
        }
    
    def _get_reversal_strategies(self, modification_type: str) -> List[str]:
        """Get potential reversal strategies for modification type."""
        strategies = {
            'base_edit': [
                "Counter base editing to revert nucleotide change",
                "Gene therapy with wild-type gene copy"
            ],
            'knockout': [
                "Gene replacement therapy",
                "Complementation with functional gene copy"
            ],
            'knockin': [
                "Targeted excision of inserted sequence",
                "Recombinase-mediated removal"
            ],
            'activation': [
                "Remove activating elements",
                "Counter-repressions"
            ],
            'interference': [
                "Discontinue interfering RNA",
                "Natural RNA degradation"
            ]
        }
        
        return strategies.get(modification_type, ["Consult genetic counselor for reversal options"])