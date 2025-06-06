"""
HLA complex analyzer for therapeutic applications.

This module handles the complexity of HLA genes including:
- Haplotype analysis
- Shared epitope recognition
- Linkage disequilibrium
- Population-specific allele frequencies
- Therapeutic target identification
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Set, Union, Any
from dataclasses import dataclass
from collections import defaultdict
import re

from ..models.therapeutic_data_classes import (
    HLAAllele, SharedEpitope, HLAHaplotype, PopulationFrequency,
    TherapeuticTarget, RiskAssessment
)
from ..models.disease_constants import HLA_ALLELE_DATABASE, SHARED_EPITOPES


@dataclass
class HLATypingResult:
    """Results from HLA typing analysis."""
    class_i_alleles: Dict[str, List[str]]  # A, B, C
    class_ii_alleles: Dict[str, List[str]]  # DR, DQ, DP
    confidence_scores: Dict[str, float]
    population_ancestry: Optional[str] = None
    
    @property
    def drb1_alleles(self) -> List[str]:
        """Get DRB1 alleles specifically."""
        return self.class_ii_alleles.get('DRB1', [])


@dataclass
class SharedEpitopeAnalysis:
    """Analysis of shared epitope presence."""
    epitope_sequence: str
    position_range: Tuple[int, int]
    alleles_containing: List[str]
    risk_score: float
    population_frequency: float


class HLAAnalyzer:
    """Comprehensive HLA analysis for therapeutic applications."""
    
    def __init__(self):
        """Initialize HLA analyzer."""
        self.logger = logging.getLogger(__name__)
        self._load_hla_database()
        self._load_shared_epitopes()
        self._load_population_frequencies()
    
    def _load_hla_database(self) -> None:
        """Load comprehensive HLA allele database."""
        # Simplified HLA-DRB1 alleles for demonstration
        self.hla_alleles = {
            'DRB1*01:01': {
                'sequence': self._get_drb1_sequence('01:01'),
                'shared_epitope': 'QKRAA',
                'risk_category': 'high',
                'population_frequencies': {
                    'European': 0.08, 'Asian': 0.01, 'African': 0.02
                }
            },
            'DRB1*04:01': {
                'sequence': self._get_drb1_sequence('04:01'),
                'shared_epitope': 'QKRAA',
                'risk_category': 'very_high',
                'population_frequencies': {
                    'European': 0.12, 'Asian': 0.03, 'African': 0.01
                }
            },
            'DRB1*04:04': {
                'sequence': self._get_drb1_sequence('04:04'),
                'shared_epitope': 'QKRAA',
                'risk_category': 'high',
                'population_frequencies': {
                    'European': 0.03, 'Asian': 0.15, 'African': 0.01
                }
            },
            'DRB1*04:05': {
                'sequence': self._get_drb1_sequence('04:05'),
                'shared_epitope': 'QRRAA',
                'risk_category': 'moderate',
                'population_frequencies': {
                    'European': 0.05, 'Asian': 0.02, 'African': 0.08
                }
            },
            'DRB1*03:01': {
                'sequence': self._get_drb1_sequence('03:01'),
                'shared_epitope': None,
                'risk_category': 'protective',
                'population_frequencies': {
                    'European': 0.15, 'Asian': 0.25, 'African': 0.20
                }
            },
            'DRB1*13:01': {
                'sequence': self._get_drb1_sequence('13:01'),
                'shared_epitope': None,
                'risk_category': 'neutral',
                'population_frequencies': {
                    'European': 0.12, 'Asian': 0.08, 'African': 0.15
                }
            },
            'DRB1*15:01': {
                'sequence': self._get_drb1_sequence('15:01'),
                'shared_epitope': None,
                'risk_category': 'protective',
                'population_frequencies': {
                    'European': 0.18, 'Asian': 0.05, 'African': 0.25
                }
            }
        }
    
    def _get_drb1_sequence(self, allele_code: str) -> str:
        """Get DRB1 protein sequence for allele (simplified)."""
        # This would normally fetch from IMGT/HLA database
        # Simplified example showing key differences at positions 70-74
        base_sequence = (
            "GDTRPRFLEQPKPWEPVPLRPRHHHLSPLSLHPLSLSLAPSCLYLFTLCALTTLVTTTLL"
            "SLLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFLFL"
        )
        
        # Position 70-74 variations (shared epitope region)
        epitope_variants = {
            '01:01': 'QKRAA',  # Shared epitope
            '04:01': 'QKRAA',  # Shared epitope
            '04:04': 'QKRAA',  # Shared epitope
            '04:05': 'QRRAA',  # Modified epitope
            '03:01': 'QRREP',  # No shared epitope
            '13:01': 'QRREG',  # No shared epitope
            '15:01': 'QRREP'   # No shared epitope
        }
        
        epitope = epitope_variants.get(allele_code, 'QRREP')
        # Insert epitope at position 70 (simplified)
        return base_sequence[:70] + epitope + base_sequence[75:]
    
    def _load_shared_epitopes(self) -> None:
        """Load shared epitope definitions."""
        self.shared_epitopes = {
            'SE1': {
                'sequence': 'QKRAA',
                'positions': (70, 74),
                'risk_weight': 1.0,
                'alleles': ['DRB1*01:01', 'DRB1*04:01', 'DRB1*04:04']
            },
            'SE2': {
                'sequence': 'QRRAA',
                'positions': (70, 74),
                'risk_weight': 0.6,
                'alleles': ['DRB1*04:05']
            },
            'SE3': {
                'sequence': 'RRRAA',
                'positions': (70, 74),
                'risk_weight': 0.4,
                'alleles': ['DRB1*04:08']
            }
        }
    
    def _load_population_frequencies(self) -> None:
        """Load population-specific allele frequencies."""
        # This would be loaded from population databases
        self.population_frequencies = {
            'European': {
                'DRB1*15:01': 0.18,
                'DRB1*03:01': 0.15,
                'DRB1*13:01': 0.12,
                'DRB1*04:01': 0.12,
                'DRB1*07:01': 0.10,
                'DRB1*11:01': 0.08,
                'DRB1*01:01': 0.08
            },
            'Asian': {
                'DRB1*03:01': 0.25,
                'DRB1*04:04': 0.15,
                'DRB1*09:01': 0.12,
                'DRB1*13:01': 0.08,
                'DRB1*15:01': 0.05,
                'DRB1*04:01': 0.03
            },
            'African': {
                'DRB1*15:01': 0.25,
                'DRB1*03:01': 0.20,
                'DRB1*13:01': 0.15,
                'DRB1*11:01': 0.12,
                'DRB1*04:05': 0.08
            }
        }
    
    def analyze_hla_typing(self, drb1_alleles: List[str], 
                          population: str = 'European') -> HLATypingResult:
        """
        Analyze HLA typing results.
        
        Args:
            drb1_alleles: List of DRB1 alleles (usually 2)
            population: Population ancestry
            
        Returns:
            HLA typing analysis results
        """
        if len(drb1_alleles) > 2:
            self.logger.warning(f"More than 2 DRB1 alleles provided: {drb1_alleles}")
            drb1_alleles = drb1_alleles[:2]
        
        # Validate alleles
        validated_alleles = []
        confidence_scores = {}
        
        for allele in drb1_alleles:
            if allele in self.hla_alleles:
                validated_alleles.append(allele)
                confidence_scores[allele] = 1.0
            else:
                # Try to find closest match
                closest = self._find_closest_allele(allele)
                if closest:
                    validated_alleles.append(closest)
                    confidence_scores[closest] = 0.8
                    self.logger.warning(f"Allele {allele} not found, using closest: {closest}")
        
        # Create typing result
        typing_result = HLATypingResult(
            class_i_alleles={},
            class_ii_alleles={'DRB1': validated_alleles},
            confidence_scores=confidence_scores,
            population_ancestry=population
        )
        
        return typing_result
    
    def _find_closest_allele(self, query_allele: str) -> Optional[str]:
        """Find closest matching allele in database."""
        # Simple matching by comparing allele groups
        if '*' in query_allele:
            group = query_allele.split('*')[1].split(':')[0]
            for allele in self.hla_alleles.keys():
                if f"*{group}:" in allele:
                    return allele
        return None
    
    def identify_shared_epitopes(self, drb1_alleles: List[str]) -> List[SharedEpitopeAnalysis]:
        """
        Identify shared epitopes in DRB1 alleles.
        
        Args:
            drb1_alleles: List of DRB1 alleles
            
        Returns:
            List of shared epitope analyses
        """
        epitope_analyses = []
        
        for allele in drb1_alleles:
            if allele not in self.hla_alleles:
                continue
            
            allele_data = self.hla_alleles[allele]
            epitope_seq = allele_data.get('shared_epitope')
            
            if epitope_seq:
                # Find which epitope this matches
                matching_epitope = None
                for ep_name, ep_data in self.shared_epitopes.items():
                    if ep_data['sequence'] == epitope_seq:
                        matching_epitope = ep_data
                        break
                
                if matching_epitope:
                    analysis = SharedEpitopeAnalysis(
                        epitope_sequence=epitope_seq,
                        position_range=matching_epitope['positions'],
                        alleles_containing=[allele],
                        risk_score=matching_epitope['risk_weight'],
                        population_frequency=self._calculate_epitope_frequency(epitope_seq)
                    )
                    epitope_analyses.append(analysis)
        
        return epitope_analyses
    
    def _calculate_epitope_frequency(self, epitope_seq: str) -> float:
        """Calculate population frequency of epitope."""
        # Sum frequencies of all alleles containing this epitope
        total_freq = 0.0
        
        for allele, data in self.hla_alleles.items():
            if data.get('shared_epitope') == epitope_seq:
                # Use European frequency as default
                freq = data['population_frequencies'].get('European', 0.0)
                total_freq += freq
        
        return total_freq
    
    def calculate_ra_risk_score(self, drb1_alleles: List[str]) -> Dict[str, Any]:
        """
        Calculate rheumatoid arthritis risk score based on HLA-DRB1.
        
        Args:
            drb1_alleles: List of DRB1 alleles
            
        Returns:
            Risk assessment with detailed breakdown
        """
        # Get shared epitopes
        epitopes = self.identify_shared_epitopes(drb1_alleles)
        
        # Calculate base risk
        base_risk = 1.0  # Population baseline
        
        # Risk from shared epitopes
        epitope_risk = 1.0
        epitope_details = []
        
        for epitope in epitopes:
            if epitope.epitope_sequence == 'QKRAA':
                epitope_risk *= 4.0  # High-risk epitope
            elif epitope.epitope_sequence == 'QRRAA':
                epitope_risk *= 2.5  # Moderate-risk epitope
            
            epitope_details.append({
                'sequence': epitope.epitope_sequence,
                'risk_multiplier': epitope.risk_score,
                'frequency': epitope.population_frequency
            })
        
        # Homozygosity effect (gene dose)
        if len(drb1_alleles) == 2 and drb1_alleles[0] == drb1_alleles[1]:
            if any(ep.epitope_sequence == 'QKRAA' for ep in epitopes):
                epitope_risk *= 1.5  # Additional risk for homozygosity
        
        # Protective alleles
        protective_alleles = []
        for allele in drb1_alleles:
            if allele in self.hla_alleles:
                risk_category = self.hla_alleles[allele]['risk_category']
                if risk_category == 'protective':
                    protective_alleles.append(allele)
                    epitope_risk *= 0.7  # Protective effect
        
        total_risk = base_risk * epitope_risk
        
        # Convert to risk categories
        if total_risk < 1.5:
            risk_category = 'low'
        elif total_risk < 3.0:
            risk_category = 'moderate'
        elif total_risk < 6.0:
            risk_category = 'high'
        else:
            risk_category = 'very_high'
        
        return {
            'total_risk_score': total_risk,
            'risk_category': risk_category,
            'base_risk': base_risk,
            'epitope_risk_multiplier': epitope_risk,
            'shared_epitopes': epitope_details,
            'protective_alleles': protective_alleles,
            'risk_interpretation': self._interpret_risk_score(total_risk),
            'population_percentile': self._calculate_risk_percentile(total_risk)
        }
    
    def _interpret_risk_score(self, risk_score: float) -> str:
        """Interpret risk score for clinical understanding."""
        if risk_score < 1.5:
            return "Below average risk for rheumatoid arthritis"
        elif risk_score < 3.0:
            return "Moderately increased risk - consider lifestyle factors"
        elif risk_score < 6.0:
            return "Significantly increased risk - enhanced monitoring recommended"
        else:
            return "Very high risk - genetic counseling and early screening advised"
    
    def _calculate_risk_percentile(self, risk_score: float) -> float:
        """Calculate what percentile this risk score represents."""
        # Simplified calculation based on risk distribution
        if risk_score < 1.0:
            return 25.0
        elif risk_score < 2.0:
            return 50.0
        elif risk_score < 4.0:
            return 75.0
        elif risk_score < 8.0:
            return 90.0
        else:
            return 95.0
    
    def identify_therapeutic_targets(self, drb1_alleles: List[str]) -> List[TherapeuticTarget]:
        """
        Identify potential therapeutic targets for HLA-based intervention.
        
        Args:
            drb1_alleles: List of DRB1 alleles
            
        Returns:
            List of therapeutic targets with strategies
        """
        targets = []
        
        # Analyze each allele for therapeutic potential
        for allele in drb1_alleles:
            if allele not in self.hla_alleles:
                continue
            
            allele_data = self.hla_alleles[allele]
            risk_category = allele_data['risk_category']
            
            if risk_category in ['high', 'very_high']:
                # High-risk allele - candidate for replacement
                replacement_targets = self._find_replacement_alleles(allele)
                
                for replacement in replacement_targets:
                    target = TherapeuticTarget(
                        target_allele=allele,
                        replacement_allele=replacement,
                        strategy_type='allele_replacement',
                        target_region='exon_2',  # Beta-1 domain
                        expected_risk_reduction=self._calculate_risk_reduction(allele, replacement),
                        technical_feasibility=self._assess_replacement_feasibility(allele, replacement),
                        safety_considerations=self._get_replacement_safety_considerations(allele, replacement)
                    )
                    targets.append(target)
            
            elif allele_data.get('shared_epitope'):
                # Epitope-specific targeting
                epitope_target = TherapeuticTarget(
                    target_allele=allele,
                    replacement_allele=None,
                    strategy_type='epitope_modification',
                    target_region='positions_70_74',
                    expected_risk_reduction=0.6,  # Moderate reduction
                    technical_feasibility=0.4,   # Challenging
                    safety_considerations=['Maintain MHC function', 'Avoid autoimmunity']
                )
                targets.append(epitope_target)
        
        return targets
    
    def _find_replacement_alleles(self, risk_allele: str) -> List[str]:
        """Find suitable replacement alleles for high-risk alleles."""
        replacements = []
        
        # Find alleles with similar population frequency but lower risk
        risk_data = self.hla_alleles[risk_allele]
        risk_freq = risk_data['population_frequencies'].get('European', 0.0)
        
        for allele, data in self.hla_alleles.items():
            if (data['risk_category'] in ['neutral', 'protective'] and
                not data.get('shared_epitope') and
                data['population_frequencies'].get('European', 0.0) >= 0.05):
                replacements.append(allele)
        
        # Prioritize by similarity and safety
        return sorted(replacements, key=lambda x: self.hla_alleles[x]['population_frequencies']['European'], reverse=True)[:3]
    
    def _calculate_risk_reduction(self, original: str, replacement: str) -> float:
        """Calculate expected risk reduction from allele replacement."""
        orig_data = self.hla_alleles[original]
        repl_data = self.hla_alleles[replacement]
        
        # Simplified calculation
        if orig_data['risk_category'] == 'very_high' and repl_data['risk_category'] == 'protective':
            return 0.8  # 80% risk reduction
        elif orig_data['risk_category'] == 'high' and repl_data['risk_category'] == 'neutral':
            return 0.6  # 60% risk reduction
        else:
            return 0.3  # 30% risk reduction
    
    def _assess_replacement_feasibility(self, original: str, replacement: str) -> float:
        """Assess technical feasibility of allele replacement."""
        # Calculate sequence similarity
        orig_seq = self.hla_alleles[original]['sequence']
        repl_seq = self.hla_alleles[replacement]['sequence']
        
        # Count differences
        differences = sum(1 for a, b in zip(orig_seq, repl_seq) if a != b)
        similarity = 1.0 - (differences / len(orig_seq))
        
        # More similar sequences are easier to replace
        if similarity > 0.95:
            return 0.8  # High feasibility
        elif similarity > 0.90:
            return 0.6  # Moderate feasibility
        else:
            return 0.3  # Low feasibility
    
    def _get_replacement_safety_considerations(self, original: str, replacement: str) -> List[str]:
        """Get safety considerations for allele replacement."""
        considerations = [
            "Maintain antigen presentation function",
            "Preserve immune repertoire diversity",
            "Monitor for new autoimmune risks"
        ]
        
        # Add specific considerations
        orig_freq = self.hla_alleles[original]['population_frequencies']['European']
        repl_freq = self.hla_alleles[replacement]['population_frequencies']['European']
        
        if repl_freq < 0.05:
            considerations.append("Low-frequency replacement allele - limited safety data")
        
        if repl_freq > orig_freq * 2:
            considerations.append("Significant frequency change - monitor population effects")
        
        return considerations
    
    def design_hla_correction_strategy(self, drb1_alleles: List[str],
                                     target_risk_reduction: float = 0.7) -> Dict[str, Any]:
        """
        Design comprehensive HLA correction strategy.
        
        Args:
            drb1_alleles: Current DRB1 alleles
            target_risk_reduction: Desired risk reduction (0.0-1.0)
            
        Returns:
            Comprehensive correction strategy
        """
        # Get current risk assessment
        current_risk = self.calculate_ra_risk_score(drb1_alleles)
        
        # Identify therapeutic targets
        targets = self.identify_therapeutic_targets(drb1_alleles)
        
        # Rank strategies by effectiveness and feasibility
        ranked_strategies = []
        
        for target in targets:
            if target.strategy_type == 'allele_replacement':
                effectiveness = target.expected_risk_reduction
                feasibility = target.technical_feasibility
                combined_score = effectiveness * feasibility
                
                ranked_strategies.append({
                    'strategy': target,
                    'effectiveness': effectiveness,
                    'feasibility': feasibility,
                    'score': combined_score
                })
        
        # Sort by combined score
        ranked_strategies.sort(key=lambda x: x['score'], reverse=True)
        
        # Select best strategy
        if ranked_strategies:
            best_strategy = ranked_strategies[0]['strategy']
            
            return {
                'current_risk': current_risk,
                'recommended_strategy': {
                    'type': best_strategy.strategy_type,
                    'target_allele': best_strategy.target_allele,
                    'replacement_allele': best_strategy.replacement_allele,
                    'expected_reduction': best_strategy.expected_risk_reduction,
                    'feasibility': best_strategy.technical_feasibility
                },
                'alternative_strategies': ranked_strategies[1:3],
                'crispr_requirements': self._get_crispr_requirements(best_strategy),
                'monitoring_plan': self._create_monitoring_plan(best_strategy),
                'success_criteria': self._define_success_criteria(target_risk_reduction)
            }
        else:
            return {
                'current_risk': current_risk,
                'recommendation': 'No suitable therapeutic targets identified',
                'alternative_approaches': [
                    'Pharmacological intervention',
                    'Lifestyle modifications',
                    'Regular monitoring'
                ]
            }
    
    def _get_crispr_requirements(self, target: TherapeuticTarget) -> Dict[str, Any]:
        """Get CRISPR design requirements for HLA modification."""
        if target.strategy_type == 'allele_replacement':
            return {
                'approach': 'Homology-directed repair',
                'template_size': 'Long-range (2-5 kb)',
                'target_efficiency': '>60%',
                'delivery_method': 'AAV or electroporation',
                'target_cells': 'Antigen-presenting cells',
                'special_requirements': [
                    'High-fidelity Cas9',
                    'Comprehensive off-target screening',
                    'Functional validation of MHC presentation'
                ]
            }
        elif target.strategy_type == 'epitope_modification':
            return {
                'approach': 'Base editing or prime editing',
                'template_size': 'Short (20-100 bp)',
                'target_efficiency': '>80%',
                'delivery_method': 'Lipid nanoparticles',
                'target_cells': 'Dendritic cells, B cells',
                'special_requirements': [
                    'Precise base editing',
                    'Minimal off-targets in HLA region',
                    'Validation of peptide binding'
                ]
            }
    
    def _create_monitoring_plan(self, target: TherapeuticTarget) -> List[str]:
        """Create monitoring plan for HLA therapeutic intervention."""
        monitoring = [
            "HLA typing confirmation at 1, 3, 6, 12 months",
            "Antigen presentation functional assays",
            "T-cell repertoire analysis",
            "Autoantibody screening",
            "Rheumatoid arthritis disease activity monitoring"
        ]
        
        if target.strategy_type == 'allele_replacement':
            monitoring.extend([
                "Integration site analysis",
                "Off-target sequencing",
                "Immunogenicity assessment"
            ])
        
        return monitoring
    
    def _define_success_criteria(self, target_reduction: float) -> Dict[str, Any]:
        """Define success criteria for HLA intervention."""
        return {
            'primary_endpoints': [
                f"Risk reduction ≥{target_reduction*100:.0f}%",
                "Successful allele modification confirmed by HLA typing",
                "Maintained MHC Class II function"
            ],
            'secondary_endpoints': [
                "No new autoimmune manifestations",
                "Stable T-cell repertoire",
                "No significant off-target effects"
            ],
            'safety_endpoints': [
                "No severe adverse events",
                "No loss of immune function",
                "No opportunistic infections"
            ],
            'timeframe': '12 months post-intervention'
        }