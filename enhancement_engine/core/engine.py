"""
Main Enhancement Engine orchestrator.

This module provides the primary interface for genetic enhancement analysis,
integrating all components:
- Gene database access
- Sequence analysis
- CRISPR guide design
- Safety assessment
- Effect simulation
- Comprehensive reporting
"""

import logging
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union, Any
from pathlib import Path

from ..utils import (
    save_json,
    load_json,
    save_csv,
    save_text,
    is_valid_email,
    current_timestamp,
)

from .database import GeneDatabase
from .sequence import SequenceAnalyzer
from .crispr import CRISPRDesigner
from .safety import SafetyAnalyzer
from .simulation import EffectSimulator

from ..models.data_classes import (
    EnhancementReport, BatchReport, GeneInfo, GuideRNA, SafetyScore,
    VariantEffect, EnhancementCategory, CasType, ProjectConfig
)
from ..models.constants import (
    ENHANCEMENT_GENES, DEFAULT_PARAMETERS, ERROR_MESSAGES
)


class EnhancementEngineError(Exception):
    """Custom exception for Enhancement Engine errors."""
    pass


class EnhancementEngine:
    """
    Main Enhancement Engine class providing comprehensive genetic enhancement analysis.
    
    This class orchestrates all enhancement analysis components to provide:
    - Gene analysis and variant assessment
    - CRISPR guide design and optimization
    - Safety evaluation and risk assessment
    - Enhancement effect simulation
    - Comprehensive reporting and recommendations
    """
    
    def __init__(self, email: str, config: Optional[ProjectConfig] = None):
        """
        Initialize Enhancement Engine.
        
        Args:
            email: Email address for NCBI access (required)
            config: Optional project configuration
        """
        if not is_valid_email(email):
            raise EnhancementEngineError(ERROR_MESSAGES["invalid_email"])
        
        self.email = email
        self.config = config or ProjectConfig(
            project_name="Enhancement Analysis",
            researcher_email=email
        )
        
        # Initialize all analysis components
        self._initialize_components()
        
        # Setup logging
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Enhancement Engine initialized for {email}")
        
        # Analysis cache for performance
        self._analysis_cache = {}
        
    def _initialize_components(self):
        """Initialize all analysis components."""
        try:
            # Core analysis components
            self.gene_database = GeneDatabase(
                email=self.email,
                cache_enabled=self.config.enable_caching,
                cache_dir=self.config.cache_directory
            )
            
            self.sequence_analyzer = SequenceAnalyzer()
            
            self.crispr_designer = CRISPRDesigner(
                cas_type=self.config.cas_type
            )
            
            self.safety_analyzer = SafetyAnalyzer()
            
            self.effect_simulator = EffectSimulator()
            
            self.logger.info("All analysis components initialized successfully")
            
        except Exception as e:
            self.logger.error(f"Failed to initialize components: {e}")
            raise EnhancementEngineError(f"Component initialization failed: {e}")
    
    def analyze_gene(self, gene_name: str, variant: str = "enhancement_variant",
                    target_tissue: str = "general") -> EnhancementReport:
        """
        Comprehensive analysis of a gene for enhancement potential.
        
        Args:
            gene_name: Gene symbol (e.g., "COMT", "BDNF")
            variant: Target variant for enhancement
            target_tissue: Target tissue/organ for enhancement
            
        Returns:
            EnhancementReport with complete analysis
        """
        try:
            self.logger.info(f"Starting analysis of {gene_name} variant {variant}")
            
            # Check cache first
            cache_key = f"{gene_name}_{variant}_{target_tissue}"
            if cache_key in self._analysis_cache:
                self.logger.debug(f"Retrieved {cache_key} from cache")
                return self._analysis_cache[cache_key]
            
            # Step 1: Gene information and sequence retrieval
            gene_info = self._get_gene_information(gene_name)
            if not gene_info:
                raise EnhancementEngineError(f"Gene {gene_name} not found or not supported")
            
            # Step 2: Sequence analysis
            sequence_analysis = self._analyze_gene_sequence(gene_info)
            
            # Step 3: Identify enhancement target position
            target_position = self._identify_enhancement_target(gene_info, variant)
            
            # Step 4: CRISPR guide design
            guide_design_results = self._design_crispr_guides(
                gene_info, target_position, sequence_analysis
            )
            
            if not guide_design_results['best_guide']:
                raise EnhancementEngineError(f"No suitable CRISPR guides found for {gene_name}")
            
            # Step 5: Safety assessment
            safety_assessment = self._assess_safety(
                guide_design_results['best_guide'], gene_name, variant
            )
            
            # Step 6: Effect simulation
            effect_simulation = self._simulate_enhancement_effects(
                gene_info, variant, target_tissue
            )
            
            # Step 7: Generate recommendations
            recommendations = self._generate_recommendations(
                gene_info, guide_design_results, safety_assessment, effect_simulation
            )
            
            # Step 8: Calculate overall confidence and feasibility
            confidence_score = self._calculate_confidence_score(
                guide_design_results, safety_assessment, effect_simulation
            )
            
            # Create comprehensive report
            enhancement_report = EnhancementReport(
                gene_name=gene_name,
                target_variant=variant,
                analysis_date=datetime.now(),
                best_guide=guide_design_results['best_guide'],
                safety_assessment=safety_assessment,
                predicted_effect=effect_simulation,
                alternative_guides=guide_design_results['alternative_guides'],
                recommendations=recommendations,
                confidence_score=confidence_score,
                warnings=self._generate_warnings(safety_assessment, effect_simulation)
            )
            
            # Cache result
            self._analysis_cache[cache_key] = enhancement_report
            
            self.logger.info(f"Analysis completed for {gene_name}: "
                           f"Feasibility {enhancement_report.feasibility_score:.1f}/100")
            
            return enhancement_report
            
        except Exception as e:
            self.logger.error(f"Analysis failed for {gene_name}: {e}")
            raise EnhancementEngineError(f"Gene analysis failed: {e}")
    
    def _get_gene_information(self, gene_name: str) -> Optional[GeneInfo]:
        """Retrieve comprehensive gene information."""
        try:
            # First check if it's a known enhancement gene
            if self.gene_database.is_enhancement_gene(gene_name):
                gene_info = self.gene_database.search_gene(gene_name)
                self.logger.debug(f"Found enhancement gene: {gene_name}")
                return gene_info
            
            # Search in general database
            gene_info = self.gene_database.search_gene(gene_name)
            if gene_info:
                self.logger.debug(f"Found gene in database: {gene_name}")
                return gene_info
            
            self.logger.warning(f"Gene not found: {gene_name}")
            return None
            
        except Exception as e:
            self.logger.error(f"Failed to retrieve gene information for {gene_name}: {e}")
            return None
    
    def _analyze_gene_sequence(self, gene_info: GeneInfo) -> Dict[str, Any]:
        """Analyze gene sequence for structure and features."""
        try:
            # Get gene sequence
            sequence = None
            if gene_info.refseq_id:
                sequence = self.gene_database.get_sequence(gene_info.refseq_id)
            
            if not sequence:
                self.logger.warning(f"No sequence available for {gene_info.symbol}")
                return {'sequence_available': False}
            
            # Perform sequence analysis
            coding_region = self.sequence_analyzer.find_coding_sequence(sequence)
            composition = self.sequence_analyzer.analyze_sequence_composition(sequence.seq)
            regulatory_elements = self.sequence_analyzer.find_regulatory_elements(sequence.seq)
            
            return {
                'sequence_available': True,
                'sequence_length': len(sequence.seq),
                'coding_region': coding_region,
                'composition': composition,
                'regulatory_elements': regulatory_elements,
                'sequence_record': sequence
            }
            
        except Exception as e:
            self.logger.error(f"Sequence analysis failed: {e}")
            return {'sequence_available': False, 'error': str(e)}
    
    def _identify_enhancement_target(self, gene_info: GeneInfo, variant: str) -> Optional[int]:
        """Identify target position for enhancement editing."""
        try:
            # Known enhancement positions for specific genes
            enhancement_positions = {
                'COMT': {
                    'Val158Met': 472,  # Position in coding sequence
                    'enhancement_variant': 472
                },
                'BDNF': {
                    'Val66Met': 196,
                    'enhancement_variant': 196
                },
                'ACTN3': {
                    'R577X': 1729,
                    'enhancement_variant': 1729
                },
                'FOXO3': {
                    'rs2802292': 1000,  # Simplified position
                    'enhancement_variant': 1000
                },
                'MSTN': {
                    'knockout': 500,  # Target exon for knockout
                    'enhancement_variant': 500
                }
            }
            
            gene_positions = enhancement_positions.get(gene_info.symbol.upper(), {})
            position = gene_positions.get(variant)
            
            if position:
                self.logger.debug(f"Found target position {position} for {gene_info.symbol} {variant}")
                return position
            
            # Default to middle of coding sequence if no specific position
            self.logger.warning(f"No specific position for {variant}, using default")
            return 1000  # Default position
            
        except Exception as e:
            self.logger.error(f"Failed to identify target position: {e}")
            return None
    
    def _design_crispr_guides(self, gene_info: GeneInfo, target_position: int,
                            sequence_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """Design and evaluate CRISPR guides for the target."""
        try:
            if not sequence_analysis['sequence_available']:
                raise EnhancementEngineError("No sequence available for guide design")
            
            sequence_record = sequence_analysis['sequence_record']
            
            # Design guides around target position
            target_region = (
                max(0, target_position - 50),
                min(len(sequence_record.seq), target_position + 50)
            )
            
            guides = self.crispr_designer.design_guides(
                sequence_record, target_region, max_guides=10
            )
            
            if not guides:
                return {
                    'best_guide': None,
                    'alternative_guides': [],
                    'design_summary': {'total_guides': 0}
                }
            
            # Filter guides by safety and efficiency thresholds
            safe_guides = [
                guide for guide in guides
                if (guide.efficiency_score.overall_efficiency >= self.config.efficiency_threshold and
                    guide.overall_safety_score >= self.config.safety_threshold)
            ]
            
            best_guide = safe_guides[0] if safe_guides else guides[0]
            alternative_guides = safe_guides[1:5] if len(safe_guides) > 1 else guides[1:5]
            
            # Generate design summary
            design_summary = self.crispr_designer.get_design_summary(guides)
            
            return {
                'best_guide': best_guide,
                'alternative_guides': alternative_guides,
                'all_guides': guides,
                'design_summary': design_summary
            }
            
        except Exception as e:
            self.logger.error(f"CRISPR guide design failed: {e}")
            return {
                'best_guide': None,
                'alternative_guides': [],
                'error': str(e)
            }
    
    def _assess_safety(self, guide_rna: GuideRNA, gene_name: str, variant: str) -> SafetyScore:
        """Comprehensive safety assessment of the enhancement strategy."""
        try:
            # Determine modification type based on variant
            modification_type = self._determine_modification_type(variant)
            
            # Perform safety analysis
            safety_score = self.safety_analyzer.analyze_safety(
                guide_rna, gene_name, modification_type
            )
            
            # Enhancement-specific safety assessment
            if gene_name.upper() in ENHANCEMENT_GENES.get('cognitive', {}):
                enhancement_category = EnhancementCategory.COGNITIVE
            elif gene_name.upper() in ENHANCEMENT_GENES.get('physical', {}):
                enhancement_category = EnhancementCategory.PHYSICAL
            elif gene_name.upper() in ENHANCEMENT_GENES.get('longevity', {}):
                enhancement_category = EnhancementCategory.LONGEVITY
            else:
                enhancement_category = EnhancementCategory.COGNITIVE  # Default
            
            enhancement_safety = self.safety_analyzer.assess_enhancement_safety(
                gene_name, enhancement_category, modification_type
            )
            
            # Adjust safety score based on enhancement-specific risks
            adjustment_factor = self._calculate_enhancement_safety_adjustment(enhancement_safety)
            safety_score.overall_score *= adjustment_factor
            
            return safety_score
            
        except Exception as e:
            self.logger.error(f"Safety assessment failed: {e}")
            # Return minimal safety score in case of error
            return SafetyScore(
                overall_score=30.0,  # Conservative low score
                off_target_score=50.0,
                genotoxicity_score=50.0,
                immunogenicity_score=50.0,
                essential_gene_score=50.0,
                confidence_level=0.3
            )
    
    def _determine_modification_type(self, variant: str) -> str:
        """Determine the type of genetic modification needed."""
        modification_mapping = {
            'knockout': 'knockout',
            'ko': 'knockout',
            'deletion': 'knockout',
            'val158met': 'base_edit',
            'val66met': 'base_edit',
            'r577x': 'base_edit',
            'overexpression': 'knockin',
            'activation': 'activation'
        }
        
        variant_lower = variant.lower()
        for key, mod_type in modification_mapping.items():
            if key in variant_lower:
                return mod_type
        
        return 'base_edit'  # Default
    
    def _calculate_enhancement_safety_adjustment(self, enhancement_safety: Dict[str, Any]) -> float:
        """Calculate safety adjustment factor for enhancement-specific risks."""
        base_adjustment = 1.0
        
        category_risks = enhancement_safety.get('category_specific_risks', {})
        
        # Adjust based on specific risk categories
        if category_risks.get('brain_specific_risks'):
            base_adjustment *= 0.9  # 10% penalty for brain modifications
        
        if category_risks.get('developmental_concerns'):
            base_adjustment *= 0.85  # 15% penalty for developmental risks
        
        if category_risks.get('cancer_risk_modulation'):
            base_adjustment *= 0.8   # 20% penalty for cancer-related risks
        
        return max(0.5, base_adjustment)  # Don't go below 50% of original score
    
    def _simulate_enhancement_effects(self, gene_info: GeneInfo, variant: str,
                                    target_tissue: str) -> VariantEffect:
        """Simulate the effects of the enhancement modification."""
        try:
            # Determine enhancement category
            enhancement_category = self._determine_enhancement_category(gene_info)
            
            # Simulate variant effect
            variant_effect = self.effect_simulator.simulate_variant_effect(
                gene_info.symbol, variant, enhancement_category
            )
            
            return variant_effect
            
        except Exception as e:
            self.logger.error(f"Effect simulation failed: {e}")
            # Return minimal effect in case of error
            from ..models.data_classes import VariantInfo, ProteinEffect, EnhancementGain
            
            return VariantEffect(
                variant=VariantInfo(name=variant),
                protein_effect=ProteinEffect(
                    stability_change=0.0,
                    activity_change=1.0,
                    confidence_score=0.1
                ),
                enhancement_gain=EnhancementGain(
                    category=EnhancementCategory.COGNITIVE,
                    primary_metric="unknown",
                    baseline_value=5.0,
                    enhanced_value=5.5,
                    improvement_factor=1.1
                ),
                side_effects=[]
            )
    
    def _determine_enhancement_category(self, gene_info: GeneInfo) -> EnhancementCategory:
        """Determine enhancement category for a gene."""
        if gene_info.enhancement_category:
            return gene_info.enhancement_category
        
        # Fallback based on gene function
        gene_symbol = gene_info.symbol.upper()
        
        if gene_symbol in ['COMT', 'BDNF', 'CACNA1C']:
            return EnhancementCategory.COGNITIVE
        elif gene_symbol in ['ACTN3', 'MSTN', 'EPO']:
            return EnhancementCategory.PHYSICAL
        elif gene_symbol in ['FOXO3', 'APOE']:
            return EnhancementCategory.LONGEVITY
        else:
            return EnhancementCategory.COGNITIVE  # Default
    
    def _generate_recommendations(self, gene_info: GeneInfo, guide_results: Dict[str, Any],
                                safety_assessment: SafetyScore, effect_simulation: VariantEffect) -> List[str]:
        """Generate actionable recommendations based on analysis results."""
        recommendations = []
        
        # Overall feasibility assessment
        feasibility_score = (
            guide_results['best_guide'].efficiency_score.overall_efficiency * 100 * 0.3 +
            safety_assessment.overall_score * 0.4 +
            effect_simulation.enhancement_gain.improvement_factor * 50 * 0.3
        )
        
        if feasibility_score >= 80:
            recommendations.append("Highly promising enhancement target - proceed with detailed planning")
        elif feasibility_score >= 60:
            recommendations.append("Viable enhancement option - consider benefits vs risks carefully")
        elif feasibility_score >= 40:
            recommendations.append("Moderate potential - additional safety measures recommended")
        else:
            recommendations.append("Low feasibility - consider alternative approaches")
        
        # Guide-specific recommendations
        best_guide = guide_results['best_guide']
        if best_guide.efficiency_score.overall_efficiency < 0.7:
            recommendations.append("Consider optimizing guide RNA design for better efficiency")
        
        if len(best_guide.high_risk_off_targets) > 0:
            recommendations.append(f"Monitor {len(best_guide.high_risk_off_targets)} high-risk off-target sites")
        
        # Safety-specific recommendations
        if safety_assessment.overall_score < 70:
            recommendations.append("Enhanced safety monitoring protocols required")
        
        if safety_assessment.essential_gene_score < 50:
            recommendations.append("Critical: Potential essential gene impacts detected")
        
        # Effect-specific recommendations
        improvement_factor = effect_simulation.enhancement_gain.improvement_factor
        if improvement_factor > 2.0:
            recommendations.append("Substantial enhancement potential - consider phased implementation")
        elif improvement_factor < 1.2:
            recommendations.append("Limited enhancement effect - evaluate cost-benefit ratio")
        
        # Risk mitigation
        if effect_simulation.side_effects:
            high_risk_effects = [se for se in effect_simulation.side_effects if se.risk_score > 5]
            if high_risk_effects:
                recommendations.append("Implement risk mitigation for identified side effects")
        
        return recommendations
    
    def _calculate_confidence_score(self, guide_results: Dict[str, Any],
                                  safety_assessment: SafetyScore, effect_simulation: VariantEffect) -> float:
        """Calculate overall confidence in the analysis."""
        confidence_factors = []
        
        # Guide design confidence
        if guide_results['best_guide']:
            guide_confidence = min(1.0, guide_results['best_guide'].efficiency_score.overall_efficiency)
            confidence_factors.append(guide_confidence)
        
        # Safety assessment confidence
        confidence_factors.append(safety_assessment.confidence_level)
        
        # Effect simulation confidence (based on available data)
        effect_confidence = 0.7  # Default moderate confidence
        if hasattr(effect_simulation.protein_effect, 'confidence_score'):
            effect_confidence = effect_simulation.protein_effect.confidence_score
        confidence_factors.append(effect_confidence)
        
        # Calculate weighted average
        overall_confidence = sum(confidence_factors) / len(confidence_factors) if confidence_factors else 0.5
        
        return max(0.1, min(1.0, overall_confidence))
    
    def _generate_warnings(self, safety_assessment: SafetyScore, effect_simulation: VariantEffect) -> List[str]:
        """Generate important warnings for the user."""
        warnings = []
        
        # Safety warnings
        if safety_assessment.overall_score < 50:
            warnings.append("WARNING: High safety risk detected - extensive evaluation required")
        
        if safety_assessment.essential_gene_score < 30:
            warnings.append("CRITICAL: Potential damage to essential genes")
        
        # Effect warnings
        high_risk_side_effects = [
            se for se in effect_simulation.side_effects 
            if se.severity == "severe" and se.probability > 0.3
        ]
        
        if high_risk_side_effects:
            warnings.append(f"WARNING: {len(high_risk_side_effects)} severe side effects predicted")
        
        # Irreversible changes warning
        irreversible_effects = [
            se for se in effect_simulation.side_effects 
            if not se.reversible and se.probability > 0.2
        ]
        
        if irreversible_effects:
            warnings.append("WARNING: Potentially irreversible modifications detected")
        
        return warnings
    
    def batch_analysis(self, gene_list: List[str], 
                      variants: Optional[List[str]] = None) -> BatchReport:
        """
        Analyze multiple genes for enhancement potential.
        
        Args:
            gene_list: List of gene symbols to analyze
            variants: Optional list of variants (must match gene_list length)
            
        Returns:
            BatchReport with results for all genes
        """
        try:
            self.logger.info(f"Starting batch analysis of {len(gene_list)} genes")
            
            if variants and len(variants) != len(gene_list):
                raise EnhancementEngineError("Variants list must match gene list length")
            
            if not variants:
                variants = ["enhancement_variant"] * len(gene_list)
            
            batch_report = BatchReport(
                analysis_date=datetime.now(),
                total_genes=len(gene_list),
                successful_analyses=0
            )
            
            for i, gene in enumerate(gene_list):
                try:
                    variant = variants[i]
                    self.logger.debug(f"Analyzing {gene} with variant {variant}")
                    
                    enhancement_report = self.analyze_gene(gene, variant)
                    batch_report.results[gene] = enhancement_report
                    batch_report.successful_analyses += 1
                    
                    self.logger.info(f"Successfully analyzed {gene} "
                                   f"(feasibility: {enhancement_report.feasibility_score:.1f})")
                    
                except Exception as e:
                    error_msg = f"Analysis failed: {str(e)}"
                    batch_report.failed_analyses[gene] = error_msg
                    self.logger.error(f"Failed to analyze {gene}: {e}")
            
            self.logger.info(f"Batch analysis completed: {batch_report.successful_analyses}/"
                           f"{batch_report.total_genes} successful")
            
            return batch_report
            
        except Exception as e:
            self.logger.error(f"Batch analysis failed: {e}")
            raise EnhancementEngineError(f"Batch analysis failed: {e}")
    
    def compare_genes(self, gene_list: List[str]) -> Dict[str, Any]:
        """
        Compare multiple genes for enhancement potential.
        
        Args:
            gene_list: List of gene symbols to compare
            
        Returns:
            Comparison analysis results
        """
        try:
            # Perform batch analysis
            batch_report = self.batch_analysis(gene_list)
            
            if not batch_report.results:
                return {'error': 'No successful analyses to compare'}
            
            # Extract metrics for comparison
            comparison_data = {}
            for gene, report in batch_report.results.items():
                comparison_data[gene] = {
                    'feasibility_score': report.feasibility_score,
                    'safety_score': report.safety_assessment.overall_score,
                    'enhancement_factor': report.predicted_effect.enhancement_gain.improvement_factor,
                    'confidence': report.confidence_score,
                    'guide_efficiency': report.best_guide.efficiency_score.overall_efficiency,
                    'off_target_count': len(report.best_guide.off_targets),
                    'side_effect_count': len(report.predicted_effect.side_effects),
                    'enhancement_category': report.predicted_effect.enhancement_gain.category.value
                }
            
            # Rank genes by different criteria
            rankings = {
                'by_feasibility': sorted(comparison_data.keys(), 
                                       key=lambda x: comparison_data[x]['feasibility_score'], 
                                       reverse=True),
                'by_safety': sorted(comparison_data.keys(), 
                                  key=lambda x: comparison_data[x]['safety_score'], 
                                  reverse=True),
                'by_enhancement': sorted(comparison_data.keys(), 
                                       key=lambda x: comparison_data[x]['enhancement_factor'], 
                                       reverse=True),
                'by_confidence': sorted(comparison_data.keys(), 
                                      key=lambda x: comparison_data[x]['confidence'], 
                                      reverse=True)
            }
            
            # Find best overall candidate
            best_gene = rankings['by_feasibility'][0] if rankings['by_feasibility'] else None
            
            return {
                'comparison_data': comparison_data,
                'rankings': rankings,
                'best_overall_candidate': best_gene,
                'summary_statistics': self._calculate_comparison_statistics(comparison_data),
                'recommendations': self._generate_comparison_recommendations(comparison_data, rankings)
            }
            
        except Exception as e:
            self.logger.error(f"Gene comparison failed: {e}")
            raise EnhancementEngineError(f"Gene comparison failed: {e}")
    
    def _calculate_comparison_statistics(self, comparison_data: Dict[str, Dict[str, float]]) -> Dict[str, float]:
        """Calculate summary statistics for gene comparison."""
        if not comparison_data:
            return {}
        
        import numpy as np
        
        feasibility_scores = [data['feasibility_score'] for data in comparison_data.values()]
        safety_scores = [data['safety_score'] for data in comparison_data.values()]
        enhancement_factors = [data['enhancement_factor'] for data in comparison_data.values()]
        
        return {
            'mean_feasibility': np.mean(feasibility_scores),
            'std_feasibility': np.std(feasibility_scores),
            'mean_safety': np.mean(safety_scores),
            'std_safety': np.std(safety_scores),
            'mean_enhancement': np.mean(enhancement_factors),
            'std_enhancement': np.std(enhancement_factors),
            'high_feasibility_count': sum(1 for score in feasibility_scores if score > 70),
            'safe_options_count': sum(1 for score in safety_scores if score > 70),
            'significant_enhancement_count': sum(1 for factor in enhancement_factors if factor > 1.5)
        }
    
    def _generate_comparison_recommendations(self, comparison_data: Dict[str, Dict[str, float]],
                                           rankings: Dict[str, List[str]]) -> List[str]:
        """Generate recommendations based on gene comparison."""
        recommendations = []
        
        if not comparison_data:
            return ["No data available for recommendations"]
        
        # Best overall recommendation
        best_gene = rankings['by_feasibility'][0]
        best_score = comparison_data[best_gene]['feasibility_score']
        
        if best_score > 80:
            recommendations.append(f"Highly recommend {best_gene} - excellent feasibility score ({best_score:.1f})")
        elif best_score > 60:
            recommendations.append(f"Consider {best_gene} as primary candidate - good feasibility ({best_score:.1f})")
        else:
            recommendations.append("All options show limited feasibility - consider alternative approaches")
        
        # Safety recommendations
        safest_gene = rankings['by_safety'][0]
        if safest_gene != best_gene:
            safety_score = comparison_data[safest_gene]['safety_score']
            recommendations.append(f"For maximum safety, consider {safest_gene} (safety score: {safety_score:.1f})")
        
        # Enhancement potential
        highest_enhancement = rankings['by_enhancement'][0]
        enhancement_factor = comparison_data[highest_enhancement]['enhancement_factor']
        if enhancement_factor > 2.0:
            recommendations.append(f"{highest_enhancement} offers highest enhancement potential ({enhancement_factor:.1f}x)")
        
        # Multi-gene strategy
        top_3_genes = rankings['by_feasibility'][:3]
        if len(top_3_genes) >= 2:
            recommendations.append(f"Consider multi-gene approach with top candidates: {', '.join(top_3_genes)}")
        
        return recommendations
    
    def save_project(self, filename: str, include_cache: bool = False):
        """
        Save current project and analysis results.
        
        Args:
            filename: Output filename
            include_cache: Whether to include analysis cache
        """
        try:
            project_data = {
                'config': self.config.to_dict(),
                'timestamp': datetime.now().isoformat(),
                'version': '0.1.0'
            }
            
            if include_cache:
                # Convert cache to serializable format
                serializable_cache = {}
                for key, report in self._analysis_cache.items():
                    serializable_cache[key] = {
                        'gene_name': report.gene_name,
                        'target_variant': report.target_variant,
                        'feasibility_score': report.feasibility_score,
                        'summary': report.summary
                    }
                project_data['analysis_cache'] = serializable_cache
            
            save_json(project_data, filename)
            
            self.logger.info(f"Project saved to {filename}")
            
        except Exception as e:
            self.logger.error(f"Failed to save project: {e}")
            raise EnhancementEngineError(f"Project save failed: {e}")
    
    def load_project(self, filename: str):
        """
        Load project and analysis results.
        
        Args:
            filename: Input filename
        """
        try:
            project_data = load_json(filename)
            
            # Update configuration
            if 'config' in project_data:
                config_dict = project_data['config']
                self.config = ProjectConfig(
                    project_name=config_dict.get('project_name', 'Loaded Project'),
                    researcher_email=config_dict.get('researcher_email', self.email),
                    cas_type=CasType(config_dict.get('cas_type', 'cas9')),
                    safety_threshold=config_dict.get('safety_threshold', 70.0),
                    efficiency_threshold=config_dict.get('efficiency_threshold', 0.5)
                )
            
            self.logger.info(f"Project loaded from {filename}")
            
        except Exception as e:
            self.logger.error(f"Failed to load project: {e}")
            raise EnhancementEngineError(f"Project load failed: {e}")
    
    def get_enhancement_genes(self, category: Optional[str] = None) -> Dict[str, GeneInfo]:
        """
        Get available enhancement genes.
        
        Args:
            category: Optional category filter ('cognitive', 'physical', 'longevity')
            
        Returns:
            Dictionary of gene symbol -> GeneInfo
        """
        try:
            if category:
                category_enum = EnhancementCategory(category)
                return self.gene_database.get_enhancement_genes(category_enum)
            else:
                return self.gene_database.get_enhancement_genes()
                
        except Exception as e:
            self.logger.error(f"Failed to get enhancement genes: {e}")
            return {}
    
    def get_database_stats(self) -> Dict[str, Any]:
        """Get database and system statistics."""
        try:
            stats = self.gene_database.get_database_stats()
            
            # Add analysis statistics
            stats['analysis_cache_size'] = len(self._analysis_cache)
            stats['config'] = self.config.to_dict()
            
            return stats
            
        except Exception as e:
            self.logger.error(f"Failed to get database stats: {e}")
            return {'error': str(e)}
    
    def clear_cache(self):
        """Clear analysis cache."""
        self._analysis_cache.clear()
        if hasattr(self.gene_database, 'clear_cache'):
            self.gene_database.clear_cache()
        self.logger.info("Analysis cache cleared")
    
    def validate_gene_input(self, gene_name: str) -> Dict[str, Any]:
        """
        Validate gene input and provide suggestions.
        
        Args:
            gene_name: Gene symbol to validate
            
        Returns:
            Validation results and suggestions
        """
        try:
            validation_result = {
                'valid': False,
                'gene_info': None,
                'suggestions': [],
                'is_enhancement_gene': False,
                'available_variants': []
            }
            
            # Check if gene exists
            gene_info = self.gene_database.search_gene(gene_name)
            if gene_info:
                validation_result['valid'] = True
                validation_result['gene_info'] = gene_info
            
            # Check if it's a known enhancement gene
            if self.gene_database.is_enhancement_gene(gene_name):
                validation_result['is_enhancement_gene'] = True
                
                # Get available variants for enhancement genes
                enhancement_variants = {
                    'COMT': ['Val158Met', 'Met158Met'],
                    'BDNF': ['Val66Met', 'Val66Val'],
                    'ACTN3': ['R577X', 'R577R'],
                    'FOXO3': ['rs2802292_T', 'rs2802292_TT'],
                    'MSTN': ['knockout', 'heterozygous_ko']
                }
                
                validation_result['available_variants'] = enhancement_variants.get(
                    gene_name.upper(), ['enhancement_variant']
                )
            
            # Generate suggestions for invalid inputs
            if not validation_result['valid']:
                # Simple fuzzy matching for suggestions
                enhancement_genes = list(ENHANCEMENT_GENES.get('cognitive', {}).keys())
                enhancement_genes.extend(ENHANCEMENT_GENES.get('physical', {}).keys())
                enhancement_genes.extend(ENHANCEMENT_GENES.get('longevity', {}).keys())
                
                suggestions = []
                gene_upper = gene_name.upper()
                
                for gene in enhancement_genes:
                    if gene_upper in gene or gene in gene_upper:
                        suggestions.append(gene)
                
                validation_result['suggestions'] = suggestions[:5]  # Top 5 suggestions
            
            return validation_result
            
        except Exception as e:
            self.logger.error(f"Gene validation failed: {e}")
            return {
                'valid': False,
                'error': str(e),
                'suggestions': []
            }
    
    def estimate_analysis_time(self, gene_list: List[str]) -> Dict[str, float]:
        """
        Estimate analysis time for given genes.
        
        Args:
            gene_list: List of genes to analyze
            
        Returns:
            Time estimates in seconds
        """
        # Base time estimates per analysis step
        time_estimates = {
            'gene_lookup': 2.0,      # seconds
            'sequence_analysis': 5.0,
            'crispr_design': 15.0,
            'safety_assessment': 10.0,
            'effect_simulation': 8.0,
            'report_generation': 3.0
        }
        
        single_gene_time = sum(time_estimates.values())
        
        # Account for caching (50% time reduction for cached results)
        cached_genes = sum(1 for gene in gene_list if f"{gene}_enhancement_variant_general" in self._analysis_cache)
        uncached_genes = len(gene_list) - cached_genes
        
        total_time = (uncached_genes * single_gene_time + 
                     cached_genes * single_gene_time * 0.5)
        
        return {
            'total_genes': len(gene_list),
            'cached_genes': cached_genes,
            'uncached_genes': uncached_genes,
            'estimated_time_seconds': total_time,
            'estimated_time_minutes': total_time / 60,
            'per_gene_breakdown': time_estimates
        }
    
    def get_analysis_summary(self) -> Dict[str, Any]:
        """Get summary of all analyses performed."""
        try:
            if not self._analysis_cache:
                return {'message': 'No analyses performed yet'}
            
            # Analyze cached results
            feasibility_scores = []
            safety_scores = []
            enhancement_factors = []
            categories = []
            
            for report in self._analysis_cache.values():
                feasibility_scores.append(report.feasibility_score)
                safety_scores.append(report.safety_assessment.overall_score)
                enhancement_factors.append(report.predicted_effect.enhancement_gain.improvement_factor)
                categories.append(report.predicted_effect.enhancement_gain.category.value)
            
            import numpy as np
            from collections import Counter
            
            summary = {
                'total_analyses': len(self._analysis_cache),
                'feasibility_stats': {
                    'mean': np.mean(feasibility_scores),
                    'std': np.std(feasibility_scores),
                    'min': np.min(feasibility_scores),
                    'max': np.max(feasibility_scores)
                },
                'safety_stats': {
                    'mean': np.mean(safety_scores),
                    'std': np.std(safety_scores),
                    'min': np.min(safety_scores),
                    'max': np.max(safety_scores)
                },
                'enhancement_stats': {
                    'mean': np.mean(enhancement_factors),
                    'std': np.std(enhancement_factors),
                    'min': np.min(enhancement_factors),
                    'max': np.max(enhancement_factors)
                },
                'category_distribution': dict(Counter(categories)),
                'high_feasibility_count': sum(1 for score in feasibility_scores if score > 70),
                'safe_options_count': sum(1 for score in safety_scores if score > 70),
                'significant_enhancement_count': sum(1 for factor in enhancement_factors if factor > 1.5)
            }
            
            return summary
            
        except Exception as e:
            self.logger.error(f"Failed to generate analysis summary: {e}")
            return {'error': str(e)}
    
    def export_results(self, output_format: str = 'json', filename: Optional[str] = None) -> str:
        """
        Export analysis results in specified format.
        
        Args:
            output_format: 'json', 'csv', or 'html'
            filename: Optional output filename
            
        Returns:
            Output filename or data string
        """
        try:
            if not self._analysis_cache:
                raise EnhancementEngineError("No analysis results to export")
            
            timestamp = current_timestamp()
            
            if output_format.lower() == 'json':
                # Export as JSON
                export_data = {
                    'metadata': {
                        'export_date': datetime.now().isoformat(),
                        'total_analyses': len(self._analysis_cache),
                        'config': self.config.to_dict()
                    },
                    'results': {}
                }
                
                for key, report in self._analysis_cache.items():
                    export_data['results'][key] = {
                        'gene_name': report.gene_name,
                        'target_variant': report.target_variant,
                        'feasibility_score': report.feasibility_score,
                        'safety_score': report.safety_assessment.overall_score,
                        'enhancement_factor': report.predicted_effect.enhancement_gain.improvement_factor,
                        'confidence': report.confidence_score,
                        'summary': report.summary,
                        'recommendations': report.recommendations
                    }
                
                output_filename = filename or f"enhancement_results_{timestamp}.json"
                save_json(export_data, output_filename)
                
                return output_filename
            
            elif output_format.lower() == 'csv':
                # Export as CSV
                output_filename = filename or f"enhancement_results_{timestamp}.csv"

                rows = []
                for key, report in self._analysis_cache.items():
                    rows.append({
                        'Gene': report.gene_name,
                        'Variant': report.target_variant,
                        'Feasibility_Score': f"{report.feasibility_score:.1f}",
                        'Safety_Score': f"{report.safety_assessment.overall_score:.1f}",
                        'Enhancement_Factor': f"{report.predicted_effect.enhancement_gain.improvement_factor:.2f}",
                        'Confidence': f"{report.confidence_score:.2f}",
                        'Category': report.predicted_effect.enhancement_gain.category.value,
                        'Summary': report.summary,
                    })

                save_csv(rows, output_filename, [
                    'Gene', 'Variant', 'Feasibility_Score', 'Safety_Score',
                    'Enhancement_Factor', 'Confidence', 'Category', 'Summary'
                ])
                
                return output_filename
            
            elif output_format.lower() == 'html':
                # Export as HTML report
                html_content = self._generate_html_report()

                output_filename = filename or f"enhancement_report_{timestamp}.html"
                save_text(html_content, output_filename)
                
                return output_filename
            
            else:
                raise EnhancementEngineError(f"Unsupported export format: {output_format}")
                
        except Exception as e:
            self.logger.error(f"Export failed: {e}")
            raise EnhancementEngineError(f"Export failed: {e}")
    
    def _generate_html_report(self) -> str:
        """Generate HTML report of analysis results."""
        html_template = """
<!DOCTYPE html>
<html>
<head>
    <title>Enhancement Engine Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
        .summary { margin: 20px 0; }
        .gene-result { border: 1px solid #ddd; margin: 10px 0; padding: 15px; border-radius: 5px; }
        .high-feasibility { border-left: 5px solid #4CAF50; }
        .medium-feasibility { border-left: 5px solid #FF9800; }
        .low-feasibility { border-left: 5px solid #f44336; }
        .metric { display: inline-block; margin: 5px 10px; padding: 5px; background-color: #f9f9f9; border-radius: 3px; }
        .recommendations { background-color: #e3f2fd; padding: 10px; border-radius: 5px; margin-top: 10px; }
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Enhancement Engine Analysis Report</h1>
        <p><strong>Generated:</strong> {timestamp}</p>
        <p><strong>Total Analyses:</strong> {total_analyses}</p>
        <p><strong>Project:</strong> {project_name}</p>
    </div>
    
    <div class="summary">
        <h2>📊 Summary Statistics</h2>
        {summary_stats}
    </div>
    
    <div class="results">
        <h2>🎯 Gene Analysis Results</h2>
        {gene_results}
    </div>
</body>
</html>
        """
        
        # Generate summary statistics
        summary = self.get_analysis_summary()
        summary_html = f"""
        <div class="metric">Average Feasibility: {summary.get('feasibility_stats', {}).get('mean', 0):.1f}</div>
        <div class="metric">Average Safety: {summary.get('safety_stats', {}).get('mean', 0):.1f}</div>
        <div class="metric">High Feasibility Options: {summary.get('high_feasibility_count', 0)}</div>
        <div class="metric">Safe Options: {summary.get('safe_options_count', 0)}</div>
        """
        
        # Generate gene results
        gene_results_html = ""
        for key, report in self._analysis_cache.items():
            feasibility = report.feasibility_score
            
            if feasibility >= 70:
                css_class = "gene-result high-feasibility"
            elif feasibility >= 50:
                css_class = "gene-result medium-feasibility"
            else:
                css_class = "gene-result low-feasibility"
            
            gene_results_html += f"""
            <div class="{css_class}">
                <h3>{report.gene_name} - {report.target_variant}</h3>
                <div class="metric">Feasibility: {feasibility:.1f}/100</div>
                <div class="metric">Safety: {report.safety_assessment.overall_score:.1f}/100</div>
                <div class="metric">Enhancement: {report.predicted_effect.enhancement_gain.improvement_factor:.1f}x</div>
                <div class="metric">Confidence: {report.confidence_score:.2f}</div>
                <p><strong>Summary:</strong> {report.summary}</p>
                <div class="recommendations">
                    <strong>Recommendations:</strong>
                    <ul>{"".join(f"<li>{rec}</li>" for rec in report.recommendations[:3])}</ul>
                </div>
            </div>
            """
        
        return html_template.format(
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            total_analyses=len(self._analysis_cache),
            project_name=self.config.project_name,
            summary_stats=summary_html,
            gene_results=gene_results_html
        )


# Convenience functions for quick analysis
def quick_analysis(gene_name: str, email: str, variant: str = "enhancement_variant") -> EnhancementReport:
    """
    Quick analysis function for single gene.
    
    Args:
        gene_name: Gene symbol
        email: Email for NCBI access
        variant: Target variant
        
    Returns:
        EnhancementReport
    """
    engine = EnhancementEngine(email)
    return engine.analyze_gene(gene_name, variant)


def compare_enhancement_genes(gene_list: List[str], email: str) -> Dict[str, Any]:
    """
    Quick comparison of multiple enhancement genes.
    
    Args:
        gene_list: List of gene symbols
        email: Email for NCBI access
        
    Returns:
        Comparison results
    """
    engine = EnhancementEngine(email)
    return engine.compare_genes(gene_list)


def get_enhancement_recommendations(category: str, email: str) -> List[str]:
    """
    Get enhancement recommendations for a category.
    
    Args:
        category: Enhancement category ('cognitive', 'physical', 'longevity')
        email: Email for NCBI access
        
    Returns:
        List of recommended genes
    """
    engine = EnhancementEngine(email)
    genes = engine.get_enhancement_genes(category)
    
    # Simple ranking by known enhancement potential
    rankings = {
        'cognitive': ['COMT', 'BDNF', 'CACNA1C'],
        'physical': ['ACTN3', 'MSTN', 'EPO'],
        'longevity': ['FOXO3', 'APOE']
    }
    
    recommended = rankings.get(category, list(genes.keys())[:5])
    return [gene for gene in recommended if gene in genes]