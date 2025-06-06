"""
Therapeutic CRISPR design for disease correction.

This module provides specialized CRISPR design for therapeutic applications,
including base editing, prime editing, and gene replacement strategies.
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any, Set
from collections import defaultdict

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError("BioPython is required. Install with: pip install biopython")

from .crispr import CRISPRDesigner, GuideEfficiencyScorer, OffTargetPredictor
from ..models.data_classes import GuideRNA, PAMSite, CasType, EfficiencyScore
from ..models.therapeutic_data_classes import (
    TherapeuticTarget, CorrectionStrategy, TherapeuticApproach,
    DeliveryMethod, TherapeuticEfficacy
)
from ..models.disease_constants import (
    BASE_EDITING_CONSTRAINTS, PRIME_EDITING_CONSTRAINTS,
    THERAPEUTIC_SAFETY_THRESHOLDS
)


class TherapeuticCRISPRError(Exception):
    """Custom exception for therapeutic CRISPR design errors."""
    pass


class BaseEditingDesigner:
    """Designs base editing strategies for point mutation correction."""
    
    def __init__(self):
        """Initialize base editing designer."""
        self.logger = logging.getLogger(__name__)
        self.efficiency_scorer = GuideEfficiencyScorer()
        self.off_target_predictor = OffTargetPredictor()
    
    def design_correction(self, target_sequence: str, mutation_position: int,
                         original_base: str, target_base: str,
                         editor_type: str = "BE4max") -> Optional[CorrectionStrategy]:
        """
        Design base editing correction strategy.
        
        Args:
            target_sequence: DNA sequence containing the mutation
            mutation_position: Position of mutation (0-based)
            original_base: Current (mutated) base
            target_base: Desired (corrected) base
            editor_type: Base editor type to use
            
        Returns:
            CorrectionStrategy object or None if not possible
        """
        try:
            # Validate correction possibility
            if not self._is_correctable(original_base, target_base, editor_type):
                self.logger.warning(f"Cannot correct {original_base}>{target_base} with {editor_type}")
                return None
            
            # Get editor constraints
            editor_info = self._get_editor_constraints(editor_type)
            if not editor_info:
                raise TherapeuticCRISPRError(f"Unknown editor type: {editor_type}")
            
            # Find suitable guide RNAs
            guides = self._find_base_editing_guides(
                target_sequence, mutation_position, editor_info
            )
            
            if not guides:
                self.logger.warning("No suitable guides found for base editing")
                return None
            
            # Select best guide
            best_guide = self._select_best_guide(guides, mutation_position, editor_info)
            
            # Create correction strategy
            correction_strategy = CorrectionStrategy(
                original_sequence=original_base,
                corrected_sequence=target_base,
                editing_method="base_editing",
                guide_rnas=[best_guide],
                correction_window=editor_info["editing_window"]
            )
            
            self.logger.info(f"Designed {editor_type} correction for {original_base}>{target_base}")
            return correction_strategy
            
        except Exception as e:
            self.logger.error(f"Base editing design failed: {e}")
            raise TherapeuticCRISPRError(f"Failed to design base editing: {e}")
    
    def _is_correctable(self, original_base: str, target_base: str, editor_type: str) -> bool:
        """Check if base change is possible with given editor."""
        editor_info = self._get_editor_constraints(editor_type)
        if not editor_info:
            return False
        
        # Check if this is a valid base conversion
        valid_conversions = {
            "C>T": ["BE3", "BE4max", "AID_BE3"],
            "T>C": [],  # Not directly possible with current editors
            "A>G": ["ABE7.10", "ABE8e"],
            "G>A": []   # Not directly possible
        }
        
        conversion = f"{original_base}>{target_base}"
        return editor_type in valid_conversions.get(conversion, [])
    
    def _get_editor_constraints(self, editor_type: str) -> Optional[Dict[str, Any]]:
        """Get constraints for specific base editor."""
        # Check cytosine base editors
        if editor_type in BASE_EDITING_CONSTRAINTS["cytosine_base_editors"]:
            return BASE_EDITING_CONSTRAINTS["cytosine_base_editors"][editor_type]
        
        # Check adenine base editors
        if editor_type in BASE_EDITING_CONSTRAINTS["adenine_base_editors"]:
            return BASE_EDITING_CONSTRAINTS["adenine_base_editors"][editor_type]
        
        return None
    
    def _find_base_editing_guides(self, sequence: str, target_pos: int,
                                 editor_info: Dict[str, Any]) -> List[GuideRNA]:
        """Find guide RNAs that place target in editing window."""
        guides = []
        window_start, window_end = editor_info["editing_window"]
        
        # Search for guides that place target position in editing window
        for guide_start in range(max(0, target_pos - window_end - 3),
                                min(len(sequence) - 23, target_pos - window_start + 4)):
            
            guide_seq = sequence[guide_start:guide_start + 20]
            pam_seq = sequence[guide_start + 20:guide_start + 23]
            
            # Check PAM compatibility (assuming NGG for now)
            if not pam_seq.endswith("GG"):
                continue
            
            # Check if target is in editing window
            target_in_guide = target_pos - guide_start
            if window_start <= target_in_guide <= window_end:
                
                # Score the guide
                efficiency = self.efficiency_scorer.score_guide(guide_seq, sequence)
                
                # Predict off-targets
                off_targets = self.off_target_predictor.predict_off_targets(guide_seq)
                
                # Create PAM site
                pam_site = PAMSite(
                    sequence=pam_seq,
                    position=guide_start + 20,
                    strand="+",
                    cas_type=CasType.BASE_EDITOR
                )
                
                # Create guide RNA
                guide = GuideRNA(
                    sequence=guide_seq,
                    pam_site=pam_site,
                    target_position=guide_start,
                    efficiency_score=efficiency,
                    gc_content=self._calculate_gc_content(guide_seq),
                    off_targets=off_targets,
                    cas_type=CasType.BASE_EDITOR
                )
                
                guides.append(guide)
        
        return guides
    
    def _select_best_guide(self, guides: List[GuideRNA], target_pos: int,
                          editor_info: Dict[str, Any]) -> GuideRNA:
        """Select best guide based on efficiency and safety."""
        if not guides:
            raise TherapeuticCRISPRError("No guides available for selection")
        
        scored_guides = []
        
        for guide in guides:
            # Calculate combined score
            efficiency_score = guide.efficiency_score.overall_efficiency
            safety_score = guide.overall_safety_score / 100
            
            # Prefer guides with target in center of editing window
            window_start, window_end = editor_info["editing_window"]
            window_center = (window_start + window_end) / 2
            target_in_guide = target_pos - guide.target_position
            position_penalty = abs(target_in_guide - window_center) / window_center
            
            combined_score = (efficiency_score * 0.4 + 
                            safety_score * 0.4 + 
                            (1 - position_penalty) * 0.2)
            
            scored_guides.append((guide, combined_score))
        
        # Sort by score and return best
        scored_guides.sort(key=lambda x: x[1], reverse=True)
        return scored_guides[0][0]
    
    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content of sequence."""
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100 if sequence else 0


class PrimeEditingDesigner:
    """Designs prime editing strategies for complex corrections."""
    
    def __init__(self):
        """Initialize prime editing designer."""
        self.logger = logging.getLogger(__name__)
        self.efficiency_scorer = GuideEfficiencyScorer()
        self.off_target_predictor = OffTargetPredictor()
    
    def design_correction(self, target_sequence: str, target_position: int,
                         original_sequence: str, corrected_sequence: str,
                         editor_type: str = "PE3") -> Optional[CorrectionStrategy]:
        """
        Design prime editing correction strategy.
        
        Args:
            target_sequence: DNA sequence containing the target
            target_position: Position where edit should occur
            original_sequence: Current sequence to be replaced
            corrected_sequence: Desired sequence after editing
            editor_type: Prime editor type
            
        Returns:
            CorrectionStrategy object or None if not possible
        """
        try:
            # Check if edit is within prime editing constraints
            constraints = PRIME_EDITING_CONSTRAINTS.get(editor_type)
            if not constraints:
                raise TherapeuticCRISPRError(f"Unknown prime editor: {editor_type}")
            
            if not self._validate_edit_size(original_sequence, corrected_sequence, constraints):
                return None
            
            # Design pegRNA (prime editing guide RNA)
            pegrna = self._design_pegrna(target_sequence, target_position, 
                                       original_sequence, corrected_sequence)
            
            if not pegrna:
                return None
            
            # Design nicking guide (for PE3)
            nicking_guide = self._design_nicking_guide(target_sequence, target_position)
            
            guides = [pegrna]
            if nicking_guide:
                guides.append(nicking_guide)
            
            correction_strategy = CorrectionStrategy(
                original_sequence=original_sequence,
                corrected_sequence=corrected_sequence,
                editing_method="prime_editing",
                guide_rnas=guides,
                template_sequence=self._design_rt_template(corrected_sequence)
            )
            
            self.logger.info(f"Designed prime editing correction: {len(original_sequence)}>{len(corrected_sequence)} bp")
            return correction_strategy
            
        except Exception as e:
            self.logger.error(f"Prime editing design failed: {e}")
            raise TherapeuticCRISPRError(f"Failed to design prime editing: {e}")
    
    def _validate_edit_size(self, original: str, corrected: str, 
                          constraints: Dict[str, Any]) -> bool:
        """Validate that edit is within size constraints."""
        size_diff = len(corrected) - len(original)
        
        if size_diff > 0:  # Insertion
            return size_diff <= constraints["max_insertion"]
        elif size_diff < 0:  # Deletion
            return abs(size_diff) <= constraints["max_deletion"]
        else:  # Replacement
            return len(corrected) <= constraints["max_replacement"]
    
    def _design_pegrna(self, sequence: str, position: int, 
                      original: str, corrected: str) -> Optional[GuideRNA]:
        """Design pegRNA for prime editing."""
        # Simplified pegRNA design - would need more sophisticated implementation
        # This is a placeholder for the actual complex pegRNA design algorithm
        
        # Find suitable spacer sequence
        spacer_start = max(0, position - 30)
        spacer_end = min(len(sequence), position + 30)
        
        for i in range(spacer_start, spacer_end - 20):
            spacer = sequence[i:i + 20]
            pam = sequence[i + 20:i + 23]
            
            if pam.endswith("GG"):  # NGG PAM
                # Score efficiency
                efficiency = self.efficiency_scorer.score_guide(spacer, sequence)
                
                # Create guide
                pam_site = PAMSite(
                    sequence=pam,
                    position=i + 20,
                    strand="+",
                    cas_type=CasType.PRIME_EDITOR
                )
                
                guide = GuideRNA(
                    sequence=spacer,
                    pam_site=pam_site,
                    target_position=i,
                    efficiency_score=efficiency,
                    gc_content=self._calculate_gc_content(spacer),
                    off_targets=[],  # Would predict in real implementation
                    cas_type=CasType.PRIME_EDITOR
                )
                
                return guide
        
        return None
    
    def _design_nicking_guide(self, sequence: str, position: int) -> Optional[GuideRNA]:
        """Design nicking guide for PE3."""
        # Simplified nicking guide design
        # Should be 40-90 bp away from pegRNA cut site
        
        nicking_start = position + 40
        nicking_end = min(len(sequence) - 23, position + 90)
        
        for i in range(nicking_start, nicking_end):
            spacer = sequence[i:i + 20]
            pam = sequence[i + 20:i + 23]
            
            if pam.endswith("GG"):
                efficiency = self.efficiency_scorer.score_guide(spacer, sequence)
                
                pam_site = PAMSite(
                    sequence=pam,
                    position=i + 20,
                    strand="+",
                    cas_type=CasType.CAS9
                )
                
                guide = GuideRNA(
                    sequence=spacer,
                    pam_site=pam_site,
                    target_position=i,
                    efficiency_score=efficiency,
                    gc_content=self._calculate_gc_content(spacer),
                    off_targets=[],
                    cas_type=CasType.CAS9
                )
                
                return guide
        
        return None
    
    def _design_rt_template(self, corrected_sequence: str) -> str:
        """Design reverse transcription template."""
        # Simplified RT template design
        # Should include the corrected sequence plus flanking regions
        return corrected_sequence  # Placeholder
    
    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content of sequence."""
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100 if sequence else 0


class GeneReplacementDesigner:
    """Designs gene replacement strategies for complex corrections."""
    
    def __init__(self):
        """Initialize gene replacement designer."""
        self.logger = logging.getLogger(__name__)
        self.crispr_designer = CRISPRDesigner()
    
    def design_allele_replacement(self, target_gene: str, 
                                harmful_allele: str, protective_allele: str,
                                replacement_strategy: str = "HDR") -> Optional[CorrectionStrategy]:
        """
        Design complete allele replacement strategy.
        
        Args:
            target_gene: Gene to modify
            harmful_allele: Current harmful allele
            protective_allele: Desired protective allele
            replacement_strategy: Method for replacement
            
        Returns:
            CorrectionStrategy object
        """
        try:
            # For HLA genes, this would involve complex allele replacement
            if target_gene.startswith("HLA-"):
                return self._design_hla_replacement(harmful_allele, protective_allele)
            
            # For other genes, design standard HDR
            return self._design_hdr_replacement(target_gene, harmful_allele, protective_allele)
            
        except Exception as e:
            self.logger.error(f"Gene replacement design failed: {e}")
            raise TherapeuticCRISPRError(f"Failed to design replacement: {e}")
    
    def _design_hla_replacement(self, harmful_allele: str, protective_allele: str) -> CorrectionStrategy:
        """Design HLA allele replacement strategy."""
        # This is a complex procedure requiring:
        # 1. Targeted deletion of harmful allele
        # 2. Insertion of protective allele
        # 3. Careful preservation of regulatory elements
        
        # Placeholder implementation
        correction_strategy = CorrectionStrategy(
            original_sequence=harmful_allele,
            corrected_sequence=protective_allele,
            editing_method="HDR_allele_replacement",
            guide_rnas=[],  # Would design multiple guides
            template_sequence=protective_allele
        )
        
        self.logger.info(f"Designed HLA replacement: {harmful_allele} -> {protective_allele}")
        return correction_strategy
    
    def _design_hdr_replacement(self, gene: str, original: str, corrected: str) -> CorrectionStrategy:
        """Design homology-directed repair replacement."""
        # Placeholder for HDR template design
        correction_strategy = CorrectionStrategy(
            original_sequence=original,
            corrected_sequence=corrected,
            editing_method="HDR",
            guide_rnas=[],  # Would design flanking guides
            template_sequence=corrected
        )
        
        return correction_strategy


class TherapeuticCRISPRDesigner:
    """Main therapeutic CRISPR designer integrating all approaches."""
    
    def __init__(self):
        """Initialize therapeutic CRISPR designer."""
        self.base_editor = BaseEditingDesigner()
        self.prime_editor = PrimeEditingDesigner()
        self.gene_replacer = GeneReplacementDesigner()
        self.logger = logging.getLogger(__name__)
    
    def design_therapeutic_intervention(self, therapeutic_target: TherapeuticTarget,
                                      target_sequence: str) -> List[CorrectionStrategy]:
        """
        Design comprehensive therapeutic intervention.
        
        Args:
            therapeutic_target: Target definition
            target_sequence: DNA sequence context
            
        Returns:
            List of possible correction strategies
        """
        strategies = []
        
        try:
            approach = therapeutic_target.therapeutic_approach
            
            if approach == TherapeuticApproach.CORRECTION:
                # Try base editing first for point mutations
                base_strategy = self._try_base_editing(therapeutic_target, target_sequence)
                if base_strategy:
                    strategies.append(base_strategy)
                
                # Try prime editing for more complex edits
                prime_strategy = self._try_prime_editing(therapeutic_target, target_sequence)
                if prime_strategy:
                    strategies.append(prime_strategy)
            
            elif approach == TherapeuticApproach.REPLACEMENT:
                replacement_strategy = self._try_gene_replacement(therapeutic_target)
                if replacement_strategy:
                    strategies.append(replacement_strategy)
            
            elif approach == TherapeuticApproach.KNOCKOUT:
                knockout_strategy = self._design_knockout(therapeutic_target, target_sequence)
                if knockout_strategy:
                    strategies.append(knockout_strategy)
            
            elif approach == TherapeuticApproach.SILENCING:
                silencing_strategy = self._design_silencing(therapeutic_target, target_sequence)
                if silencing_strategy:
                    strategies.append(silencing_strategy)
            
            self.logger.info(f"Designed {len(strategies)} therapeutic strategies")
            return strategies
            
        except Exception as e:
            self.logger.error(f"Therapeutic design failed: {e}")
            raise TherapeuticCRISPRError(f"Failed to design intervention: {e}")
    
    def _try_base_editing(self, target: TherapeuticTarget, sequence: str) -> Optional[CorrectionStrategy]:
        """Attempt base editing approach."""
        # Parse variant information to determine base change needed
        variant = target.target_variant
        
        # For PTPN22 R620W: need T>C correction (revert C1858T mutation)
        if "R620W" in variant:
            return self.base_editor.design_correction(
                sequence, target.target_position, "T", "C", "ABE8e"
            )
        
        return None
    
    def _try_prime_editing(self, target: TherapeuticTarget, sequence: str) -> Optional[CorrectionStrategy]:
        """Attempt prime editing approach."""
        if target.correction_sequence:
            original = sequence[target.target_position:target.target_position + len(target.target_variant)]
            return self.prime_editor.design_correction(
                sequence, target.target_position, original, target.correction_sequence
            )
        return None
    
    def _try_gene_replacement(self, target: TherapeuticTarget) -> Optional[CorrectionStrategy]:
        """Attempt gene replacement approach."""
        if target.correction_sequence:
            return self.gene_replacer.design_allele_replacement(
                target.gene_symbol, target.target_variant, target.correction_sequence
            )
        return None
    
    def _design_knockout(self, target: TherapeuticTarget, sequence: str) -> Optional[CorrectionStrategy]:
        """Design knockout strategy."""
        # Use standard CRISPR to knockout harmful gene function
        guides = self.gene_replacer.crispr_designer.design_guides(sequence, max_guides=3)
        
        if guides:
            return CorrectionStrategy(
                original_sequence="functional",
                corrected_sequence="knockout",
                editing_method="knockout",
                guide_rnas=guides[:2]  # Use top 2 guides for double knockout
            )
        return None
    
    def _design_silencing(self, target: TherapeuticTarget, sequence: str) -> Optional[CorrectionStrategy]:
        """Design gene silencing strategy using CRISPRi."""
        # Design guides targeting promoter region for gene silencing
        promoter_region = sequence[:500]  # Simplified promoter region
        guides = self.gene_replacer.crispr_designer.design_guides(promoter_region, max_guides=5)
        
        if guides:
            # Modify guides for dCas9-KRAB (CRISPRi)
            for guide in guides:
                guide.cas_type = CasType.CAS9  # Would be dCas9 in practice
            
            return CorrectionStrategy(
                original_sequence="active_expression",
                corrected_sequence="silenced_expression",
                editing_method="CRISPRi",
                guide_rnas=guides[:3]  # Use top 3 guides for effective silencing
            )
        return None
    
    def evaluate_strategy_feasibility(self, strategy: CorrectionStrategy) -> float:
        """
        Evaluate feasibility of correction strategy.
        
        Args:
            strategy: Correction strategy to evaluate
            
        Returns:
            Feasibility score (0-1)
        """
        base_feasibility = {
            "base_editing": 0.8,
            "prime_editing": 0.6,
            "HDR": 0.4,
            "knockout": 0.9,
            "CRISPRi": 0.7
        }
        
        method_score = base_feasibility.get(strategy.editing_method, 0.5)
        complexity_penalty = strategy.complexity_score * 0.3
        
        return max(0.0, method_score - complexity_penalty)
    
    def predict_therapeutic_efficacy(self, strategies: List[CorrectionStrategy],
                                   delivery_method: DeliveryMethod,
                                   target_tissue: str) -> TherapeuticEfficacy:
        """
        Predict therapeutic efficacy of correction strategies.
        
        Args:
            strategies: List of correction strategies
            delivery_method: Delivery method
            target_tissue: Target tissue
            
        Returns:
            TherapeuticEfficacy object
        """
        if not strategies:
            return TherapeuticEfficacy(
                correction_efficiency=0.0,
                tissue_penetration=0.0,
                persistence_months=0.0,
                clinical_improvement=0.0,
                remission_probability=0.0
            )
        
        # Use best strategy for prediction
        best_strategy = max(strategies, key=self.evaluate_strategy_feasibility)
        
        # Base efficacy by method
        method_efficacies = {
            "base_editing": 0.7,
            "prime_editing": 0.5,
            "HDR": 0.3,
            "knockout": 0.8,
            "CRISPRi": 0.6
        }
        
        base_efficacy = method_efficacies.get(best_strategy.editing_method, 0.5)
        
        # Tissue penetration by delivery method
        penetration_scores = {
            DeliveryMethod.INTRA_ARTICULAR: 0.8,
            DeliveryMethod.INTRAVENOUS: 0.6,
            DeliveryMethod.LOCAL_INJECTION: 0.9,
            DeliveryMethod.INTRATHECAL: 0.7
        }
        
        tissue_penetration = penetration_scores.get(delivery_method, 0.5)
        
        # Estimate clinical improvement based on correction efficiency
        clinical_improvement = base_efficacy * 0.6  # Conservative estimate
        
        return TherapeuticEfficacy(
            correction_efficiency=base_efficacy,
            tissue_penetration=tissue_penetration,
            persistence_months=6.0,  # Default 6 months
            clinical_improvement=clinical_improvement,
            remission_probability=clinical_improvement * 0.8,
            time_to_effect_days=30
        )