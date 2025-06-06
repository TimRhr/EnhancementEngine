"""
CRISPR guide RNA design and optimization module for Enhancement Engine.

This module provides comprehensive CRISPR-Cas system design including:
- Guide RNA design and scoring
- PAM site identification
- Off-target prediction
- Base editing design
- Prime editing design
"""

import re
import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any, Set
from collections import defaultdict
import itertools

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    try:  # Older BioPython
        from Bio.SeqUtils import GC
    except ImportError:  # Newer BioPython versions
        from Bio.SeqUtils import gc_fraction

        def GC(seq) -> float:
            """Return GC percentage for *seq* as a float."""
            return gc_fraction(seq) * 100

except ImportError:
    raise ImportError("BioPython is required. Install with: pip install biopython")

from ..models.data_classes import (
    GuideRNA,
    PAMSite,
    OffTarget,
    EfficiencyScore,
    CasType,
    VariantPosition,
    RiskLevel,
)
from ..models.constants import (
    PAM_PATTERNS,
    CAS_TYPES,
    DEFAULT_PARAMETERS,
    SCORING_THRESHOLDS,
)


class CRISPRDesignError(Exception):
    """Custom exception for CRISPR design errors."""

    pass


class PAMFinder:
    """Finds PAM sites for different Cas systems."""

    def __init__(self, cas_type: CasType = CasType.CAS9):
        """
        Initialize PAM finder.

        Args:
            cas_type: Type of Cas system
        """
        self.cas_type = cas_type
        self.logger = logging.getLogger(__name__)

        # Get PAM pattern for this Cas type
        if cas_type.value in PAM_PATTERNS:
            self.pam_info = PAM_PATTERNS[cas_type.value]
        else:
            raise CRISPRDesignError(f"Unsupported Cas type: {cas_type}")

    def find_pam_sites(
        self, sequence: Union[str, Seq], strand: str = "both"
    ) -> List[PAMSite]:
        """
        Find all PAM sites in a sequence.

        Args:
            sequence: DNA sequence to search
            strand: Search strand ("+" or "-" or "both")

        Returns:
            List of PAMSite objects
        """
        if isinstance(sequence, Seq):
            sequence = str(sequence)

        sequence = sequence.upper()
        pam_sites = []

        if strand in ["+", "both"]:
            pam_sites.extend(self._find_pam_on_strand(sequence, "+"))

        if strand in ["-", "both"]:
            # Search reverse complement
            rev_comp = str(Seq(sequence).reverse_complement())
            minus_sites = self._find_pam_on_strand(rev_comp, "-")
            # Adjust positions for original sequence
            for site in minus_sites:
                site.position = len(sequence) - site.position - len(site.sequence)
            pam_sites.extend(minus_sites)

        return pam_sites

    def _find_pam_on_strand(self, sequence: str, strand: str) -> List[PAMSite]:
        """Find PAM sites on a specific strand."""
        pam_sites = []
        pam_pattern = self.pam_info["pattern"]

        if not pam_pattern:  # e.g., Cas13 doesn't use PAM
            return pam_sites

        # Convert pattern to regex
        regex_pattern = self._pam_to_regex(pam_pattern)

        for match in re.finditer(regex_pattern, sequence):
            pam_site = PAMSite(
                sequence=match.group(),
                position=match.start(),
                strand=strand,
                cas_type=self.cas_type,
            )

            if pam_site.is_valid:
                pam_sites.append(pam_site)

        return pam_sites

    def _pam_to_regex(self, pam_pattern: str) -> str:
        """Convert PAM pattern to regex."""
        # Standard IUPAC nucleotide codes
        iupac_codes = {
            "N": "[ATCG]",
            "R": "[AG]",  # puRine
            "Y": "[CT]",  # pYrimidine
            "S": "[GC]",  # Strong (3 H bonds)
            "W": "[AT]",  # Weak (2 H bonds)
            "K": "[GT]",  # Keto
            "M": "[AC]",  # aMino
            "B": "[CGT]",  # not A
            "D": "[AGT]",  # not C
            "H": "[ACT]",  # not G
            "V": "[ACG]",  # not T
        }

        regex = ""
        for char in pam_pattern:
            if char in iupac_codes:
                regex += iupac_codes[char]
            else:
                regex += char

        return regex


class GuideEfficiencyScorer:
    """Scores guide RNA efficiency using multiple algorithms."""

    def __init__(self):
        """Initialize efficiency scorer."""
        self.logger = logging.getLogger(__name__)

    def score_guide(
        self, guide_sequence: str, context_sequence: str = ""
    ) -> EfficiencyScore:
        """
        Score guide RNA efficiency.

        Args:
            guide_sequence: 20-nt guide sequence
            context_sequence: Surrounding sequence context (optional)

        Returns:
            EfficiencyScore object
        """
        # Calculate individual scores
        on_target = self._calculate_on_target_score(guide_sequence)
        doench = self._calculate_doench_score(guide_sequence, context_sequence)
        azimuth = self._calculate_azimuth_score(guide_sequence)

        # Combined score (weighted average)
        scores = [s for s in [on_target, doench, azimuth] if s is not None]
        combined = np.mean(scores) if scores else on_target

        return EfficiencyScore(
            on_target_score=on_target,
            doench_score=doench,
            azimuth_score=azimuth,
            combined_score=combined,
        )

    def _calculate_on_target_score(self, guide_sequence: str) -> float:
        """Calculate basic on-target score based on sequence features."""
        if len(guide_sequence) != 20:
            return 0.0

        score = 1.0
        guide = guide_sequence.upper()

        # GC content penalty (optimal range: 40-60%)
        try:
            gc_content = GC(Seq(guide))
            if gc_content < 30 or gc_content > 80:
                score *= 0.5
            elif gc_content < 40 or gc_content > 60:
                score *= 0.8
        except Exception:
            # Fallback GC calculation
            gc_count = guide.count("G") + guide.count("C")
            gc_content = (gc_count / len(guide)) * 100
            if gc_content < 30 or gc_content > 80:
                score *= 0.5
            elif gc_content < 40 or gc_content > 60:
                score *= 0.8

        # Poly-nucleotide stretches penalty
        # FIX: Verhindere max() auf leerer Liste
        for base in "ATCG":
            matches = list(re.finditer(f"{base}+", guide))
            if matches:  # KRITISCHE ÄNDERUNG: Prüfe ob matches existieren
                max_stretch = max(len(m.group()) for m in matches)
                if max_stretch >= 4:
                    score *= 0.7 ** (max_stretch - 3)
            # Wenn keine matches: max_stretch = 0, keine Penalty

        # Position-specific nucleotide preferences
        position_weights = {
            # Based on simplified Doench 2016 model
            0: {"G": 0.1, "A": 0.05, "T": 0.05, "C": 0.0},
            19: {"G": 0.15, "C": 0.1, "A": 0.05, "T": 0.0},
        }

        for pos, weights in position_weights.items():
            if pos < len(guide):
                nucleotide = guide[pos]
                score += weights.get(nucleotide, 0)

        return max(0.0, min(1.0, score))

    def _calculate_doench_score(
        self, guide_sequence: str, context: str
    ) -> Optional[float]:
        """Simplified Doench 2016 scoring (requires 30-nt context)."""
        if len(context) < 30:
            return None

        # This is a greatly simplified version
        # Real implementation would use the full Doench model

        # Find guide position in context
        guide_pos = context.find(guide_sequence)
        if guide_pos == -1:
            return None

        # Extract 30-mer centered on cut site
        cut_site = guide_pos + 17  # 3 bp upstream of PAM
        start = max(0, cut_site - 15)
        end = min(len(context), cut_site + 15)

        if end - start < 20:
            return None

        # Simplified scoring based on nucleotide preferences
        score = 0.5  # Base score

        # Position-specific weights (simplified)
        for i, base in enumerate(context[start:end]):
            pos = i - 15  # Relative to cut site
            if pos == -3:  # Important position
                if base == "C":
                    score += 0.2
            elif pos == -2:
                if base == "G":
                    score += 0.1

        return max(0.0, min(1.0, score))

    def _calculate_azimuth_score(self, guide_sequence: str) -> float:
        """Simplified Azimuth scoring."""
        # Very simplified version of Azimuth algorithm
        score = 0.5
        guide = guide_sequence.upper()

        # Nucleotide composition features
        gc_content = GC(Seq(guide)) / 100
        score += (0.5 - abs(gc_content - 0.5)) * 0.4

        # Simple position-specific scoring
        if len(guide) >= 20:
            if guide[0] == "G":
                score += 0.1
            if guide[19] == "G":
                score += 0.1
            if guide[5:15].count("T") >= 3:
                score -= 0.1

        return max(0.0, min(1.0, score))


class OffTargetPredictor:
    """Predicts off-target binding sites for guide RNAs."""

    def __init__(self, genome_sequence: Optional[str] = None):
        """
        Initialize off-target predictor.

        Args:
            genome_sequence: Reference genome sequence (optional)
        """
        self.genome_sequence = genome_sequence
        self.logger = logging.getLogger(__name__)

    def predict_off_targets(
        self, guide_sequence: str, max_mismatches: int = 4, pam_sequence: str = "NGG"
    ) -> List[OffTarget]:
        """
        Predict off-target sites for a guide RNA.

        Args:
            guide_sequence: 20-nt guide sequence
            max_mismatches: Maximum allowed mismatches
            pam_sequence: PAM sequence pattern

        Returns:
            List of OffTarget objects
        """
        if not self.genome_sequence:
            # Generate synthetic off-targets for demonstration
            return self._generate_synthetic_off_targets(guide_sequence, max_mismatches)

        off_targets = []
        guide = guide_sequence.upper()

        # Search genome for potential off-target sites
        for i in range(len(self.genome_sequence) - len(guide) - 3):
            # Check for PAM site
            pam_site = self.genome_sequence[i + len(guide) : i + len(guide) + 3]
            if not self._matches_pam(pam_site, pam_sequence):
                continue

            # Check guide similarity
            target_seq = self.genome_sequence[i : i + len(guide)]
            mismatches, bulges = self._count_mismatches(guide, target_seq)

            if mismatches <= max_mismatches:
                score = self._calculate_off_target_score(
                    guide, target_seq, mismatches, bulges
                )

                off_target = OffTarget(
                    sequence=target_seq,
                    chromosome="chr1",  # Simplified
                    position=i,
                    strand="+",
                    mismatches=mismatches,
                    bulges=bulges,
                    score=score,
                )

                off_targets.append(off_target)

        return sorted(off_targets, key=lambda x: x.score, reverse=True)

    def _generate_synthetic_off_targets(
        self, guide_sequence: str, max_mismatches: int
    ) -> List[OffTarget]:
        """Generate synthetic off-targets for testing."""
        off_targets = []
        guide = guide_sequence.upper()

        # Generate variants with different numbers of mismatches
        for num_mismatches in range(0, max_mismatches + 1):
            for _ in range(
                max(1, 5 - num_mismatches)
            ):  # More targets with fewer mismatches
                off_target_seq = self._generate_variant_sequence(guide, num_mismatches)

                score = self._calculate_off_target_score(
                    guide, off_target_seq, num_mismatches, 0
                )

                off_target = OffTarget(
                    sequence=off_target_seq,
                    chromosome=f"chr{np.random.randint(1, 23)}",
                    position=np.random.randint(1000000, 100000000),
                    strand=np.random.choice(["+", "-"]),
                    mismatches=num_mismatches,
                    bulges=0,
                    score=score,
                    essential_gene=np.random.random() < 0.1,  # 10% in essential genes
                )

                off_targets.append(off_target)

        return sorted(off_targets, key=lambda x: x.score, reverse=True)

    def _generate_variant_sequence(self, original: str, num_mismatches: int) -> str:
        """Generate a sequence variant with specified number of mismatches."""
        variant = list(original)
        bases = ["A", "T", "C", "G"]

        # Randomly choose positions to mutate
        positions = np.random.choice(len(original), size=num_mismatches, replace=False)

        for pos in positions:
            # Choose a different base
            current_base = variant[pos]
            new_base = np.random.choice([b for b in bases if b != current_base])
            variant[pos] = new_base

        return "".join(variant)

    def _matches_pam(self, sequence: str, pam_pattern: str) -> bool:
        """Check if sequence matches PAM pattern."""
        if len(sequence) != len(pam_pattern):
            return False

        for seq_char, pat_char in zip(sequence.upper(), pam_pattern.upper()):
            if pat_char == "N":
                continue
            elif pat_char != seq_char:
                return False

        return True

    def _count_mismatches(self, seq1: str, seq2: str) -> Tuple[int, int]:
        """Count mismatches and bulges between two sequences."""
        if len(seq1) != len(seq2):
            return 20, abs(len(seq1) - len(seq2))  # High penalty for length difference

        mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
        bulges = 0  # Simplified - no bulge detection

        return mismatches, bulges

    def _calculate_off_target_score(
        self, guide: str, target: str, mismatches: int, bulges: int
    ) -> float:
        """Calculate off-target binding score."""
        # Simplified CFD-like scoring
        if mismatches == 0 and bulges == 0:
            return 1.0

        # Position-weighted mismatch penalties
        score = 1.0

        for i, (g, t) in enumerate(zip(guide, target)):
            if g != t:
                # Higher penalty for mismatches in seed region (positions 1-12)
                if i < 12:
                    score *= 0.1
                else:
                    score *= 0.5

        # Bulge penalty
        score *= 0.5**bulges

        return max(0.0, score)


class BaseEditor:
    """Designs base editing strategies."""

    def __init__(self, editor_type: str = "BE3"):
        """
        Initialize base editor.

        Args:
            editor_type: Type of base editor (BE3, ABE, etc.)
        """
        self.editor_type = editor_type
        self.logger = logging.getLogger(__name__)

        # Base editing windows and preferences
        self.editing_windows = {
            "BE3": (4, 8),  # C->T editing window
            "ABE": (4, 8),  # A->G editing window
            "BE4max": (4, 8),
            "AID": (1, 17),  # Wider window
        }

    def design_base_edit(
        self, target_sequence: str, target_position: int, desired_edit: str
    ) -> Optional[GuideRNA]:
        """
        Design base editing guide RNA.

        Args:
            target_sequence: Target sequence context
            target_position: Position to edit (0-based)
            desired_edit: Desired base change (e.g., "C>T")

        Returns:
            GuideRNA object or None if not possible
        """
        if ">" not in desired_edit:
            self.logger.error("Invalid edit format. Use 'X>Y' format.")
            return None

        from_base, to_base = desired_edit.split(">")
        from_base, to_base = from_base.upper(), to_base.upper()

        # Check if edit is compatible with editor type
        if not self._is_compatible_edit(from_base, to_base):
            self.logger.warning(
                f"Edit {desired_edit} not compatible with {self.editor_type}"
            )
            return None

        # Find suitable guide RNAs
        editing_window = self.editing_windows.get(self.editor_type, (4, 8))

        # Search for guides that place target in editing window
        for guide_start in range(
            max(0, target_position - editing_window[1]),
            min(len(target_sequence) - 20, target_position - editing_window[0] + 1),
        ):

            guide_seq = target_sequence[guide_start : guide_start + 20]

            # Check PAM site
            pam_start = guide_start + 20
            if pam_start + 3 <= len(target_sequence):
                pam_seq = target_sequence[pam_start : pam_start + 3]

                if self._matches_pam(pam_seq, "NGG"):
                    # Calculate position in editing window
                    edit_pos_in_guide = target_position - guide_start

                    if editing_window[0] <= edit_pos_in_guide <= editing_window[1]:
                        # Score the guide
                        scorer = GuideEfficiencyScorer()
                        efficiency = scorer.score_guide(guide_seq, target_sequence)

                        # Create PAM site object
                        pam_site = PAMSite(
                            sequence=pam_seq,
                            position=pam_start,
                            strand="+",
                            cas_type=CasType.BASE_EDITOR,
                        )

                        # Predict off-targets
                        off_target_predictor = OffTargetPredictor()
                        off_targets = off_target_predictor.predict_off_targets(
                            guide_seq
                        )

                        guide_rna = GuideRNA(
                            sequence=guide_seq,
                            pam_site=pam_site,
                            target_position=guide_start,
                            efficiency_score=efficiency,
                            gc_content=GC(Seq(guide_seq)),
                            off_targets=off_targets,
                            cas_type=CasType.BASE_EDITOR,
                        )

                        return guide_rna

        return None

    def _is_compatible_edit(self, from_base: str, to_base: str) -> bool:
        """Check if base edit is compatible with editor type."""
        compatible_edits = {
            "BE3": [("C", "T")],
            "ABE": [("A", "G")],
            "BE4max": [("C", "T")],
            "AID": [("C", "T")],
        }

        allowed = compatible_edits.get(self.editor_type, [])
        return (from_base, to_base) in allowed

    def _matches_pam(self, sequence: str, pam_pattern: str) -> bool:
        """Check if sequence matches PAM pattern."""
        if len(sequence) != len(pam_pattern):
            return False

        for seq_char, pat_char in zip(sequence.upper(), pam_pattern.upper()):
            if pat_char == "N":
                continue
            elif pat_char != seq_char:
                return False

        return True


class CRISPRDesigner:
    """Main CRISPR design class integrating all components."""

    def __init__(self, cas_type: CasType = CasType.CAS9):
        """
        Initialize CRISPR designer.

        Args:
            cas_type: Type of Cas system to use
        """
        self.cas_type = cas_type
        self.pam_finder = PAMFinder(cas_type)
        self.efficiency_scorer = GuideEfficiencyScorer()
        self.off_target_predictor = OffTargetPredictor()
        self.base_editor = BaseEditor()

        self.logger = logging.getLogger(__name__)

        # Get parameters for this Cas type
        self.cas_info = CAS_TYPES.get(cas_type.value, CAS_TYPES["cas9"])
        self.guide_length = self.cas_info["guide_length"]

    def design_guides(
        self,
        sequence: Union[str, SeqRecord],
        target_region: Optional[Tuple[int, int]] = None,
        max_guides: int = 10,
    ) -> List[GuideRNA]:
        """
        Design guide RNAs for a target sequence.

        Args:
            sequence: Target DNA sequence
            target_region: Optional target region (start, end)
            max_guides: Maximum number of guides to return

        Returns:
            List of GuideRNA objects sorted by score
        """
        if isinstance(sequence, SeqRecord):
            sequence = str(sequence.seq)

        sequence = sequence.upper()

        if target_region:
            start, end = target_region
            target_seq = sequence[start:end]
            offset = start
        else:
            target_seq = sequence
            offset = 0

        # Find PAM sites
        pam_sites = self.pam_finder.find_pam_sites(target_seq)

        guides = []

        for pam_site in pam_sites:
            # Extract guide sequence based on Cas type
            guide_seq = self._extract_guide_sequence(target_seq, pam_site)

            if not guide_seq or len(guide_seq) != self.guide_length:
                continue

            # Score efficiency
            context_start = max(0, pam_site.position - 30)
            context_end = min(len(target_seq), pam_site.position + 30)
            context = target_seq[context_start:context_end]

            efficiency = self.efficiency_scorer.score_guide(guide_seq, context)

            # Skip low-efficiency guides
            if (
                efficiency.overall_efficiency
                < DEFAULT_PARAMETERS["guide_design"]["min_efficiency_score"]
            ):
                continue

            # Predict off-targets
            off_targets = self.off_target_predictor.predict_off_targets(
                guide_seq,
                max_mismatches=DEFAULT_PARAMETERS["guide_design"]["max_mismatches"],
            )

            # Filter by safety criteria
            high_risk_off_targets = [ot for ot in off_targets if ot.mismatches <= 2]
            if (
                len(high_risk_off_targets)
                > DEFAULT_PARAMETERS["guide_design"]["max_off_targets"]
            ):
                continue

            # Calculate specificity score
            specificity = self._calculate_specificity_score(off_targets)

            # Adjust PAM site position for original sequence
            pam_site.position += offset

            guide_rna = GuideRNA(
                sequence=guide_seq,
                pam_site=pam_site,
                target_position=pam_site.position - self.guide_length + offset,
                efficiency_score=efficiency,
                gc_content=GC(Seq(guide_seq)),
                off_targets=off_targets,
                specificity_score=specificity,
                cas_type=self.cas_type,
            )

            guides.append(guide_rna)

        # Sort by combined score (efficiency + specificity)
        guides.sort(
            key=lambda g: g.efficiency_score.overall_efficiency * g.specificity_score,
            reverse=True,
        )

        return guides[:max_guides]

    def _extract_guide_sequence(
        self, sequence: str, pam_site: PAMSite
    ) -> Optional[str]:
        """Extract guide sequence based on PAM position and Cas type."""
        pam_pos = pam_site.position

        pam_position = self.pam_finder.pam_info.get("position")

        if pam_position == "3prime":
            # Guide is upstream of PAM (e.g., Cas9)
            guide_start = pam_pos - self.guide_length
            guide_end = pam_pos
        elif pam_position == "5prime":
            # Guide is downstream of PAM (e.g., Cas12a)
            guide_start = pam_pos + len(pam_site.sequence)
            guide_end = guide_start + self.guide_length
        else:
            return None

        if guide_start < 0 or guide_end > len(sequence):
            return None

        guide_seq = sequence[guide_start:guide_end]

        # Reverse complement for minus strand
        if pam_site.strand == "-":
            guide_seq = str(Seq(guide_seq).reverse_complement())

        return guide_seq

    def _calculate_specificity_score(self, off_targets: List[OffTarget]) -> float:
        """Calculate specificity score based on off-targets."""
        if not off_targets:
            return 1.0

        # Weight off-targets by severity
        total_penalty = 0.0

        for ot in off_targets:
            if ot.mismatches == 0:
                penalty = 1.0  # Perfect match
            elif ot.mismatches == 1:
                penalty = 0.5
            elif ot.mismatches == 2:
                penalty = 0.2
            else:
                penalty = 0.05

            # Extra penalty for essential genes
            if ot.essential_gene:
                penalty *= 2.0

            total_penalty += penalty

        # Convert to score (0-1)
        specificity = 1.0 / (1.0 + total_penalty)
        return specificity

    def optimize_guide(
        self, target_sequence: str, target_position: int
    ) -> Optional[GuideRNA]:
        """
        Find the optimal guide RNA for a specific target position.

        Args:
            target_sequence: Target sequence
            target_position: Specific position to target

        Returns:
            Best GuideRNA or None
        """
        # Define search window around target position
        window_size = 100
        start = max(0, target_position - window_size)
        end = min(len(target_sequence), target_position + window_size)

        guides = self.design_guides(target_sequence, (start, end))

        if not guides:
            return None

        # Find guide closest to target position
        best_guide = min(guides, key=lambda g: abs(g.target_position - target_position))

        return best_guide

    def design_base_editor(
        self, target_sequence: str, target_position: int, desired_edit: str
    ) -> Optional[GuideRNA]:
        """
        Design base editing guide RNA.

        Args:
            target_sequence: Target sequence
            target_position: Position to edit
            desired_edit: Desired base change (e.g., "C>T")

        Returns:
            GuideRNA for base editing or None
        """
        return self.base_editor.design_base_edit(
            target_sequence, target_position, desired_edit
        )

    def batch_design(self, targets: List[Tuple[str, str]]) -> Dict[str, List[GuideRNA]]:
        """
        Design guides for multiple targets.

        Args:
            targets: List of (target_name, sequence) tuples

        Returns:
            Dictionary mapping target names to guide lists
        """
        results = {}

        for target_name, sequence in targets:
            try:
                guides = self.design_guides(sequence)
                results[target_name] = guides
                self.logger.info(f"Designed {len(guides)} guides for {target_name}")
            except Exception as e:
                self.logger.error(f"Failed to design guides for {target_name}: {e}")
                results[target_name] = []

        return results

    def get_design_summary(self, guides: List[GuideRNA]) -> Dict[str, Any]:
        """
        Generate summary statistics for designed guides.

        Args:
            guides: List of designed guides

        Returns:
            Summary statistics dictionary
        """
        if not guides:
            return {"total_guides": 0}

        efficiencies = [g.efficiency_score.overall_efficiency for g in guides]
        specificities = [g.specificity_score for g in guides]
        gc_contents = [g.gc_content for g in guides]
        off_target_counts = [len(g.off_targets) for g in guides]

        return {
            "total_guides": len(guides),
            "mean_efficiency": np.mean(efficiencies),
            "mean_specificity": np.mean(specificities),
            "mean_gc_content": np.mean(gc_contents),
            "mean_off_targets": np.mean(off_target_counts),
            "high_efficiency_guides": sum(1 for e in efficiencies if e > 0.7),
            "high_specificity_guides": sum(1 for s in specificities if s > 0.8),
            "safe_guides": sum(1 for g in guides if g.overall_safety_score > 70),
        }
