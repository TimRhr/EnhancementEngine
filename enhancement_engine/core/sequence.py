"""
Sequence analysis module for Enhancement Engine.

This module provides comprehensive DNA, RNA, and protein sequence analysis
including:
- Coding sequence identification
- Variant position mapping
- Codon usage analysis
- Secondary structure prediction
- Conservation analysis
"""

import re
import logging
from typing import Dict, List, Optional, Tuple, Union, Any, Set
from collections import Counter
import numpy as np

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqUtils import GC, molecular_weight
    from Bio import SeqIO
except ImportError:
    raise ImportError("BioPython is required. Install with: pip install biopython")

from ..models.data_classes import (
    CodingRegion, VariantPosition, VariantInfo, GeneInfo, 
    ProteinEffect, EnhancementCategory
)
from ..models.constants import DEFAULT_PARAMETERS


class SequenceAnalysisError(Exception):
    """Custom exception for sequence analysis errors."""
    pass


class CodonAnalyzer:
    """Analyzes codon usage and optimization."""
    
    # Standard genetic code
    GENETIC_CODE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    # Human codon usage frequencies (simplified)
    HUMAN_CODON_USAGE = {
        'F': {'TTT': 0.45, 'TTC': 0.55},
        'L': {'TTA': 0.07, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.41},
        'S': {'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.05, 'AGT': 0.15, 'AGC': 0.24},
        'Y': {'TAT': 0.43, 'TAC': 0.57},
        'C': {'TGT': 0.45, 'TGC': 0.55},
        'W': {'TGG': 1.00},
        'P': {'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11},
        'H': {'CAT': 0.41, 'CAC': 0.59},
        'Q': {'CAA': 0.25, 'CAG': 0.75},
        'R': {'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, 'AGA': 0.20, 'AGG': 0.20},
        'I': {'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16},
        'M': {'ATG': 1.00},
        'T': {'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12},
        'N': {'AAT': 0.46, 'AAC': 0.54},
        'K': {'AAA': 0.42, 'AAG': 0.58},
        'V': {'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47},
        'A': {'GCT': 0.26, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11},
        'D': {'GAT': 0.46, 'GAC': 0.54},
        'E': {'GAA': 0.42, 'GAG': 0.58},
        'G': {'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25},
        '*': {'TAA': 0.28, 'TAG': 0.20, 'TGA': 0.52}
    }
    
    def __init__(self):
        """Initialize codon analyzer."""
        self.logger = logging.getLogger(__name__)
    
    def translate_sequence(self, dna_seq: Union[str, Seq]) -> str:
        """
        Translate DNA sequence to amino acids.
        
        Args:
            dna_seq: DNA sequence
            
        Returns:
            Amino acid sequence
        """
        if isinstance(dna_seq, str):
            dna_seq = Seq(dna_seq)
        
        return str(dna_seq.translate())
    
    def get_codon_usage(self, sequence: Union[str, Seq]) -> Dict[str, float]:
        """
        Calculate codon usage frequencies.
        
        Args:
            sequence: Coding DNA sequence
            
        Returns:
            Dictionary of codon -> frequency
        """
        if isinstance(sequence, str):
            sequence = Seq(sequence)
        
        # Ensure sequence length is multiple of 3
        seq_len = len(sequence) - (len(sequence) % 3)
        sequence = sequence[:seq_len]
        
        # Count codons
        codon_counts = Counter()
        for i in range(0, len(sequence), 3):
            codon = str(sequence[i:i+3])
            if len(codon) == 3:
                codon_counts[codon] += 1
        
        # Convert to frequencies
        total_codons = sum(codon_counts.values())
        if total_codons == 0:
            return {}
        
        return {codon: count / total_codons for codon, count in codon_counts.items()}
    
    def calculate_codon_adaptation_index(self, sequence: Union[str, Seq]) -> float:
        """
        Calculate Codon Adaptation Index (CAI) for human codon usage.
        
        Args:
            sequence: Coding DNA sequence
            
        Returns:
            CAI score (0-1, higher is better)
        """
        codon_usage = self.get_codon_usage(sequence)
        
        if not codon_usage:
            return 0.0
        
        cai_values = []
        
        for codon, freq in codon_usage.items():
            if codon in self.GENETIC_CODE:
                aa = self.GENETIC_CODE[codon]
                if aa in self.HUMAN_CODON_USAGE and aa != '*':
                    # Get relative adaptiveness for this codon
                    aa_codons = self.HUMAN_CODON_USAGE[aa]
                    max_freq = max(aa_codons.values())
                    relative_adaptiveness = aa_codons.get(codon, 0) / max_freq
                    cai_values.append(relative_adaptiveness)
        
        if not cai_values:
            return 0.0
        
        # Geometric mean of relative adaptiveness values
        return np.exp(np.mean(np.log(np.maximum(cai_values, 1e-10))))
    
    def optimize_codon_usage(self, aa_sequence: str) -> str:
        """
        Optimize DNA sequence for human codon usage.
        
        Args:
            aa_sequence: Amino acid sequence
            
        Returns:
            Codon-optimized DNA sequence
        """
        optimized_codons = []
        
        for aa in aa_sequence:
            if aa == '*':
                # Use most common stop codon
                optimized_codons.append('TGA')
            elif aa in self.HUMAN_CODON_USAGE:
                # Choose most frequent codon for this amino acid
                aa_codons = self.HUMAN_CODON_USAGE[aa]
                best_codon = max(aa_codons.keys(), key=aa_codons.get)
                optimized_codons.append(best_codon)
            else:
                self.logger.warning(f"Unknown amino acid: {aa}")
                optimized_codons.append('NNN')
        
        return ''.join(optimized_codons)
    
    def find_alternative_codons(self, codon: str) -> List[Tuple[str, float]]:
        """
        Find alternative codons for the same amino acid.
        
        Args:
            codon: Original codon
            
        Returns:
            List of (codon, frequency) tuples
        """
        if codon not in self.GENETIC_CODE:
            return []
        
        aa = self.GENETIC_CODE[codon]
        if aa not in self.HUMAN_CODON_USAGE:
            return []
        
        aa_codons = self.HUMAN_CODON_USAGE[aa]
        return [(cod, freq) for cod, freq in aa_codons.items()]


class VariantAnalyzer:
    """Analyzes genetic variants and their effects."""
    
    def __init__(self):
        """Initialize variant analyzer."""
        self.codon_analyzer = CodonAnalyzer()
        self.logger = logging.getLogger(__name__)
    
    def find_variant_position(self, sequence: SeqRecord, variant: VariantInfo) -> Optional[VariantPosition]:
        """
        Find the position of a variant in a sequence.
        
        Args:
            sequence: Gene sequence
            variant: Variant information
            
        Returns:
            VariantPosition object or None if not found
        """
        # This is a simplified implementation
        # In practice, you'd need more sophisticated variant mapping
        
        if not variant.position:
            self.logger.warning(f"No position information for variant {variant.name}")
            return None
        
        try:
            # Convert genomic position to sequence position (simplified)
            seq_position = variant.position  # This would need proper coordinate conversion
            
            if seq_position >= len(sequence.seq):
                self.logger.warning(f"Variant position {seq_position} outside sequence length {len(sequence.seq)}")
                return None
            
            # Find coding sequence (simplified - look for first ATG)
            coding_start = str(sequence.seq).find('ATG')
            if coding_start == -1:
                self.logger.warning("No start codon found in sequence")
                return None
            
            # Calculate position relative to coding sequence
            if seq_position < coding_start:
                self.logger.warning("Variant in non-coding region")
                return None
            
            coding_position = seq_position - coding_start
            codon_number = coding_position // 3
            codon_position = coding_position % 3
            
            # Get codons
            codon_start = coding_start + (codon_number * 3)
            ref_codon = str(sequence.seq[codon_start:codon_start + 3])
            
            # Create alternative codon
            alt_codon = list(ref_codon)
            if codon_position < len(alt_codon):
                alt_codon[codon_position] = variant.alt_allele
            alt_codon = ''.join(alt_codon)
            
            # Translate codons
            ref_aa = self.codon_analyzer.GENETIC_CODE.get(ref_codon, 'X')
            alt_aa = self.codon_analyzer.GENETIC_CODE.get(alt_codon, 'X')
            
            return VariantPosition(
                gene_position=seq_position,
                codon_position=codon_position,
                codon_number=codon_number,
                amino_acid_position=codon_number + 1,
                reference_codon=ref_codon,
                alternative_codon=alt_codon,
                reference_aa=ref_aa,
                alternative_aa=alt_aa
            )
            
        except Exception as e:
            self.logger.error(f"Failed to find variant position: {e}")
            return None
    
    def predict_protein_effect(self, variant_pos: VariantPosition) -> ProteinEffect:
        """
        Predict the effect of a variant on protein function.
        
        Args:
            variant_pos: Variant position information
            
        Returns:
            ProteinEffect object with predictions
        """
        # Simplified protein effect prediction
        # In practice, you'd use tools like PolyPhen, SIFT, etc.
        
        stability_change = 0.0
        activity_change = 1.0
        structure_disrupted = False
        confidence = 0.5
        
        # Basic rules for effect prediction
        if variant_pos.is_synonymous:
            # Synonymous variants typically have minimal effect
            stability_change = np.random.normal(0, 0.2)
            activity_change = np.random.normal(1.0, 0.05)
            confidence = 0.8
        
        elif variant_pos.variant_type == "nonsense":
            # Stop codon variants are typically deleterious
            stability_change = -3.0
            activity_change = 0.0
            structure_disrupted = True
            confidence = 0.9
        
        elif variant_pos.variant_type == "missense":
            # Missense variants have variable effects
            # Consider amino acid properties
            ref_aa = variant_pos.reference_aa
            alt_aa = variant_pos.alternative_aa
            
            # Simplified amino acid property scoring
            aa_properties = {
                'hydrophobic': set('AILMFPWV'),
                'polar': set('NQST'),
                'charged': set('DEKR'),
                'aromatic': set('FWY'),
                'small': set('AGSV')
            }
            
            property_change = 0
            for prop, aa_set in aa_properties.items():
                if (ref_aa in aa_set) != (alt_aa in aa_set):
                    property_change += 1
            
            # More property changes = more severe effect
            stability_change = np.random.normal(-property_change * 0.5, 0.5)
            activity_change = max(0.1, 1.0 - property_change * 0.2)
            
            if property_change >= 3:
                structure_disrupted = True
            
            confidence = 0.6 + property_change * 0.1
        
        return ProteinEffect(
            stability_change=stability_change,
            activity_change=activity_change,
            structure_disrupted=structure_disrupted,
            confidence_score=min(1.0, confidence)
        )


class SequenceAnalyzer:
    """Main sequence analysis class."""
    
    def __init__(self):
        """Initialize sequence analyzer."""
        self.codon_analyzer = CodonAnalyzer()
        self.variant_analyzer = VariantAnalyzer()
        self.logger = logging.getLogger(__name__)
    
    def find_coding_sequence(self, sequence: SeqRecord) -> Optional[CodingRegion]:
        """
        Identify coding sequence within a gene sequence.
        
        Args:
            sequence: Gene sequence record
            
        Returns:
            CodingRegion object or None if not found
        """
        seq_str = str(sequence.seq).upper()
        
        # Look for start codon (ATG)
        start_positions = []
        for i, base in enumerate(seq_str):
            if seq_str[i:i+3] == 'ATG':
                start_positions.append(i)
        
        if not start_positions:
            self.logger.warning("No start codon found")
            return None
        
        # For each start position, look for stop codon in frame
        for start_pos in start_positions:
            for i in range(start_pos, len(seq_str) - 2, 3):
                codon = seq_str[i:i+3]
                if codon in ['TAA', 'TAG', 'TGA']:
                    # Found stop codon
                    end_pos = i + 2
                    
                    # Find exons (simplified - assume no introns for now)
                    exons = [(start_pos, end_pos)]
                    
                    return CodingRegion(
                        start=start_pos,
                        end=end_pos,
                        exons=exons,
                        introns=[]
                    )
        
        self.logger.warning("No complete coding sequence found")
        return None
    
    def analyze_sequence_composition(self, sequence: Union[str, Seq]) -> Dict[str, Any]:
        """
        Analyze basic sequence composition.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with composition analysis
        """
        if isinstance(sequence, str):
            seq_obj = Seq(sequence)
        else:
            seq_obj = sequence
        
        seq_str = str(seq_obj).upper()
        
        # Basic composition
        composition = Counter(seq_str)
        total_bases = len(seq_str)
        
        if total_bases == 0:
            return {}
        
        # Calculate percentages
        composition_pct = {base: count / total_bases * 100 
                          for base, count in composition.items()}
        
        # GC content
        gc_content = GC(seq_obj)
        
        # Molecular weight
        try:
            mol_weight = molecular_weight(seq_obj, seq_type='DNA')
        except:
            mol_weight = 0
        
        # Find repeats
        repeats = self._find_repeats(seq_str)
        
        return {
            'length': total_bases,
            'composition': composition_pct,
            'gc_content': gc_content,
            'molecular_weight': mol_weight,
            'repeats': repeats,
            'complexity_score': self._calculate_complexity(seq_str)
        }
    
    def _find_repeats(self, sequence: str, min_length: int = 4, 
                     max_length: int = 20) -> List[Dict[str, Any]]:
        """Find repetitive elements in sequence."""
        repeats = []
        
        for repeat_len in range(min_length, max_length + 1):
            for i in range(len(sequence) - repeat_len + 1):
                repeat_unit = sequence[i:i + repeat_len]
                
                # Count consecutive occurrences
                count = 1
                pos = i + repeat_len
                while pos + repeat_len <= len(sequence):
                    if sequence[pos:pos + repeat_len] == repeat_unit:
                        count += 1
                        pos += repeat_len
                    else:
                        break
                
                if count >= 3:  # At least 3 repeats
                    repeats.append({
                        'sequence': repeat_unit,
                        'position': i,
                        'length': repeat_len,
                        'count': count,
                        'total_length': repeat_len * count
                    })
        
        # Remove overlapping repeats (keep longest)
        repeats.sort(key=lambda x: x['total_length'], reverse=True)
        non_overlapping = []
        used_positions = set()
        
        for repeat in repeats:
            start = repeat['position']
            end = start + repeat['total_length']
            
            if not any(pos in used_positions for pos in range(start, end)):
                non_overlapping.append(repeat)
                used_positions.update(range(start, end))
        
        return non_overlapping
    
    def _calculate_complexity(self, sequence: str, window_size: int = 50) -> float:
        """Calculate sequence complexity using Shannon entropy."""
        if len(sequence) < window_size:
            window_size = len(sequence)
        
        if window_size == 0:
            return 0.0
        
        complexities = []
        
        for i in range(0, len(sequence) - window_size + 1, window_size // 2):
            window = sequence[i:i + window_size]
            base_counts = Counter(window)
            
            # Calculate Shannon entropy
            entropy = 0
            for count in base_counts.values():
                if count > 0:
                    p = count / len(window)
                    entropy -= p * np.log2(p)
            
            # Normalize by maximum possible entropy
            max_entropy = np.log2(min(4, len(set(window))))
            complexity = entropy / max_entropy if max_entropy > 0 else 0
            complexities.append(complexity)
        
        return np.mean(complexities) if complexities else 0.0
    
    def find_regulatory_elements(self, sequence: Union[str, Seq]) -> Dict[str, List[Dict]]:
        """
        Find potential regulatory elements in sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary of regulatory element types and positions
        """
        if isinstance(sequence, Seq):
            sequence = str(sequence)
        
        sequence = sequence.upper()
        
        # Common regulatory motifs (simplified patterns)
        motifs = {
            'TATA_box': r'TATAAA',
            'CAAT_box': r'CCAAT',
            'GC_box': r'GGGCGG',
            'initiator': r'[CT][CT]A[ACGT][ACGT][CT][CT]',
            'CpG_island': r'CG.{0,20}CG.{0,20}CG',  # Simplified CpG pattern
            'kozak_consensus': r'[AG]CCATGG',  # Start codon context
            'poly_A_signal': r'AATAAA'
        }
        
        regulatory_elements = {}
        
        for element_type, pattern in motifs.items():
            matches = []
            for match in re.finditer(pattern, sequence):
                matches.append({
                    'position': match.start(),
                    'sequence': match.group(),
                    'score': 1.0  # Simplified scoring
                })
            regulatory_elements[element_type] = matches
        
        return regulatory_elements
    
    def predict_secondary_structure(self, sequence: Union[str, Seq]) -> Dict[str, Any]:
        """
        Predict DNA secondary structure (simplified).
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with structure predictions
        """
        if isinstance(sequence, Seq):
            sequence = str(sequence)
        
        sequence = sequence.upper()
        
        # Find potential hairpin structures
        hairpins = self._find_hairpins(sequence)
        
        # Find inverted repeats
        inverted_repeats = self._find_inverted_repeats(sequence)
        
        return {
            'hairpins': hairpins,
            'inverted_repeats': inverted_repeats,
            'structure_complexity': len(hairpins) + len(inverted_repeats)
        }
    
    def _find_hairpins(self, sequence: str, min_stem: int = 6, 
                      max_loop: int = 10) -> List[Dict[str, Any]]:
        """Find potential hairpin structures."""
        hairpins = []
        
        for i in range(len(sequence) - min_stem * 2 - max_loop):
            for stem_len in range(min_stem, min(20, (len(sequence) - i) // 2)):
                for loop_len in range(3, max_loop + 1):
                    left_stem = sequence[i:i + stem_len]
                    loop_start = i + stem_len
                    loop_end = loop_start + loop_len
                    right_stem = sequence[loop_end:loop_end + stem_len]
                    
                    if loop_end + stem_len > len(sequence):
                        break
                    
                    # Check if stems are complementary
                    complement = str(Seq(left_stem).reverse_complement())
                    
                    # Calculate complementarity
                    matches = sum(1 for a, b in zip(complement, right_stem) if a == b)
                    complementarity = matches / len(complement)
                    
                    if complementarity >= 0.8:  # At least 80% complementary
                        hairpins.append({
                            'position': i,
                            'stem_length': stem_len,
                            'loop_length': loop_len,
                            'complementarity': complementarity,
                            'stability_score': complementarity * stem_len
                        })
        
        return hairpins
    
    def _find_inverted_repeats(self, sequence: str, min_length: int = 8) -> List[Dict[str, Any]]:
        """Find inverted repeat sequences."""
        inverted_repeats = []
        
        for i in range(len(sequence) - min_length + 1):
            for length in range(min_length, min(50, len(sequence) - i + 1)):
                subseq = sequence[i:i + length]
                rev_comp = str(Seq(subseq).reverse_complement())
                
                # Look for this reverse complement elsewhere in sequence
                for j in range(i + length, len(sequence) - length + 1):
                    target = sequence[j:j + length]
                    
                    if target == rev_comp:
                        inverted_repeats.append({
                            'position1': i,
                            'position2': j,
                            'length': length,
                            'sequence': subseq,
                            'distance': j - i - length
                        })
        
        return inverted_repeats
    
    def compare_sequences(self, seq1: Union[str, Seq], seq2: Union[str, Seq]) -> Dict[str, Any]:
        """
        Compare two sequences for similarity.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Comparison metrics
        """
        if isinstance(seq1, Seq):
            seq1 = str(seq1)
        if isinstance(seq2, Seq):
            seq2 = str(seq2)
        
        seq1, seq2 = seq1.upper(), seq2.upper()
        
        # Basic alignment (simplified)
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        
        if min_len == 0:
            return {'identity': 0.0, 'similarity': 0.0}
        
        # Count matches in overlapping region
        matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
        
        identity = matches / min_len
        similarity = matches / max_len  # Account for length differences
        
        return {
            'identity': identity,
            'similarity': similarity,
            'matches': matches,
            'length_diff': abs(len(seq1) - len(seq2)),
            'alignment_score': identity * min_len
        }