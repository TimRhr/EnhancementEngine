"""
Disease risk calculation and assessment module.

This module provides comprehensive disease risk calculation based on genetic
variants, including population genetics, penetrance, and personalized risk scores.
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
import math

from ..models.data_classes import GeneInfo, VariantInfo, RiskLevel
from ..models.therapeutic_data_classes import (
    DiseaseRisk, DiseaseCategory, PatientStratification
)
from .disease_db import DiseaseDatabaseClient
from ..models.disease_constants import (
    DISEASE_GENES, POPULATION_FREQUENCIES, DISEASE_SCORING_SYSTEMS
)


class RiskCalculationError(Exception):
    """Custom exception for risk calculation errors."""
    pass


@dataclass
class PopulationRiskProfile:
    """Population-specific risk profile."""
    ethnicity: str
    baseline_prevalence: float
    age_adjusted_risk: Dict[str, float]  # Age group -> risk multiplier
    sex_specific_risk: Dict[str, float]  # Sex -> risk multiplier
    environmental_factors: Dict[str, float] = None
    
    def __post_init__(self):
        if self.environmental_factors is None:
            self.environmental_factors = {}


@dataclass
class GeneticRiskScore:
    """Comprehensive genetic risk score."""
    total_score: float
    individual_contributions: Dict[str, float]  # Gene -> contribution
    confidence_interval: Tuple[float, float]
    risk_category: str
    lifetime_risk: float
    age_specific_risks: Dict[int, float]  # Age -> risk


class PopulationGeneticsCalculator:
    """Calculates population-level genetic risk parameters."""
    
    def __init__(self):
        """Initialize population genetics calculator."""
        self.logger = logging.getLogger(__name__)
        
        # Load population-specific data
        self.population_baselines = self._load_population_baselines()
        self.hardy_weinberg_calculator = HardyWeinbergCalculator()
    
    def _load_population_baselines(self) -> Dict[str, PopulationRiskProfile]:
        """Load baseline risk profiles for different populations."""
        return {
            "european": PopulationRiskProfile(
                ethnicity="european",
                baseline_prevalence=0.008,  # 0.8% for RA
                age_adjusted_risk={
                    "20-30": 0.3,
                    "30-40": 0.7,
                    "40-50": 1.0,
                    "50-60": 1.5,
                    "60-70": 2.0,
                    "70+": 2.5
                },
                sex_specific_risk={
                    "female": 2.5,  # Females 2.5x higher risk
                    "male": 1.0
                }
            ),
            "east_asian": PopulationRiskProfile(
                ethnicity="east_asian",
                baseline_prevalence=0.003,  # 0.3% for RA
                age_adjusted_risk={
                    "20-30": 0.2,
                    "30-40": 0.6,
                    "40-50": 1.0,
                    "50-60": 1.3,
                    "60-70": 1.8,
                    "70+": 2.2
                },
                sex_specific_risk={
                    "female": 3.0,
                    "male": 1.0
                }
            ),
            "african": PopulationRiskProfile(
                ethnicity="african",
                baseline_prevalence=0.005,  # 0.5% for RA
                age_adjusted_risk={
                    "20-30": 0.4,
                    "30-40": 0.8,
                    "40-50": 1.0,
                    "50-60": 1.4,
                    "60-70": 1.9,
                    "70+": 2.3
                },
                sex_specific_risk={
                    "female": 2.2,
                    "male": 1.0
                }
            )
        }
    
    def calculate_population_attributable_risk(self, gene: str, variant: str,
                                             population: str = "european") -> float:
        """
        Calculate population attributable risk for a genetic variant.
        
        Args:
            gene: Gene symbol
            variant: Variant identifier
            population: Population group
            
        Returns:
            Population attributable risk (0-1)
        """
        try:
            # Get variant information
            gene_data = DISEASE_GENES.get("autoimmune", {}).get(gene)
            if not gene_data:
                return 0.0
            
            variant_data = gene_data.get("pathogenic_variants", {}).get(variant)
            if not variant_data:
                return 0.0
            
            # Get population frequency
            freq_key = f"{gene}_{variant}"
            allele_freq = POPULATION_FREQUENCIES.get(population, {}).get(freq_key, 0.0)
            
            if allele_freq == 0:
                return 0.0
            
            # Get odds ratio from disease associations
            disease_data = gene_data.get("disease_associations", {}).get("rheumatoid_arthritis", {})
            odds_ratio = disease_data.get("odds_ratio", 1.0)
            
            # Calculate PAR using standard formula: p(OR-1) / [1 + p(OR-1)]
            # where p = allele frequency, OR = odds ratio
            if odds_ratio <= 1.0:
                return 0.0
            
            par = (allele_freq * (odds_ratio - 1)) / (1 + allele_freq * (odds_ratio - 1))
            
            self.logger.debug(f"PAR for {gene} {variant} in {population}: {par:.3f}")
            return par
            
        except Exception as e:
            self.logger.error(f"PAR calculation failed: {e}")
            return 0.0
    
    def estimate_penetrance(self, gene: str, variant: str, age: int, sex: str,
                          population: str = "european") -> float:
        """
        Estimate age and sex-specific penetrance.
        
        Args:
            gene: Gene symbol
            variant: Variant identifier
            age: Patient age
            sex: Patient sex
            population: Population group
            
        Returns:
            Penetrance estimate (0-1)
        """
        try:
            # Get base penetrance
            gene_data = DISEASE_GENES.get("autoimmune", {}).get(gene, {})
            disease_data = gene_data.get("disease_associations", {}).get("rheumatoid_arthritis", {})
            base_penetrance = disease_data.get("penetrance", 0.1)
            
            # Get population profile
            pop_profile = self.population_baselines.get(population)
            if not pop_profile:
                return base_penetrance
            
            # Age adjustment
            age_group = self._get_age_group(age)
            age_multiplier = pop_profile.age_adjusted_risk.get(age_group, 1.0)
            
            # Sex adjustment
            sex_multiplier = pop_profile.sex_specific_risk.get(sex.lower(), 1.0)
            
            # Calculate adjusted penetrance
            adjusted_penetrance = base_penetrance * age_multiplier * sex_multiplier
            
            # Cap at 1.0
            return min(1.0, adjusted_penetrance)
            
        except Exception as e:
            self.logger.error(f"Penetrance estimation failed: {e}")
            return 0.1  # Default low penetrance
    
    def _get_age_group(self, age: int) -> str:
        """Categorize age into groups."""
        if age < 30:
            return "20-30"
        elif age < 40:
            return "30-40"
        elif age < 50:
            return "40-50"
        elif age < 60:
            return "50-60"
        elif age < 70:
            return "60-70"
        else:
            return "70+"


class HardyWeinbergCalculator:
    """Calculates Hardy-Weinberg equilibrium and deviations."""
    
    def __init__(self):
        """Initialize Hardy-Weinberg calculator."""
        self.logger = logging.getLogger(__name__)
    
    def calculate_genotype_frequencies(self, allele_freq: float) -> Dict[str, float]:
        """
        Calculate expected genotype frequencies under Hardy-Weinberg equilibrium.
        
        Args:
            allele_freq: Frequency of variant allele (q)
            
        Returns:
            Dictionary of genotype frequencies
        """
        p = 1 - allele_freq  # Wild-type allele frequency
        q = allele_freq      # Variant allele frequency
        
        return {
            "homozygous_wildtype": p * p,
            "heterozygous": 2 * p * q,
            "homozygous_variant": q * q
        }
    
    def test_hardy_weinberg_deviation(self, observed_genotypes: Dict[str, int],
                                    expected_freq: Dict[str, float]) -> Dict[str, Any]:
        """
        Test for deviation from Hardy-Weinberg equilibrium.
        
        Args:
            observed_genotypes: Observed genotype counts
            expected_freq: Expected genotype frequencies
            
        Returns:
            Test results including chi-square statistic and p-value
        """
        try:
            # Calculate total count
            total = sum(observed_genotypes.values())
            
            # Calculate expected counts
            expected_counts = {
                genotype: freq * total 
                for genotype, freq in expected_freq.items()
            }
            
            # Chi-square test
            chi_square = 0
            for genotype in observed_genotypes:
                observed = observed_genotypes[genotype]
                expected = expected_counts.get(genotype, 0)
                
                if expected > 0:
                    chi_square += ((observed - expected) ** 2) / expected
            
            # Degrees of freedom = genotypes - alleles = 3 - 2 = 1
            df = 1
            
            # Simple p-value approximation (would use scipy.stats.chi2 in practice)
            p_value = math.exp(-chi_square / 2) if chi_square < 10 else 0.001
            
            return {
                "chi_square": chi_square,
                "p_value": p_value,
                "degrees_freedom": df,
                "significant_deviation": p_value < 0.05,
                "expected_counts": expected_counts,
                "observed_counts": observed_genotypes
            }
            
        except Exception as e:
            self.logger.error(f"Hardy-Weinberg test failed: {e}")
            return {"error": str(e)}


class MultiGeneRiskCalculator:
    """Calculates combined risk from multiple genetic variants."""
    
    def __init__(self):
        """Initialize multi-gene risk calculator."""
        self.logger = logging.getLogger(__name__)
        self.population_calculator = PopulationGeneticsCalculator()
    
    def calculate_polygenic_risk_score(self, genotype_data: Dict[str, str],
                                     patient_data: Dict[str, Any],
                                     disease: str = "rheumatoid_arthritis") -> GeneticRiskScore:
        """
        Calculate polygenic risk score from multiple variants.
        
        Args:
            genotype_data: Gene -> genotype mapping
            patient_data: Patient demographics and characteristics
            disease: Target disease
            
        Returns:
            GeneticRiskScore object
        """
        try:
            individual_contributions = {}
            total_log_odds = 0
            
            # Process each gene
            for gene, genotype in genotype_data.items():
                gene_contribution = self._calculate_gene_contribution(
                    gene, genotype, disease, patient_data
                )
                
                if gene_contribution is not None:
                    individual_contributions[gene] = gene_contribution
                    total_log_odds += math.log(gene_contribution) if gene_contribution > 0 else 0
            
            # Convert log odds to probability
            total_odds = math.exp(total_log_odds)
            risk_score = total_odds / (1 + total_odds)
            
            # Calculate confidence interval
            confidence_interval = self._calculate_confidence_interval(
                individual_contributions, patient_data
            )
            
            # Categorize risk
            risk_category = self._categorize_risk(risk_score)
            
            # Calculate lifetime risk
            population = patient_data.get("ethnicity", "european")
            baseline_risk = self.population_calculator.population_baselines[population].baseline_prevalence
            lifetime_risk = min(1.0, baseline_risk * total_odds)
            
            # Calculate age-specific risks
            age_specific_risks = self._calculate_age_specific_risks(
                total_odds, patient_data, population
            )
            
            return GeneticRiskScore(
                total_score=risk_score,
                individual_contributions=individual_contributions,
                confidence_interval=confidence_interval,
                risk_category=risk_category,
                lifetime_risk=lifetime_risk,
                age_specific_risks=age_specific_risks
            )
            
        except Exception as e:
            self.logger.error(f"Polygenic risk score calculation failed: {e}")
            raise RiskCalculationError(f"Failed to calculate polygenic risk: {e}")
    
    def _calculate_gene_contribution(self, gene: str, genotype: str, disease: str,
                                   patient_data: Dict[str, Any]) -> Optional[float]:
        """Calculate individual gene contribution to risk."""
        try:
            # Get gene data
            gene_data = None
            for category in DISEASE_GENES.values():
                if gene in category:
                    gene_data = category[gene]
                    break
            
            if not gene_data:
                return None
            
            # Get disease association
            disease_data = gene_data.get("disease_associations", {}).get(disease, {})
            if not disease_data:
                return None
            
            # Determine risk based on genotype
            odds_ratio = disease_data.get("odds_ratio", 1.0)
            
            # Adjust odds ratio based on genotype
            if "heterozygous" in genotype.lower() or "het" in genotype.lower():
                # Heterozygous - intermediate risk
                adjusted_or = 1 + (odds_ratio - 1) * 0.6
            elif "homozygous_variant" in genotype.lower() or "hom" in genotype.lower():
                # Homozygous variant - full risk
                adjusted_or = odds_ratio
            else:
                # Wild-type - no increased risk
                adjusted_or = 1.0
            
            # Adjust for population and demographics
            population = patient_data.get("ethnicity", "european")
            age = patient_data.get("age", 40)
            sex = patient_data.get("sex", "female")
            
            penetrance = self.population_calculator.estimate_penetrance(
                gene, genotype, age, sex, population
            )
            
            # Final contribution is odds ratio weighted by penetrance
            return adjusted_or * penetrance
            
        except Exception as e:
            self.logger.error(f"Gene contribution calculation failed for {gene}: {e}")
            return None
    
    def _calculate_confidence_interval(self, contributions: Dict[str, float],
                                     patient_data: Dict[str, Any]) -> Tuple[float, float]:
        """Calculate confidence interval for risk score."""
        # Simplified confidence interval calculation
        # In practice, would use more sophisticated methods
        
        if not contributions:
            return (0.0, 0.0)
        
        mean_contribution = np.mean(list(contributions.values()))
        std_contribution = np.std(list(contributions.values()))
        
        # 95% confidence interval
        margin_error = 1.96 * std_contribution / math.sqrt(len(contributions))
        
        lower_bound = max(0.0, mean_contribution - margin_error)
        upper_bound = min(1.0, mean_contribution + margin_error)
        
        return (lower_bound, upper_bound)
    
    def _categorize_risk(self, risk_score: float) -> str:
        """Categorize risk score into clinical categories."""
        if risk_score < 0.05:
            return "Very Low Risk"
        elif risk_score < 0.1:
            return "Low Risk"
        elif risk_score < 0.2:
            return "Moderate Risk"
        elif risk_score < 0.4:
            return "High Risk"
        else:
            return "Very High Risk"
    
    def _calculate_age_specific_risks(self, total_odds: float, patient_data: Dict[str, Any],
                                    population: str) -> Dict[int, float]:
        """Calculate risk at different ages."""
        baseline_risk = self.population_calculator.population_baselines[population].baseline_prevalence
        age_risks = {}
        
        # Calculate risk for ages 20-80
        for age in range(20, 81, 10):
            # Age-specific adjustment
            age_group = self.population_calculator._get_age_group(age)
            pop_profile = self.population_calculator.population_baselines[population]
            age_multiplier = pop_profile.age_adjusted_risk.get(age_group, 1.0)
            
            # Age-specific risk
            age_specific_risk = min(1.0, baseline_risk * total_odds * age_multiplier)
            age_risks[age] = age_specific_risk
        
        return age_risks


class DiseaseRiskCalculator:
    """Main disease risk calculator integrating all components."""

    def __init__(self, disease_client: Optional["DiseaseDatabaseClient"] = None):
        """Initialize disease risk calculator."""
        self.population_calculator = PopulationGeneticsCalculator()
        self.multigene_calculator = MultiGeneRiskCalculator()
        self.hardy_weinberg = HardyWeinbergCalculator()
        self.disease_client = disease_client
        self.logger = logging.getLogger(__name__)
    
    def calculate_disease_risk(self, gene: str, variant: str, patient_data: Dict[str, Any],
                             disease: str = "rheumatoid_arthritis") -> DiseaseRisk:
        """
        Calculate comprehensive disease risk for a single variant.
        
        Args:
            gene: Gene symbol
            variant: Variant identifier
            patient_data: Patient demographics
            disease: Target disease
            
        Returns:
            DiseaseRisk object
        """
        try:
            # Get variant data from constants
            gene_data = None
            for category in DISEASE_GENES.values():
                if gene in category:
                    gene_data = category[gene]
                    break

            disease_data = {}
            if gene_data:
                disease_data = gene_data.get("disease_associations", {}).get(disease, {})

            # Fallback to dynamic client if no data
            if not gene_data or not disease_data:
                if self.disease_client:
                    dynamic = self.disease_client.fetch_disease_info(gene, variant, disease)
                    if dynamic:
                        odds_ratio = dynamic.get("odds_ratio", 1.0)
                        population_frequency = dynamic.get("allele_frequency", 0.0)
                        penetrance = dynamic.get("penetrance", 0.1)
                    else:
                        raise RiskCalculationError(f"No data for {gene} {variant}")
                else:
                    raise RiskCalculationError(f"Gene {gene} not found in disease database")
            else:
                odds_ratio = disease_data.get("odds_ratio", 1.0)
                population_frequency = disease_data.get("population_frequency", 0.0)
                penetrance = disease_data.get("penetrance", 0.1)
            
            # Adjust for patient-specific factors
            population = patient_data.get("ethnicity", "european")
            age = patient_data.get("age", 40)
            sex = patient_data.get("sex", "female")
            
            # Calculate adjusted penetrance
            adjusted_penetrance = self.population_calculator.estimate_penetrance(
                gene, variant, age, sex, population
            )
            
            # Calculate absolute risk increase
            baseline_risk = self.population_calculator.population_baselines[population].baseline_prevalence
            absolute_risk_increase = baseline_risk * (odds_ratio - 1) * adjusted_penetrance
            
            # Get population-specific frequencies
            ethnicity_frequencies = {}
            for pop in POPULATION_FREQUENCIES:
                freq_key = f"{gene}_{variant}"
                freq = POPULATION_FREQUENCIES[pop].get(freq_key, 0.0)
                ethnicity_frequencies[pop] = freq
            
            # Calculate confidence interval (simplified)
            ci_lower = odds_ratio * 0.8
            ci_upper = odds_ratio * 1.2
            
            disease_risk = DiseaseRisk(
                disease_name=disease,
                odds_ratio=odds_ratio,
                population_frequency=population_frequency,
                absolute_risk_increase=absolute_risk_increase,
                confidence_interval=(ci_lower, ci_upper),
                ethnicity_specific=ethnicity_frequencies,
                age_dependent=True,
                sex_specific=(sex.lower() == "female"),
                penetrance=adjusted_penetrance
            )
            
            self.logger.info(f"Calculated {disease} risk for {gene} {variant}: OR={odds_ratio:.2f}")
            return disease_risk
            
        except Exception as e:
            self.logger.error(f"Disease risk calculation failed: {e}")
            raise RiskCalculationError(f"Failed to calculate disease risk: {e}")
    
    def calculate_combined_risk(self, genotype_data: Dict[str, str], patient_data: Dict[str, Any],
                              disease: str = "rheumatoid_arthritis") -> GeneticRiskScore:
        """
        Calculate combined risk from multiple genetic variants.
        
        Args:
            genotype_data: Gene -> genotype mapping
            patient_data: Patient demographics
            disease: Target disease
            
        Returns:
            GeneticRiskScore object
        """
        return self.multigene_calculator.calculate_polygenic_risk_score(
            genotype_data, patient_data, disease
        )
    
    def simulate_intervention_benefit(self, baseline_risk: DiseaseRisk,
                                    correction_efficiency: float) -> DiseaseRisk:
        """
        Simulate risk reduction from therapeutic intervention.
        
        Args:
            baseline_risk: Original disease risk
            correction_efficiency: Fraction of risk alleles corrected (0-1)
            
        Returns:
            Reduced disease risk after intervention
        """
        try:
            # Calculate residual risk after correction
            residual_risk_fraction = 1 - correction_efficiency
            
            # New odds ratio accounting for correction
            corrected_or = 1 + (baseline_risk.odds_ratio - 1) * residual_risk_fraction
            
            # New absolute risk increase
            corrected_ari = baseline_risk.absolute_risk_increase * residual_risk_fraction
            
            # Create corrected risk object
            corrected_risk = DiseaseRisk(
                disease_name=f"{baseline_risk.disease_name}_post_intervention",
                odds_ratio=corrected_or,
                population_frequency=baseline_risk.population_frequency,
                absolute_risk_increase=corrected_ari,
                confidence_interval=(
                    baseline_risk.confidence_interval[0] * residual_risk_fraction,
                    baseline_risk.confidence_interval[1] * residual_risk_fraction
                ),
                ethnicity_specific=baseline_risk.ethnicity_specific,
                age_dependent=baseline_risk.age_dependent,
                sex_specific=baseline_risk.sex_specific,
                penetrance=baseline_risk.penetrance * residual_risk_fraction
            )
            
            self.logger.info(f"Simulated intervention: OR {baseline_risk.odds_ratio:.2f} -> {corrected_or:.2f}")
            return corrected_risk
            
        except Exception as e:
            self.logger.error(f"Intervention simulation failed: {e}")
            raise RiskCalculationError(f"Failed to simulate intervention: {e}")
    
    def validate_population_frequencies(self, population: str, observed_data: Dict[str, Dict[str, int]]) -> Dict[str, Any]:
        """
        Validate population frequencies against observed data.
        
        Args:
            population: Population group
            observed_data: Gene -> {genotype: count} mapping
            
        Returns:
            Validation results
        """
        validation_results = {}
        
        for gene, genotype_counts in observed_data.items():
            try:
                # Get expected frequency for this gene's main variant
                main_variant = list(DISEASE_GENES.get("autoimmune", {}).get(gene, {}).get("pathogenic_variants", {}).keys())[0]
                freq_key = f"{gene}_{main_variant}"
                expected_allele_freq = POPULATION_FREQUENCIES.get(population, {}).get(freq_key, 0.0)
                
                if expected_allele_freq > 0:
                    # Calculate expected genotype frequencies
                    expected_freq = self.hardy_weinberg.calculate_genotype_frequencies(expected_allele_freq)
                    
                    # Test for Hardy-Weinberg deviation
                    hw_test = self.hardy_weinberg.test_hardy_weinberg_deviation(
                        genotype_counts, expected_freq
                    )
                    
                    validation_results[gene] = {
                        "expected_allele_frequency": expected_allele_freq,
                        "expected_genotype_frequencies": expected_freq,
                        "hardy_weinberg_test": hw_test,
                        "sample_size": sum(genotype_counts.values())
                    }
                
            except Exception as e:
                self.logger.error(f"Validation failed for {gene}: {e}")
                validation_results[gene] = {"error": str(e)}
        
        return validation_results