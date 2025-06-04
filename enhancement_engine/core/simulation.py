"""
Effect simulation module for Enhancement Engine.

This module provides comprehensive simulation of enhancement effects including:
- Phenotype prediction from genetic modifications
- Enhancement gain calculation
- Population-level effect modeling
- Dose-response relationships
- Long-term outcome simulation
"""

import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass
from collections import defaultdict
import json

try:
    from scipy import stats
    from scipy.optimize import curve_fit
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    logging.warning("SciPy not available. Some advanced statistical features will be limited.")

from ..models.data_classes import (
    EnhancementGain, VariantEffect, SideEffect, EnhancementCategory,
    VariantInfo, ProteinEffect, GeneInfo
)
from ..models.constants import (
    ENHANCEMENT_GENES, DEFAULT_PARAMETERS, SCORING_THRESHOLDS
)


class SimulationError(Exception):
    """Custom exception for simulation errors."""
    pass


@dataclass
class PopulationModel:
    """Population parameters for effect simulation."""
    mean_baseline: float
    std_baseline: float
    genetic_background_variance: float
    environmental_variance: float
    age_effect_coefficient: float = 0.0
    sex_effect_coefficient: float = 0.0


@dataclass
class DoseResponseModel:
    """Dose-response relationship parameters."""
    ec50: float  # Half-maximal effective concentration
    hill_coefficient: float  # Hill slope
    max_effect: float  # Maximum possible effect
    baseline_effect: float  # Effect at zero dose


class PhenotypePredictor:
    """Predicts phenotypic outcomes from genetic modifications."""
    
    def __init__(self):
        """Initialize phenotype predictor."""
        self.logger = logging.getLogger(__name__)
        
        # Load phenotype models for enhancement genes
        self.phenotype_models = self._load_phenotype_models()
        
        # Population baseline values
        self.population_baselines = self._load_population_baselines()
    
    def _load_phenotype_models(self) -> Dict[str, Dict[str, Any]]:
        """Load phenotype prediction models for enhancement genes."""
        models = {
            'COMT': {
                'cognitive_performance': {
                    'val_val_baseline': 7.2,
                    'val_met_effect': 0.6,
                    'met_met_effect': 1.4,
                    'variance': 0.8,
                    'age_interaction': -0.02,  # Decreases with age
                    'stress_interaction': 0.3   # Better under stress
                },
                'working_memory': {
                    'val_val_baseline': 6.8,
                    'val_met_effect': 0.8,
                    'met_met_effect': 1.5,
                    'variance': 0.6
                },
                'stress_sensitivity': {
                    'val_val_baseline': 8.1,
                    'val_met_effect': -1.2,  # Negative = better (less sensitive)
                    'met_met_effect': -2.7,
                    'variance': 1.0
                }
            },
            'BDNF': {
                'learning_rate': {
                    'val_val_baseline': 6.5,
                    'val_met_effect': -0.4,  # Met allele reduces function
                    'met_met_effect': -0.8,
                    'variance': 0.7,
                    'age_interaction': -0.015
                },
                'memory_consolidation': {
                    'val_val_baseline': 7.3,
                    'val_met_effect': -0.6,
                    'met_met_effect': -1.2,
                    'variance': 0.8
                },
                'neuroplasticity': {
                    'val_val_baseline': 7.0,
                    'val_met_effect': -0.5,
                    'met_met_effect': -1.0,
                    'variance': 0.9
                }
            },
            'ACTN3': {
                'sprint_performance': {
                    'rr_baseline': 7.8,  # R/R genotype
                    'rx_effect': -0.8,   # R/X genotype
                    'xx_effect': -1.6,   # X/X genotype (no functional protein)
                    'variance': 1.2
                },
                'power_output': {
                    'rr_baseline': 8.2,
                    'rx_effect': -0.6,
                    'xx_effect': -1.2,
                    'variance': 1.0
                },
                'muscle_fiber_composition': {
                    'rr_baseline': 6.5,  # % fast-twitch fibers
                    'rx_effect': -1.0,
                    'xx_effect': -2.0,
                    'variance': 1.5
                }
            },
            'FOXO3': {
                'cellular_stress_resistance': {
                    'cc_baseline': 6.0,   # Common allele
                    'ct_effect': 0.8,     # Protective allele heterozygous
                    'tt_effect': 1.6,     # Protective allele homozygous
                    'variance': 0.9
                },
                'longevity_score': {
                    'cc_baseline': 75.0,  # Expected lifespan
                    'ct_effect': 3.2,
                    'tt_effect': 6.5,
                    'variance': 8.0
                },
                'oxidative_stress_resistance': {
                    'cc_baseline': 5.8,
                    'ct_effect': 1.0,
                    'tt_effect': 2.0,
                    'variance': 1.1
                }
            },
            'MSTN': {
                'muscle_mass': {
                    'wt_baseline': 6.5,   # Wild-type
                    'het_ko_effect': 1.5, # Heterozygous knockout
                    'hom_ko_effect': 3.0, # Homozygous knockout
                    'variance': 1.2
                },
                'strength': {
                    'wt_baseline': 7.0,
                    'het_ko_effect': 1.8,
                    'hom_ko_effect': 3.5,
                    'variance': 1.4
                },
                'endurance': {
                    'wt_baseline': 6.8,
                    'het_ko_effect': -0.5,  # Might reduce endurance
                    'hom_ko_effect': -1.0,
                    'variance': 1.0
                }
            }
        }
        
        return models
    
    def _load_population_baselines(self) -> Dict[EnhancementCategory, PopulationModel]:
        """Load population baseline parameters."""
        return {
            EnhancementCategory.COGNITIVE: PopulationModel(
                mean_baseline=100.0,  # IQ-like scale
                std_baseline=15.0,
                genetic_background_variance=8.0,
                environmental_variance=10.0,
                age_effect_coefficient=-0.1,  # Slight decline with age
                sex_effect_coefficient=0.0     # No systematic difference
            ),
            EnhancementCategory.PHYSICAL: PopulationModel(
                mean_baseline=50.0,   # Fitness percentile
                std_baseline=20.0,
                genetic_background_variance=12.0,
                environmental_variance=15.0,
                age_effect_coefficient=-0.3,  # Decline with age
                sex_effect_coefficient=5.0     # Males slightly higher on average
            ),
            EnhancementCategory.LONGEVITY: PopulationModel(
                mean_baseline=78.5,   # Life expectancy
                std_baseline=8.0,
                genetic_background_variance=6.0,
                environmental_variance=10.0,
                age_effect_coefficient=0.0,   # N/A for longevity
                sex_effect_coefficient=3.0    # Females live longer
            )
        }
    
    def predict_phenotype(self, gene_name: str, variant: str, 
                         individual_params: Optional[Dict[str, Any]] = None) -> Dict[str, float]:
        """
        Predict phenotype for a specific gene variant.
        
        Args:
            gene_name: Gene symbol
            variant: Variant name or genotype
            individual_params: Individual-specific parameters (age, sex, etc.)
            
        Returns:
            Dictionary of phenotype predictions
        """
        if gene_name not in self.phenotype_models:
            raise SimulationError(f"No phenotype model available for gene: {gene_name}")
        
        gene_model = self.phenotype_models[gene_name]
        predictions = {}
        
        # Default individual parameters
        params = {
            'age': 30,
            'sex': 'M',
            'stress_level': 5.0,
            'fitness_level': 5.0
        }
        if individual_params:
            params.update(individual_params)
        
        for phenotype, model in gene_model.items():
            prediction = self._calculate_phenotype_value(model, variant, params)
            predictions[phenotype] = prediction
        
        return predictions
    
    def _calculate_phenotype_value(self, model: Dict[str, float], variant: str, 
                                 params: Dict[str, Any]) -> float:
        """Calculate phenotype value for a specific model and variant."""
        # Get baseline and effect based on variant
        if 'val_val_baseline' in model:  # COMT-style model
            if variant.lower() in ['val158met', 'val_met', 'heterozygous']:
                baseline = model['val_val_baseline']
                effect = model['val_met_effect']
            elif variant.lower() in ['met158met', 'met_met', 'homozygous_variant']:
                baseline = model['val_val_baseline']
                effect = model['met_met_effect']
            else:  # val_val or wild-type
                baseline = model['val_val_baseline']
                effect = 0.0
                
        elif 'rr_baseline' in model:  # ACTN3-style model
            if variant.lower() in ['r577x', 'rx', 'heterozygous']:
                baseline = model['rr_baseline']
                effect = model['rx_effect']
            elif variant.lower() in ['x577x', 'xx', 'homozygous_variant']:
                baseline = model['rr_baseline']
                effect = model['xx_effect']
            else:  # RR or wild-type
                baseline = model['rr_baseline']
                effect = 0.0
                
        elif 'wt_baseline' in model:  # Knockout-style model (MSTN)
            if variant.lower() in ['heterozygous_ko', 'het_ko', '+/-']:
                baseline = model['wt_baseline']
                effect = model['het_ko_effect']
            elif variant.lower() in ['homozygous_ko', 'hom_ko', '-/-']:
                baseline = model['wt_baseline']
                effect = model['hom_ko_effect']
            else:  # Wild-type
                baseline = model['wt_baseline']
                effect = 0.0
                
        else:  # Generic model
            baseline = model.get('baseline', 5.0)
            effect = model.get('effect', 0.0)
        
        # Apply interactions
        total_effect = effect
        
        if 'age_interaction' in model:
            age_effect = model['age_interaction'] * (params['age'] - 30)  # Relative to age 30
            total_effect += age_effect
        
        if 'stress_interaction' in model and 'stress_level' in params:
            stress_effect = model['stress_interaction'] * (params['stress_level'] - 5.0)
            total_effect += stress_effect
        
        # Add random noise
        variance = model.get('variance', 0.5)
        noise = np.random.normal(0, variance)
        
        final_value = baseline + total_effect + noise
        
        return max(0.0, final_value)  # Ensure non-negative values


class EnhancementCalculator:
    """Calculates enhancement gains and improvements."""
    
    def __init__(self):
        """Initialize enhancement calculator."""
        self.phenotype_predictor = PhenotypePredictor()
        self.logger = logging.getLogger(__name__)
    
    def calculate_enhancement_gain(self, gene_name: str, baseline_variant: str,
                                 enhanced_variant: str, enhancement_category: EnhancementCategory,
                                 individual_params: Optional[Dict[str, Any]] = None) -> EnhancementGain:
        """
        Calculate enhancement gain from baseline to enhanced variant.
        
        Args:
            gene_name: Gene symbol
            baseline_variant: Baseline genotype
            enhanced_variant: Enhanced genotype
            enhancement_category: Category of enhancement
            individual_params: Individual-specific parameters
            
        Returns:
            EnhancementGain object
        """
        try:
            # Predict phenotypes for both variants
            baseline_phenotypes = self.phenotype_predictor.predict_phenotype(
                gene_name, baseline_variant, individual_params
            )
            enhanced_phenotypes = self.phenotype_predictor.predict_phenotype(
                gene_name, enhanced_variant, individual_params
            )
            
            # Select primary metric based on gene and category
            primary_metric = self._select_primary_metric(gene_name, enhancement_category)
            
            if primary_metric not in baseline_phenotypes:
                raise SimulationError(f"Primary metric {primary_metric} not found for {gene_name}")
            
            baseline_value = baseline_phenotypes[primary_metric]
            enhanced_value = enhanced_phenotypes[primary_metric]
            
            # Calculate improvement factor
            if baseline_value > 0:
                improvement_factor = enhanced_value / baseline_value
            else:
                improvement_factor = 1.0
            
            # Calculate confidence interval
            confidence_interval = self._calculate_confidence_interval(
                baseline_value, enhanced_value, gene_name, primary_metric
            )
            
            # Estimate population percentile
            population_percentile = self._estimate_population_percentile(
                enhanced_value, enhancement_category
            )
            
            return EnhancementGain(
                category=enhancement_category,
                primary_metric=primary_metric,
                baseline_value=baseline_value,
                enhanced_value=enhanced_value,
                improvement_factor=improvement_factor,
                confidence_interval=confidence_interval,
                population_percentile=population_percentile
            )
            
        except Exception as e:
            self.logger.error(f"Failed to calculate enhancement gain: {e}")
            raise SimulationError(f"Enhancement calculation failed: {e}")
    
    def _select_primary_metric(self, gene_name: str, category: EnhancementCategory) -> str:
        """Select primary metric for enhancement calculation."""
        primary_metrics = {
            'COMT': 'cognitive_performance',
            'BDNF': 'learning_rate', 
            'ACTN3': 'sprint_performance',
            'FOXO3': 'longevity_score',
            'MSTN': 'muscle_mass'
        }
        
        return primary_metrics.get(gene_name, 'overall_score')
    
    def _calculate_confidence_interval(self, baseline: float, enhanced: float,
                                     gene_name: str, metric: str) -> Tuple[float, float]:
        """Calculate confidence interval for enhancement effect."""
        # Get variance from model
        model = self.phenotype_predictor.phenotype_models.get(gene_name, {})
        metric_model = model.get(metric, {})
        variance = metric_model.get('variance', 1.0)
        
        # Calculate standard error
        se = variance / np.sqrt(100)  # Assume n=100 for effect size estimation
        
        # 95% confidence interval
        margin_of_error = 1.96 * se
        effect_size = enhanced - baseline
        
        ci_lower = effect_size - margin_of_error
        ci_upper = effect_size + margin_of_error
        
        return (baseline + ci_lower, baseline + ci_upper)
    
    def _estimate_population_percentile(self, value: float, 
                                      category: EnhancementCategory) -> float:
        """Estimate population percentile for a given value."""
        pop_model = self.phenotype_predictor.population_baselines.get(category)
        
        if not pop_model:
            return 50.0  # Default to median
        
        # Calculate z-score
        z_score = (value - pop_model.mean_baseline) / pop_model.std_baseline
        
        # Convert to percentile
        if SCIPY_AVAILABLE:
            percentile = stats.norm.cdf(z_score) * 100
        else:
            # Approximate using error function
            percentile = 50 + 50 * np.tanh(z_score / np.sqrt(2))
        
        return max(0.0, min(100.0, percentile))


class PopulationSimulator:
    """Simulates enhancement effects at population level."""
    
    def __init__(self):
        """Initialize population simulator."""
        self.enhancement_calculator = EnhancementCalculator()
        self.logger = logging.getLogger(__name__)
    
    def simulate_population_effect(self, gene_name: str, enhancement_variant: str,
                                 enhancement_category: EnhancementCategory,
                                 population_size: int = 10000,
                                 baseline_allele_frequency: float = 0.5) -> Dict[str, Any]:
        """
        Simulate enhancement effect across a population.
        
        Args:
            gene_name: Gene symbol
            enhancement_variant: Variant to introduce
            enhancement_category: Enhancement category
            population_size: Number of individuals to simulate
            baseline_allele_frequency: Frequency of baseline allele
            
        Returns:
            Population simulation results
        """
        results = {
            'population_size': population_size,
            'baseline_distribution': [],
            'enhanced_distribution': [],
            'improvement_factors': [],
            'population_statistics': {}
        }
        
        # Generate diverse population
        individuals = self._generate_population_diversity(population_size)
        
        for individual in individuals:
            try:
                # Determine baseline genotype based on allele frequency
                baseline_variant = self._assign_baseline_genotype(baseline_allele_frequency)
                
                # Calculate enhancement gain for this individual
                enhancement_gain = self.enhancement_calculator.calculate_enhancement_gain(
                    gene_name, baseline_variant, enhancement_variant, 
                    enhancement_category, individual
                )
                
                results['baseline_distribution'].append(enhancement_gain.baseline_value)
                results['enhanced_distribution'].append(enhancement_gain.enhanced_value)
                results['improvement_factors'].append(enhancement_gain.improvement_factor)
                
            except Exception as e:
                self.logger.warning(f"Failed to simulate individual: {e}")
                continue
        
        # Calculate population statistics
        results['population_statistics'] = self._calculate_population_statistics(results)
        
        return results
    
    def _generate_population_diversity(self, size: int) -> List[Dict[str, Any]]:
        """Generate diverse population with various characteristics."""
        individuals = []
        
        for _ in range(size):
            individual = {
                'age': np.random.normal(40, 15),  # Age distribution
                'sex': np.random.choice(['M', 'F']),
                'stress_level': np.random.normal(5.0, 2.0),
                'fitness_level': np.random.normal(5.0, 2.0),
                'genetic_background': np.random.normal(0, 1),  # Random genetic background
                'environmental_factors': np.random.normal(0, 1)
            }
            
            # Ensure reasonable bounds
            individual['age'] = max(18, min(80, individual['age']))
            individual['stress_level'] = max(1, min(10, individual['stress_level']))
            individual['fitness_level'] = max(1, min(10, individual['fitness_level']))
            
            individuals.append(individual)
        
        return individuals
    
    def _assign_baseline_genotype(self, allele_frequency: float) -> str:
        """Assign baseline genotype based on allele frequency."""
        # Simple Hardy-Weinberg assignment
        p = allele_frequency  # Frequency of "normal" allele
        q = 1 - p            # Frequency of "variant" allele
        
        rand = np.random.random()
        
        if rand < p * p:
            return "homozygous_normal"
        elif rand < p * p + 2 * p * q:
            return "heterozygous"
        else:
            return "homozygous_variant"
    
    def _calculate_population_statistics(self, results: Dict[str, Any]) -> Dict[str, float]:
        """Calculate statistics for population simulation."""
        baseline_values = np.array(results['baseline_distribution'])
        enhanced_values = np.array(results['enhanced_distribution'])
        improvement_factors = np.array(results['improvement_factors'])
        
        if len(baseline_values) == 0:
            return {}
        
        stats = {
            'baseline_mean': np.mean(baseline_values),
            'baseline_std': np.std(baseline_values),
            'enhanced_mean': np.mean(enhanced_values),
            'enhanced_std': np.std(enhanced_values),
            'mean_improvement_factor': np.mean(improvement_factors),
            'improvement_factor_std': np.std(improvement_factors),
            'percent_improved': np.sum(improvement_factors > 1.0) / len(improvement_factors) * 100,
            'median_improvement': np.median(improvement_factors),
            'max_improvement': np.max(improvement_factors),
            'min_improvement': np.min(improvement_factors)
        }
        
        # Calculate percentiles
        stats['improvement_25th_percentile'] = np.percentile(improvement_factors, 25)
        stats['improvement_75th_percentile'] = np.percentile(improvement_factors, 75)
        
        return stats


class DoseResponseSimulator:
    """Simulates dose-response relationships for enhancement interventions."""
    
    def __init__(self):
        """Initialize dose-response simulator."""
        self.logger = logging.getLogger(__name__)
    
    def simulate_dose_response(self, gene_name: str, enhancement_type: str,
                             dose_range: Tuple[float, float] = (0.0, 1.0),
                             num_points: int = 20) -> Dict[str, Any]:
        """
        Simulate dose-response relationship.
        
        Args:
            gene_name: Gene symbol
            enhancement_type: Type of enhancement intervention
            dose_range: Range of doses to simulate
            num_points: Number of dose points
            
        Returns:
            Dose-response simulation results
        """
        doses = np.linspace(dose_range[0], dose_range[1], num_points)
        
        # Get dose-response model for gene
        model = self._get_dose_response_model(gene_name, enhancement_type)
        
        responses = []
        for dose in doses:
            response = self._calculate_dose_response(dose, model)
            responses.append(response)
        
        responses = np.array(responses)
        
        # Calculate derived metrics
        results = {
            'doses': doses.tolist(),
            'responses': responses.tolist(),
            'model_parameters': model.__dict__,
            'ec50': model.ec50,
            'max_effect': model.max_effect,
            'hill_coefficient': model.hill_coefficient,
            'therapeutic_window': self._calculate_therapeutic_window(doses, responses),
            'optimal_dose': self._find_optimal_dose(doses, responses)
        }
        
        return results
    
    def _get_dose_response_model(self, gene_name: str, enhancement_type: str) -> DoseResponseModel:
        """Get dose-response model parameters for gene and enhancement type."""
        # Default models for different genes and intervention types
        models = {
            'COMT': {
                'base_edit': DoseResponseModel(
                    ec50=0.3, hill_coefficient=2.0, max_effect=1.5, baseline_effect=1.0
                ),
                'knockout': DoseResponseModel(
                    ec50=0.5, hill_coefficient=1.5, max_effect=2.0, baseline_effect=1.0
                )
            },
            'BDNF': {
                'base_edit': DoseResponseModel(
                    ec50=0.4, hill_coefficient=1.8, max_effect=1.3, baseline_effect=1.0
                ),
                'overexpression': DoseResponseModel(
                    ec50=0.2, hill_coefficient=2.5, max_effect=1.8, baseline_effect=1.0
                )
            },
            'ACTN3': {
                'knockout': DoseResponseModel(
                    ec50=0.6, hill_coefficient=1.2, max_effect=1.4, baseline_effect=1.0
                )
            },
            'MSTN': {
                'knockout': DoseResponseModel(
                    ec50=0.4, hill_coefficient=1.8, max_effect=2.5, baseline_effect=1.0
                ),
                'inhibition': DoseResponseModel(
                    ec50=0.3, hill_coefficient=2.2, max_effect=2.0, baseline_effect=1.0
                )
            }
        }
        
        # Get model or use default
        gene_models = models.get(gene_name, {})
        model = gene_models.get(enhancement_type)
        
        if not model:
            # Default model
            model = DoseResponseModel(
                ec50=0.5, hill_coefficient=2.0, max_effect=1.5, baseline_effect=1.0
            )
        
        return model
    
    def _calculate_dose_response(self, dose: float, model: DoseResponseModel) -> float:
        """Calculate response for a given dose using Hill equation."""
        # Hill equation: Response = baseline + (max_effect - baseline) * dose^n / (EC50^n + dose^n)
        numerator = (model.max_effect - model.baseline_effect) * (dose ** model.hill_coefficient)
        denominator = (model.ec50 ** model.hill_coefficient) + (dose ** model.hill_coefficient)
        
        if denominator == 0:
            return model.baseline_effect
        
        response = model.baseline_effect + (numerator / denominator)
        
        # Add some noise
        noise = np.random.normal(0, 0.05)
        response += noise
        
        return max(0, response)
    
    def _calculate_therapeutic_window(self, doses: np.ndarray, responses: np.ndarray) -> Dict[str, float]:
        """Calculate therapeutic window metrics."""
        # Find effective dose range (e.g., 20-80% of max effect)
        max_response = np.max(responses)
        min_effective = 0.2 * max_response
        max_safe = 0.8 * max_response
        
        effective_indices = np.where((responses >= min_effective) & (responses <= max_safe))[0]
        
        if len(effective_indices) == 0:
            return {'width': 0.0, 'min_dose': 0.0, 'max_dose': 0.0}
        
        min_dose = doses[effective_indices[0]]
        max_dose = doses[effective_indices[-1]]
        width = max_dose - min_dose
        
        return {
            'width': width,
            'min_effective_dose': min_dose,
            'max_safe_dose': max_dose
        }
    
    def _find_optimal_dose(self, doses: np.ndarray, responses: np.ndarray) -> float:
        """Find optimal dose (maximum efficacy with minimal side effects)."""
        # For now, use dose that gives 80% of maximum response
        max_response = np.max(responses)
        target_response = 0.8 * max_response
        
        # Find closest dose
        diff = np.abs(responses - target_response)
        optimal_idx = np.argmin(diff)
        
        return doses[optimal_idx]


class EffectSimulator:
    """Main effect simulation class integrating all simulation components."""
    
    def __init__(self):
        """Initialize effect simulator."""
        self.phenotype_predictor = PhenotypePredictor()
        self.enhancement_calculator = EnhancementCalculator()
        self.population_simulator = PopulationSimulator()
        self.dose_response_simulator = DoseResponseSimulator()
        
        self.logger = logging.getLogger(__name__)
    
    def simulate_variant_effect(self, gene_name: str, variant: str,
                              enhancement_category: EnhancementCategory,
                              simulation_params: Optional[Dict[str, Any]] = None) -> VariantEffect:
        """
        Comprehensive simulation of variant effect.
        
        Args:
            gene_name: Gene symbol
            variant: Variant identifier
            enhancement_category: Enhancement category
            simulation_params: Simulation parameters
            
        Returns:
            VariantEffect object with comprehensive predictions
        """
        try:
            # Default simulation parameters
            params = {
                'individual_params': {'age': 30, 'sex': 'M'},
                'baseline_variant': 'wild_type',
                'population_size': 1000,
                'include_side_effects': True
            }
            if simulation_params:
                params.update(simulation_params)
            
            # Create variant info object
            variant_info = VariantInfo(name=variant)
            
            # Calculate enhancement gain
            enhancement_gain = self.enhancement_calculator.calculate_enhancement_gain(
                gene_name, params['baseline_variant'], variant, enhancement_category,
                params['individual_params']
            )
            
            # Predict protein effect (simplified)
            protein_effect = self._predict_protein_effect(gene_name, variant)
            
            # Predict side effects
            side_effects = []
            if params['include_side_effects']:
                side_effects = self._predict_side_effects(gene_name, variant, enhancement_gain)
            
            # Estimate population frequency (if available)
            population_frequency = self._estimate_population_frequency(gene_name, variant)
            
            variant_effect = VariantEffect(
                variant=variant_info,
                protein_effect=protein_effect,
                enhancement_gain=enhancement_gain,
                side_effects=side_effects,
                population_frequency=population_frequency
            )
            
            self.logger.info(f"Simulated effect for {gene_name} {variant}: "
                           f"{enhancement_gain.improvement_factor:.2f}x improvement")
            
            return variant_effect
            
        except Exception as e:
            self.logger.error(f"Failed to simulate variant effect: {e}")
            raise SimulationError(f"Variant effect simulation failed: {e}")
    
    def _predict_protein_effect(self, gene_name: str, variant: str) -> ProteinEffect:
        """Predict protein-level effects of variant."""
        # Simplified protein effect prediction
        # In practice, this would use structural modeling tools
        
        effect_predictions = {
            'COMT': {
                'Val158Met': ProteinEffect(
                    stability_change=-0.5,  # Slightly less stable
                    activity_change=0.7,    # Reduced activity (beneficial for cognition)
                    structure_disrupted=False,
                    confidence_score=0.8
                )
            },
            'BDNF': {
                'Val66Met': ProteinEffect(
                    stability_change=-0.8,
                    activity_change=0.6,    # Reduced activity
                    structure_disrupted=False,
                    confidence_score=0.7
                )
            },
            'ACTN3': {
                'R577X': ProteinEffect(
                    stability_change=-5.0,  # Complete loss of function
                    activity_change=0.0,    # No activity
                    structure_disrupted=True,
                    confidence_score=0.9
                )
            },
            'MSTN': {
                'knockout': ProteinEffect(
                    stability_change=-10.0,  # Complete loss
                    activity_change=0.0,
                    structure_disrupted=True,
                    confidence_score=0.95
                )
            }
        }
        
        gene_effects = effect_predictions.get(gene_name, {})
        effect = gene_effects.get(variant)
        
        if not effect:
            # Default prediction
            effect = ProteinEffect(
                stability_change=np.random.normal(-1.0, 0.5),
                activity_change=np.random.uniform(0.5, 1.5),
                structure_disrupted=np.random.random() < 0.1,
                confidence_score=0.5
            )
        
        return effect
    
    def _predict_side_effects(self, gene_name: str, variant: str, 
                            enhancement_gain: EnhancementGain) -> List[SideEffect]:
        """Predict potential side effects of enhancement."""
        side_effects = []
        
        # Gene-specific side effect profiles
        side_effect_profiles = {
            'COMT': [
                SideEffect(
                    description="Reduced performance in low-stress environments",
                    probability=0.6,
                    severity="mild",
                    reversible=False,
                    onset_time="immediate",
                    evidence_level="human_case"
                ),
                SideEffect(
                    description="Altered response to dopaminergic medications",
                    probability=0.3,
                    severity="moderate",
                    reversible=True,
                    onset_time="immediate",
                    evidence_level="theoretical"
                )
            ],
            'BDNF': [
                SideEffect(
                    description="Potential mood regulation changes",
                    probability=0.4,
                    severity="mild",
                    reversible=True,
                    onset_time="weeks",
                    evidence_level="animal_studies"
                )
            ],
            'ACTN3': [
                SideEffect(
                    description="Reduced endurance capacity",
                    probability=0.7,
                    severity="mild",
                    reversible=False,
                    onset_time="immediate",
                    evidence_level="human_case"
                ),
                SideEffect(
                    description="Increased injury risk during explosive movements",
                    probability=0.2,
                    severity="moderate",
                    reversible=True,
                    onset_time="immediate",
                    evidence_level="theoretical"
                )
            ],
            'MSTN': [
                SideEffect(
                    description="Excessive muscle growth (muscle hypertrophy)",
                    probability=0.9,
                    severity="moderate",
                    reversible=False,
                    onset_time="weeks",
                    evidence_level="animal_studies"
                ),
                SideEffect(
                    description="Potential cardiac muscle effects",
                    probability=0.3,
                    severity="severe",
                    reversible=False,
                    onset_time="months",
                    evidence_level="theoretical"
                )
            ],
            'FOXO3': [
                SideEffect(
                    description="Altered stress response patterns",
                    probability=0.4,
                    severity="mild",
                    reversible=True,
                    onset_time="months",
                    evidence_level="theoretical"
                )
            ]
        }
        
        gene_side_effects = side_effect_profiles.get(gene_name, [])
        
        # Adjust probabilities based on enhancement strength
        improvement_factor = enhancement_gain.improvement_factor
        probability_multiplier = min(2.0, improvement_factor)  # Stronger effects = more side effects
        
        for side_effect in gene_side_effects:
            adjusted_side_effect = SideEffect(
                description=side_effect.description,
                probability=min(1.0, side_effect.probability * probability_multiplier),
                severity=side_effect.severity,
                reversible=side_effect.reversible,
                onset_time=side_effect.onset_time,
                evidence_level=side_effect.evidence_level
            )
            side_effects.append(adjusted_side_effect)
        
        return side_effects
    
    def _estimate_population_frequency(self, gene_name: str, variant: str) -> Optional[float]:
        """Estimate population frequency of variant."""
        # Known population frequencies for common variants
        frequencies = {
            'COMT': {
                'Val158Met': 0.6,  # Met allele frequency in Europeans
                'Met158Met': 0.25  # Met/Met genotype frequency
            },
            'BDNF': {
                'Val66Met': 0.3,   # Met allele frequency
                'Met66Met': 0.04   # Met/Met genotype frequency
            },
            'ACTN3': {
                'R577X': 0.45,     # X allele frequency
                'X577X': 0.18      # X/X genotype frequency
            },
            'FOXO3': {
                'rs2802292_T': 0.4  # Protective allele frequency
            }
        }
        
        gene_frequencies = frequencies.get(gene_name, {})
        return gene_frequencies.get(variant)
    
    def run_comprehensive_simulation(self, gene_name: str, variants: List[str],
                                   enhancement_category: EnhancementCategory,
                                   simulation_config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Run comprehensive simulation for multiple variants.
        
        Args:
            gene_name: Gene symbol
            variants: List of variants to simulate
            enhancement_category: Enhancement category
            simulation_config: Configuration parameters
            
        Returns:
            Comprehensive simulation results
        """
        config = {
            'population_size': 10000,
            'include_dose_response': True,
            'include_population_simulation': True,
            'monte_carlo_iterations': 1000,
            'confidence_level': 0.95
        }
        if simulation_config:
            config.update(simulation_config)
        
        results = {
            'gene_name': gene_name,
            'enhancement_category': enhancement_category.value,
            'simulation_config': config,
            'variant_effects': {},
            'population_simulations': {},
            'dose_response_curves': {},
            'comparative_analysis': {}
        }
        
        # Simulate each variant
        for variant in variants:
            try:
                # Individual variant effect
                variant_effect = self.simulate_variant_effect(
                    gene_name, variant, enhancement_category
                )
                results['variant_effects'][variant] = variant_effect
                
                # Population simulation
                if config['include_population_simulation']:
                    pop_sim = self.population_simulator.simulate_population_effect(
                        gene_name, variant, enhancement_category, config['population_size']
                    )
                    results['population_simulations'][variant] = pop_sim
                
                # Dose-response simulation
                if config['include_dose_response']:
                    dose_response = self.dose_response_simulator.simulate_dose_response(
                        gene_name, 'base_edit'  # Default intervention type
                    )
                    results['dose_response_curves'][variant] = dose_response
                
            except Exception as e:
                self.logger.error(f"Failed to simulate variant {variant}: {e}")
                continue
        
        # Comparative analysis
        results['comparative_analysis'] = self._perform_comparative_analysis(
            results['variant_effects']
        )
        
        return results
    
    def _perform_comparative_analysis(self, variant_effects: Dict[str, VariantEffect]) -> Dict[str, Any]:
        """Perform comparative analysis of variant effects."""
        if not variant_effects:
            return {}
        
        # Extract improvement factors
        improvements = {
            variant: effect.enhancement_gain.improvement_factor
            for variant, effect in variant_effects.items()
        }
        
        # Find best and worst variants
        best_variant = max(improvements.keys(), key=lambda x: improvements[x])
        worst_variant = min(improvements.keys(), key=lambda x: improvements[x])
        
        # Calculate statistics
        improvement_values = list(improvements.values())
        mean_improvement = np.mean(improvement_values)
        std_improvement = np.std(improvement_values)
        
        # Risk-benefit analysis
        risk_benefit_scores = {}
        for variant, effect in variant_effects.items():
            benefit = effect.enhancement_gain.improvement_factor
            risk = sum(se.risk_score for se in effect.side_effects)
            risk_benefit_scores[variant] = benefit / (1 + risk) if risk > 0 else benefit
        
        best_risk_benefit = max(risk_benefit_scores.keys(), 
                               key=lambda x: risk_benefit_scores[x])
        
        return {
            'best_improvement_variant': best_variant,
            'best_improvement_factor': improvements[best_variant],
            'worst_improvement_variant': worst_variant,
            'worst_improvement_factor': improvements[worst_variant],
            'mean_improvement_factor': mean_improvement,
            'std_improvement_factor': std_improvement,
            'best_risk_benefit_variant': best_risk_benefit,
            'risk_benefit_scores': risk_benefit_scores,
            'total_variants_analyzed': len(variant_effects)
        }
    
    def predict_long_term_outcomes(self, gene_name: str, variant: str,
                                 enhancement_category: EnhancementCategory,
                                 time_horizon_years: int = 10) -> Dict[str, Any]:
        """
        Predict long-term outcomes of enhancement.
        
        Args:
            gene_name: Gene symbol
            variant: Variant identifier
            enhancement_category: Enhancement category
            time_horizon_years: Time horizon for prediction
            
        Returns:
            Long-term outcome predictions
        """
        # Simulate variant effect
        variant_effect = self.simulate_variant_effect(gene_name, variant, enhancement_category)
        
        # Model time-dependent changes
        time_points = np.linspace(0, time_horizon_years, 20)
        
        outcomes = {
            'time_points': time_points.tolist(),
            'enhancement_trajectory': [],
            'side_effect_trajectory': [],
            'cumulative_benefit': [],
            'cumulative_risk': []
        }
        
        initial_improvement = variant_effect.enhancement_gain.improvement_factor
        cumulative_benefit = 0
        cumulative_risk = 0
        
        for t in time_points:
            # Model enhancement decay or stability
            if enhancement_category == EnhancementCategory.COGNITIVE:
                # Cognitive enhancements may have some decay
                decay_rate = 0.02  # 2% per year
                current_improvement = initial_improvement * np.exp(-decay_rate * t)
            elif enhancement_category == EnhancementCategory.PHYSICAL:
                # Physical enhancements may decay faster
                decay_rate = 0.05  # 5% per year
                current_improvement = initial_improvement * np.exp(-decay_rate * t)
            else:
                # Longevity enhancements may be stable
                current_improvement = initial_improvement
            
            outcomes['enhancement_trajectory'].append(current_improvement)
            
            # Model side effect accumulation
            side_effect_intensity = 0
            for side_effect in variant_effect.side_effects:
                if side_effect.onset_time in ['immediate', 'days']:
                    intensity = side_effect.probability * side_effect.risk_score
                elif side_effect.onset_time == 'weeks':
                    intensity = side_effect.probability * side_effect.risk_score * min(1.0, t / 0.25)
                elif side_effect.onset_time == 'months':
                    intensity = side_effect.probability * side_effect.risk_score * min(1.0, t / 1.0)
                else:  # years
                    intensity = side_effect.probability * side_effect.risk_score * min(1.0, t / 2.0)
                
                side_effect_intensity += intensity
            
            outcomes['side_effect_trajectory'].append(side_effect_intensity)
            
            # Calculate cumulative metrics
            if t > 0:
                dt = time_points[1] - time_points[0]  # Time step
                cumulative_benefit += (current_improvement - 1.0) * dt
                cumulative_risk += side_effect_intensity * dt
            
            outcomes['cumulative_benefit'].append(cumulative_benefit)
            outcomes['cumulative_risk'].append(cumulative_risk)
        
        # Calculate summary metrics
        outcomes['summary'] = {
            'final_enhancement_level': outcomes['enhancement_trajectory'][-1],
            'total_cumulative_benefit': cumulative_benefit,
            'total_cumulative_risk': cumulative_risk,
            'benefit_risk_ratio': cumulative_benefit / max(0.1, cumulative_risk),
            'stable_enhancement': outcomes['enhancement_trajectory'][-1] > 1.1,
            'acceptable_risk_profile': cumulative_risk < 5.0
        }
        
        return outcomes
    
    def optimize_enhancement_strategy(self, gene_targets: List[str],
                                    enhancement_categories: List[EnhancementCategory],
                                    constraints: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Optimize enhancement strategy across multiple genes.
        
        Args:
            gene_targets: List of genes to consider
            enhancement_categories: Enhancement categories
            constraints: Optimization constraints
            
        Returns:
            Optimized enhancement strategy
        """
        default_constraints = {
            'max_total_risk': 10.0,
            'min_improvement_factor': 1.2,
            'max_genes_modified': 3,
            'safety_threshold': 70.0
        }
        if constraints:
            default_constraints.update(constraints)
        
        # Generate all possible combinations
        strategies = []
        
        for num_genes in range(1, min(len(gene_targets), default_constraints['max_genes_modified']) + 1):
            from itertools import combinations
            
            for gene_combo in combinations(gene_targets, num_genes):
                for cat_combo in combinations(enhancement_categories, num_genes):
                    strategy = {
                        'genes': list(gene_combo),
                        'categories': list(cat_combo),
                        'total_improvement': 1.0,
                        'total_risk': 0.0,
                        'variant_effects': {}
                    }
                    
                    # Simulate each gene in the strategy
                    valid_strategy = True
                    for gene, category in zip(gene_combo, cat_combo):
                        try:
                            # Use best known variant for each gene
                            best_variants = {
                                'COMT': 'Met158Met',
                                'BDNF': 'Val66Val',  # Wild-type is better for BDNF
                                'ACTN3': 'R577R',   # Wild-type for power
                                'FOXO3': 'rs2802292_TT',
                                'MSTN': 'knockout'
                            }
                            
                            variant = best_variants.get(gene, 'enhancement_variant')
                            
                            effect = self.simulate_variant_effect(gene, variant, category)
                            strategy['variant_effects'][gene] = effect
                            
                            # Calculate combined effects (simplified multiplicative model)
                            strategy['total_improvement'] *= effect.enhancement_gain.improvement_factor
                            strategy['total_risk'] += sum(se.risk_score for se in effect.side_effects)
                            
                        except Exception as e:
                            self.logger.warning(f"Failed to simulate {gene}: {e}")
                            valid_strategy = False
                            break
                    
                    # Check constraints
                    if (valid_strategy and 
                        strategy['total_improvement'] >= default_constraints['min_improvement_factor'] and
                        strategy['total_risk'] <= default_constraints['max_total_risk']):
                        
                        strategy['benefit_risk_ratio'] = (
                            strategy['total_improvement'] / max(0.1, strategy['total_risk'])
                        )
                        strategies.append(strategy)
        
        # Sort strategies by benefit-risk ratio
        strategies.sort(key=lambda x: x['benefit_risk_ratio'], reverse=True)
        
        return {
            'optimal_strategy': strategies[0] if strategies else None,
            'alternative_strategies': strategies[1:6],  # Top 5 alternatives
            'total_strategies_evaluated': len(strategies),
            'constraints_used': default_constraints
        }