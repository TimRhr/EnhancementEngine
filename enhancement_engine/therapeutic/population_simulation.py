"""
Population simulation module for therapeutic CRISPR interventions.

This module provides comprehensive population-level modeling including:
- Epidemiological impact assessment
- Health economic modeling
- Population genetics effects
- Implementation strategy optimization
- Long-term outcome prediction
- Public health impact evaluation
"""

import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass
from datetime import datetime, timedelta
import itertools
from collections import defaultdict

from ..models.therapeutic_data_classes import (
    PopulationParameters, EpidemiologicalModel, HealthEconomicModel,
    ImplementationStrategy, PopulationOutcome, TherapeuticImpact
)
from ..models.disease_constants import DISEASE_EPIDEMIOLOGY, POPULATION_GENETICS


@dataclass
class PopulationCohort:
    """Population cohort for simulation."""
    cohort_id: str
    size: int
    age_distribution: Dict[str, float]
    sex_distribution: Dict[str, float]
    genetic_risk_distribution: Dict[str, float]
    disease_prevalence: float
    baseline_characteristics: Dict[str, Any]


@dataclass
class InterventionScenario:
    """Therapeutic intervention scenario."""
    scenario_name: str
    target_population: str
    intervention_strategy: str
    coverage_rate: float  # Proportion of eligible population receiving intervention
    implementation_timeline: Dict[str, int]  # Years
    cost_per_treatment: float
    efficacy_parameters: Dict[str, float]


class TherapeuticPopulationSimulator:
    """Main simulator for population-level therapeutic interventions."""
    
    def __init__(self):
        """Initialize population simulator."""
        self.logger = logging.getLogger(__name__)
        self._load_population_data()
        self._load_disease_models()
        self._initialize_economic_models()
    
    def _load_population_data(self) -> None:
        """Load population demographic and genetic data."""
        self.population_data = {
            'demographics': {
                'total_population': 330000000,  # US population
                'age_distribution': {
                    '0-17': 0.22, '18-34': 0.23, '35-54': 0.25, 
                    '55-74': 0.20, '75+': 0.10
                },
                'sex_distribution': {'male': 0.49, 'female': 0.51},
                'ethnicity_distribution': {
                    'caucasian': 0.60, 'hispanic': 0.18, 'african_american': 0.13,
                    'asian': 0.06, 'other': 0.03
                }
            },
            'genetic_data': {
                'rheumatoid_arthritis': {
                    'overall_prevalence': 0.01,
                    'genetic_risk_distribution': {
                        'low': 0.65, 'moderate': 0.25, 'high': 0.08, 'very_high': 0.02
                    },
                    'heritability': 0.60,
                    'environmental_factors': {
                        'smoking': {'prevalence': 0.15, 'risk_increase': 2.0},
                        'infections': {'prevalence': 0.30, 'risk_increase': 1.5}
                    }
                }
            }
        }
    
    def _load_disease_models(self) -> None:
        """Load disease progression and epidemiological models."""
        self.disease_models = {
            'rheumatoid_arthritis': {
                'incidence_rate_per_100k': 50,  # New cases per 100,000 per year
                'age_incidence_curve': {
                    '20-29': 0.1, '30-39': 0.15, '40-49': 0.25,
                    '50-59': 0.30, '60-69': 0.20
                },
                'progression_model': {
                    'remission_rate': 0.20,
                    'mild_to_moderate': 0.40,
                    'moderate_to_severe': 0.30,
                    'mortality_increase': 1.3
                },
                'economic_burden': {
                    'annual_direct_cost': 15000,  # USD per patient
                    'annual_indirect_cost': 12000,  # Lost productivity
                    'qaly_loss_per_year': 0.3
                }
            }
        }
    
    def _initialize_economic_models(self) -> None:
        """Initialize health economic models."""
        self.economic_models = {
            'cost_effectiveness_thresholds': {
                'us': 100000,  # $100k per QALY
                'uk': 50000,   # £50k per QALY
                'eu': 75000    # €75k per QALY
            },
            'discount_rates': {
                'costs': 0.03,    # 3% annual discount
                'benefits': 0.03  # 3% annual discount
            },
            'time_horizon_years': 20
        }
    
    def simulate_population_intervention(self, intervention_data: Dict[str, Any],
                                       scenario_parameters: Dict[str, Any]) -> PopulationOutcome:
        """
        Simulate population-level therapeutic intervention.
        
        Args:
            intervention_data: Therapeutic intervention details
            scenario_parameters: Simulation scenario parameters
            
        Returns:
            Population-level outcomes and impacts
        """
        # Create population cohorts
        cohorts = self._create_population_cohorts(scenario_parameters)
        
        # Create intervention scenarios
        scenarios = self._create_intervention_scenarios(intervention_data, scenario_parameters)
        
        # Run simulation for each scenario
        scenario_results = {}
        
        for scenario in scenarios:
            results = self._simulate_scenario(cohorts, scenario, scenario_parameters)
            scenario_results[scenario.scenario_name] = results
        
        # Compare scenarios
        comparison = self._compare_scenarios(scenario_results)
        
        # Calculate population-level metrics
        population_metrics = self._calculate_population_metrics(scenario_results, cohorts)
        
        # Economic analysis
        economic_analysis = self._perform_economic_analysis(scenario_results, scenario_parameters)
        
        return PopulationOutcome(
            intervention_name=intervention_data.get('name', 'CRISPR Therapeutic'),
            simulation_date=datetime.now(),
            population_size=sum(cohort.size for cohort in cohorts),
            scenario_results=scenario_results,
            scenario_comparison=comparison,
            population_metrics=population_metrics,
            economic_analysis=economic_analysis,
            implementation_recommendations=self._generate_implementation_recommendations(
                comparison, economic_analysis
            )
        )
    
    def _create_population_cohorts(self, scenario_parameters: Dict[str, Any]) -> List[PopulationCohort]:
        """Create population cohorts for simulation."""
        target_disease = scenario_parameters.get('target_disease', 'rheumatoid_arthritis')
        population_size = scenario_parameters.get('population_size', 100000)
        
        disease_data = self.disease_models.get(target_disease, {})
        genetic_data = self.population_data['genetic_data'].get(target_disease, {})
        
        cohorts = []
        
        # Create cohorts by genetic risk level
        risk_distribution = genetic_data.get('genetic_risk_distribution', {'low': 1.0})
        
        for risk_level, proportion in risk_distribution.items():
            cohort_size = int(population_size * proportion)
            
            # Calculate disease prevalence for this risk group
            base_prevalence = genetic_data.get('overall_prevalence', 0.01)
            risk_multipliers = {'low': 0.5, 'moderate': 1.0, 'high': 3.0, 'very_high': 8.0}
            risk_prevalence = base_prevalence * risk_multipliers.get(risk_level, 1.0)
            
            cohort = PopulationCohort(
                cohort_id=f"{target_disease}_{risk_level}_risk",
                size=cohort_size,
                age_distribution=self.population_data['demographics']['age_distribution'],
                sex_distribution=self.population_data['demographics']['sex_distribution'],
                genetic_risk_distribution={risk_level: 1.0},
                disease_prevalence=risk_prevalence,
                baseline_characteristics={
                    'risk_level': risk_level,
                    'intervention_eligible': risk_level in ['high', 'very_high'],
                    'baseline_qaly': 0.8 if risk_level == 'low' else 0.7
                }
            )
            cohorts.append(cohort)
        
        return cohorts
    
    def _create_intervention_scenarios(self, intervention_data: Dict[str, Any],
                                     scenario_parameters: Dict[str, Any]) -> List[InterventionScenario]:
        """Create intervention scenarios for comparison."""
        scenarios = []
        
        # Baseline scenario (no intervention)
        baseline = InterventionScenario(
            scenario_name='baseline_no_intervention',
            target_population='none',
            intervention_strategy='standard_care',
            coverage_rate=0.0,
            implementation_timeline={'start_year': 0, 'ramp_up_years': 0},
            cost_per_treatment=0.0,
            efficacy_parameters={'disease_modification': 0.0, 'safety_profile': 1.0}
        )
        scenarios.append(baseline)
        
        # High-risk only scenario
        high_risk_scenario = InterventionScenario(
            scenario_name='high_risk_targeted',
            target_population='high_and_very_high_risk',
            intervention_strategy=intervention_data.get('strategy', 'base_editing'),
            coverage_rate=0.7,  # 70% coverage
            implementation_timeline={'start_year': 1, 'ramp_up_years': 3},
            cost_per_treatment=intervention_data.get('cost_per_treatment', 200000),
            efficacy_parameters={
                'disease_modification': intervention_data.get('efficacy', 0.6),
                'safety_profile': intervention_data.get('safety_score', 0.8) / 100
            }
        )
        scenarios.append(high_risk_scenario)
        
        # Broader population scenario
        broad_scenario = InterventionScenario(
            scenario_name='moderate_risk_and_above',
            target_population='moderate_high_very_high_risk',
            intervention_strategy=intervention_data.get('strategy', 'base_editing'),
            coverage_rate=0.5,  # 50% coverage (lower due to broader population)
            implementation_timeline={'start_year': 2, 'ramp_up_years': 5},
            cost_per_treatment=intervention_data.get('cost_per_treatment', 200000),
            efficacy_parameters={
                'disease_modification': intervention_data.get('efficacy', 0.6) * 0.8,  # Reduced efficacy in moderate risk
                'safety_profile': intervention_data.get('safety_score', 0.8) / 100
            }
        )
        scenarios.append(broad_scenario)
        
        return scenarios
    
    def _simulate_scenario(self, cohorts: List[PopulationCohort],
                         scenario: InterventionScenario,
                         parameters: Dict[str, Any]) -> Dict[str, Any]:
        """Simulate a specific intervention scenario."""
        simulation_years = parameters.get('simulation_years', 20)
        results = {
            'scenario': scenario,
            'annual_outcomes': [],
            'cumulative_metrics': {},
            'treated_population': 0,
            'total_costs': 0.0,
            'total_qalys_gained': 0.0
        }
        
        # Simulate each year
        for year in range(simulation_years):
            year_results = self._simulate_year(cohorts, scenario, year, parameters)
            results['annual_outcomes'].append(year_results)
            
            # Update cumulative metrics
            results['treated_population'] += year_results['new_treatments']
            results['total_costs'] += year_results['total_costs']
            results['total_qalys_gained'] += year_results['qalys_gained']
        
        # Calculate summary metrics
        results['cumulative_metrics'] = self._calculate_cumulative_metrics(results)
        
        return results
    
    def _simulate_year(self, cohorts: List[PopulationCohort],
                      scenario: InterventionScenario, year: int,
                      parameters: Dict[str, Any]) -> Dict[str, Any]:
        """Simulate outcomes for a single year."""
        target_disease = parameters.get('target_disease', 'rheumatoid_arthritis')
        disease_model = self.disease_models[target_disease]
        
        # Calculate intervention coverage for this year
        coverage = self._calculate_year_coverage(scenario, year)
        
        # Initialize year results
        year_results = {
            'year': year,
            'new_cases': 0,
            'new_treatments': 0,
            'prevented_cases': 0,
            'total_costs': 0.0,
            'qalys_gained': 0.0,
            'adverse_events': 0
        }
        
        # Simulate each cohort
        for cohort in cohorts:
            cohort_results = self._simulate_cohort_year(
                cohort, scenario, coverage, disease_model, year
            )
            
            # Aggregate cohort results
            for metric, value in cohort_results.items():
                if metric in year_results:
                    year_results[metric] += value
        
        return year_results
    
    def _simulate_cohort_year(self, cohort: PopulationCohort,
                            scenario: InterventionScenario, coverage: float,
                            disease_model: Dict[str, Any], year: int) -> Dict[str, Any]:
        """Simulate outcomes for a cohort in a single year."""
        # Calculate new disease cases
        incidence_rate = disease_model['incidence_rate_per_100k'] / 100000
        expected_new_cases = cohort.size * incidence_rate * cohort.disease_prevalence
        
        # Determine if cohort is eligible for intervention
        eligible = self._is_cohort_eligible(cohort, scenario.target_population)
        
        # Calculate interventions
        new_treatments = 0
        prevented_cases = 0
        
        if eligible and coverage > 0:
            # Number of people receiving intervention
            new_treatments = int(expected_new_cases * coverage)
            
            # Calculate prevention efficacy
            efficacy = scenario.efficacy_parameters.get('disease_modification', 0.0)
            prevented_cases = new_treatments * efficacy
        
        # Calculate costs
        treatment_costs = new_treatments * scenario.cost_per_treatment
        
        # Calculate QALYs gained
        qaly_gain_per_prevention = self._calculate_qaly_gain_per_prevention(disease_model)
        qalys_gained = prevented_cases * qaly_gain_per_prevention
        
        # Calculate adverse events
        safety_profile = scenario.efficacy_parameters.get('safety_profile', 1.0)
        adverse_events = new_treatments * (1.0 - safety_profile)
        
        return {
            'new_cases': max(0, expected_new_cases - prevented_cases),
            'new_treatments': new_treatments,
            'prevented_cases': prevented_cases,
            'total_costs': treatment_costs,
            'qalys_gained': qalys_gained,
            'adverse_events': adverse_events
        }
    
    def _calculate_year_coverage(self, scenario: InterventionScenario, year: int) -> float:
        """Calculate intervention coverage for a specific year."""
        start_year = scenario.implementation_timeline.get('start_year', 1)
        ramp_up_years = scenario.implementation_timeline.get('ramp_up_years', 3)
        max_coverage = scenario.coverage_rate
        
        if year < start_year:
            return 0.0
        elif year < start_year + ramp_up_years:
            # Linear ramp-up
            years_since_start = year - start_year + 1
            return max_coverage * (years_since_start / ramp_up_years)
        else:
            return max_coverage
    
    def _is_cohort_eligible(self, cohort: PopulationCohort, target_population: str) -> bool:
        """Determine if cohort is eligible for intervention."""
        risk_level = list(cohort.genetic_risk_distribution.keys())[0]
        
        eligibility_criteria = {
            'none': [],
            'very_high_risk': ['very_high'],
            'high_and_very_high_risk': ['high', 'very_high'],
            'moderate_high_very_high_risk': ['moderate', 'high', 'very_high'],
            'all_risk_levels': ['low', 'moderate', 'high', 'very_high']
        }
        
        eligible_levels = eligibility_criteria.get(target_population, [])
        return risk_level in eligible_levels
    
    def _calculate_qaly_gain_per_prevention(self, disease_model: Dict[str, Any]) -> float:
        """Calculate QALY gain per prevented disease case."""
        # Annual QALY loss from disease
        annual_qaly_loss = disease_model['economic_burden'].get('qaly_loss_per_year', 0.3)
        
        # Average disease duration (assume 20 years with treatment)
        disease_duration_years = 20
        
        # Discount future QALYs
        discount_rate = self.economic_models['discount_rates']['benefits']
        
        # Calculate present value of QALY gains
        total_qaly_gain = 0.0
        for year in range(disease_duration_years):
            discounted_gain = annual_qaly_loss / ((1 + discount_rate) ** year)
            total_qaly_gain += discounted_gain
        
        return total_qaly_gain
    
    def _perform_economic_analysis(self, scenario_results: Dict[str, Any],
                                 parameters: Dict[str, Any]) -> Dict[str, Any]:
        """Perform comprehensive economic analysis."""
        economic_analysis = {}
        
        # Cost-effectiveness analysis for each scenario
        for scenario_name, results in scenario_results.items():
            if scenario_name == 'baseline_no_intervention':
                continue  # Skip baseline for cost-effectiveness
            
            # Compare to baseline
            baseline = scenario_results['baseline_no_intervention']
            
            # Incremental costs and QALYs
            incremental_costs = results['total_costs'] - baseline['total_costs']
            incremental_qalys = results['total_qalys_gained'] - baseline['total_qalys_gained']
            
            # Cost-effectiveness ratio
            if incremental_qalys > 0:
                icer = incremental_costs / incremental_qalys
            else:
                icer = float('inf')
            
            # Budget impact
            total_population = parameters.get('population_size', 100000)
            annual_budget_impact = incremental_costs / parameters.get('simulation_years', 20)
            
            economic_analysis[scenario_name] = {
                'incremental_costs': incremental_costs,
                'incremental_qalys': incremental_qalys,
                'icer': icer,
                'cost_effective': icer <= self.economic_models['cost_effectiveness_thresholds']['us'],
                'annual_budget_impact': annual_budget_impact,
                'treated_population': results['treated_population'],
                'cost_per_case_prevented': (
                    incremental_costs / sum(year['prevented_cases'] for year in results['annual_outcomes'])
                    if sum(year['prevented_cases'] for year in results['annual_outcomes']) > 0 else float('inf')
                )
            }
        
        return economic_analysis
    
    def _compare_scenarios(self, scenario_results: Dict[str, Any]) -> Dict[str, Any]:
        """Compare different intervention scenarios."""
        comparison = {
            'best_scenario': None,
            'ranking': [],
            'trade_offs': {},
            'sensitivity_analysis': {}
        }
        
        # Rank scenarios by cost-effectiveness
        scenario_scores = {}
        
        for scenario_name, results in scenario_results.items():
            if scenario_name == 'baseline_no_intervention':
                continue
            
            # Simple scoring: combine efficacy and cost-effectiveness
            qalys_gained = results['total_qalys_gained']
            total_costs = results['total_costs']
            
            # Normalize scores (higher is better)
            efficacy_score = min(100, qalys_gained * 10)  # Cap at 100
            cost_score = max(0, 100 - (total_costs / 1000000))  # $1M = 0 points
            
            combined_score = efficacy_score * 0.6 + cost_score * 0.4
            scenario_scores[scenario_name] = combined_score
        
        # Rank scenarios
        comparison['ranking'] = sorted(scenario_scores.items(), 
                                     key=lambda x: x[1], reverse=True)
        
        if comparison['ranking']:
            comparison['best_scenario'] = comparison['ranking'][0][0]
        
        return comparison
    
    def model_genetic_drift(self, intervention_data: Dict[str, Any],
                          generations: int = 10) -> Dict[str, Any]:
        """
        Model long-term genetic effects of population intervention.
        
        Args:
            intervention_data: Therapeutic intervention details
            generations: Number of generations to simulate
            
        Returns:
            Genetic drift analysis results
        """
        # Initial allele frequencies
        target_gene = intervention_data.get('target_gene', 'PTPN22')
        
        # Simplified allele frequencies (risk allele frequency)
        initial_frequencies = {
            'PTPN22': 0.15,  # 15% carry risk allele
            'HLA-DRB1': 0.25,  # 25% carry shared epitope
            'STAT4': 0.20   # 20% carry risk allele
        }
        
        base_frequency = initial_frequencies.get(target_gene, 0.15)
        
        # Simulate genetic drift
        frequencies_over_time = [base_frequency]
        
        # Parameters
        population_size = 1000000  # Effective population size
        selection_coefficient = 0.0  # No natural selection assumed
        
        # Intervention effects
        intervention_coverage = intervention_data.get('population_coverage', 0.1)
        intervention_efficacy = intervention_data.get('efficacy', 0.6)
        
        # Simulate each generation
        for generation in range(generations):
            current_freq = frequencies_over_time[-1]
            
            # Apply intervention effect (reduces fitness disadvantage)
            if generation >= 1:  # Intervention starts in generation 1
                # Intervention reduces the selection against risk allele
                effective_selection = selection_coefficient * (1 - intervention_coverage * intervention_efficacy)
            else:
                effective_selection = selection_coefficient
            
            # Genetic drift (random sampling)
            drift_variance = current_freq * (1 - current_freq) / (2 * population_size)
            drift_change = np.random.normal(0, np.sqrt(drift_variance))
            
            # Selection effect
            selection_change = effective_selection * current_freq * (1 - current_freq)
            
            # New frequency
            new_freq = current_freq + drift_change + selection_change
            new_freq = max(0.001, min(0.999, new_freq))  # Keep within bounds
            
            frequencies_over_time.append(new_freq)
        
        return {
            'target_gene': target_gene,
            'initial_frequency': base_frequency,
            'final_frequency': frequencies_over_time[-1],
            'frequency_trajectory': frequencies_over_time,
            'relative_change': (frequencies_over_time[-1] - base_frequency) / base_frequency,
            'intervention_impact': self._assess_genetic_intervention_impact(
                base_frequency, frequencies_over_time[-1], intervention_data
            )
        }
    
    def _assess_genetic_intervention_impact(self, initial_freq: float, final_freq: float,
                                          intervention_data: Dict[str, Any]) -> Dict[str, Any]:
        """Assess the impact of intervention on genetic composition."""
        frequency_change = final_freq - initial_freq
        relative_change = frequency_change / initial_freq
        
        # Assess significance
        if abs(relative_change) < 0.01:
            impact_level = 'negligible'
        elif abs(relative_change) < 0.05:
            impact_level = 'minimal'
        elif abs(relative_change) < 0.10:
            impact_level = 'moderate'
        else:
            impact_level = 'significant'
        
        return {
            'frequency_change': frequency_change,
            'relative_change': relative_change,
            'impact_level': impact_level,
            'population_health_effect': self._estimate_population_health_effect(
                frequency_change, intervention_data
            ),
            'ethical_considerations': self._identify_genetic_ethical_considerations(impact_level)
        }
    
    def generate_population_report(self, population_outcome: PopulationOutcome) -> str:
        """Generate comprehensive population simulation report."""
        report = f"""
POPULATION INTERVENTION SIMULATION REPORT
========================================

INTERVENTION: {population_outcome.intervention_name}
SIMULATION DATE: {population_outcome.simulation_date.strftime('%Y-%m-%d')}
POPULATION SIZE: {population_outcome.population_size:,}

SCENARIO COMPARISON:
Best Scenario: {population_outcome.scenario_comparison.get('best_scenario', 'N/A')}

SCENARIO RANKINGS:
{chr(10).join([f"{i+1}. {scenario}: {score:.1f} points" 
               for i, (scenario, score) in enumerate(population_outcome.scenario_comparison.get('ranking', []))])}

ECONOMIC ANALYSIS:
"""
        
        # Add economic results for each scenario
        for scenario_name, economic_data in population_outcome.economic_analysis.items():
            icer = economic_data['icer']
            if icer == float('inf'):
                icer_str = "Dominated (no QALY gain)"
            else:
                icer_str = f"${icer:,.0f}/QALY"
            
            report += f"""
{scenario_name.upper()}:
- ICER: {icer_str}
- Cost-Effective: {'Yes' if economic_data['cost_effective'] else 'No'}
- Cases Prevented: {economic_data.get('cases_prevented', 'N/A')}
- Annual Budget Impact: ${economic_data['annual_budget_impact']:,.0f}
"""
        
        report += f"""
IMPLEMENTATION RECOMMENDATIONS:
{chr(10).join([f"- {rec}" for rec in population_outcome.implementation_recommendations])}

POPULATION HEALTH IMPACT:
{self._summarize_population_health_impact(population_outcome)}
        """
        
        return report.strip()
    
    # Additional helper methods would be implemented here...
    def _calculate_cumulative_metrics(self, results: Dict[str, Any]) -> Dict[str, float]:
        """Calculate cumulative metrics from annual results."""
        annual_outcomes = results['annual_outcomes']
        
        return {
            'total_cases_prevented': sum(year['prevented_cases'] for year in annual_outcomes),
            'total_treatments_given': sum(year['new_treatments'] for year in annual_outcomes),
            'average_annual_cost': results['total_costs'] / len(annual_outcomes) if annual_outcomes else 0,
            'qalys_per_treatment': (results['total_qalys_gained'] / results['treated_population'] 
                                   if results['treated_population'] > 0 else 0)
        }