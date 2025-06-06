"""
Combination therapy module for Enhancement Engine.

This module provides comprehensive multi-gene therapeutic strategies including:
- Multi-target optimization
- Synergistic effect prediction
- Sequential vs simultaneous delivery
- Risk-benefit optimization
- Personalized combination design
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any, Set
from dataclasses import dataclass
from itertools import combinations, permutations
from collections import defaultdict

from ..models.therapeutic_data_classes import (
    CombinationStrategy,
    TherapeuticTarget,
    SynergyScore,
    DeliverySchedule,
    RiskProfile,
    CombinationResult,
)
from ..models.disease_constants import GENE_INTERACTIONS, THERAPEUTIC_SYNERGIES


@dataclass
class CombinationTarget:
    """Single target in combination therapy."""

    gene_name: str
    therapeutic_strategy: str
    priority: int  # 1 = highest priority
    delivery_requirements: Dict[str, Any]
    expected_effect_size: float
    safety_profile: Dict[str, Any]

    @property
    def risk_score(self) -> float:
        """Calculate risk score for this target."""
        return self.safety_profile.get("overall_risk", 0.5)


@dataclass
class CombinationOptimization:
    """Results from combination optimization."""

    optimal_targets: List[CombinationTarget]
    delivery_schedule: DeliverySchedule
    predicted_synergy: float
    risk_benefit_ratio: float
    implementation_complexity: str
    monitoring_requirements: List[str]


class CombinationTherapyDesigner:
    """Main class for designing combination therapeutic strategies."""

    def __init__(self):
        """Initialize combination therapy designer."""
        self.logger = logging.getLogger(__name__)
        self._load_gene_interactions()
        self._load_synergy_models()

    def _load_gene_interactions(self) -> None:
        """Load gene interaction networks."""
        self.gene_interactions = {
            # Rheumatoid arthritis gene network
            "PTPN22": {
                "HLA-DRB1": {"type": "epistatic", "strength": 0.7},
                "STAT4": {"type": "additive", "strength": 0.4},
                "TNFAIP3": {"type": "compensatory", "strength": 0.3},
            },
            "HLA-DRB1": {
                "PTPN22": {"type": "epistatic", "strength": 0.7},
                "STAT4": {"type": "modifying", "strength": 0.5},
                "PADI4": {"type": "synergistic", "strength": 0.6},
            },
            "STAT4": {
                "PTPN22": {"type": "additive", "strength": 0.4},
                "HLA-DRB1": {"type": "modifying", "strength": 0.5},
                "IRF4": {"type": "pathway", "strength": 0.8},
            },
        }

        # Pathway-level interactions
        self.pathway_interactions = {
            "immune_regulation": ["PTPN22", "STAT4", "IRF4"],
            "antigen_presentation": ["HLA-DRB1", "HLA-DQB1", "PADI4"],
            "inflammation": ["TNFAIP3", "TRAF1", "IL2RA"],
            "autoimmunity": ["PTPN22", "HLA-DRB1", "CTLA4"],
        }

    def _load_synergy_models(self) -> None:
        """Load models for predicting therapeutic synergy."""
        self.synergy_models = {
            # Gene pair synergies for RA treatment
            ("PTPN22", "HLA-DRB1"): {
                "synergy_type": "multiplicative",
                "synergy_factor": 1.8,
                "confidence": 0.8,
                "mechanism": "Combined T-cell regulation and antigen presentation",
            },
            ("HLA-DRB1", "STAT4"): {
                "synergy_type": "additive_plus",
                "synergy_factor": 1.3,
                "confidence": 0.6,
                "mechanism": "Antigen presentation and T-cell differentiation",
            },
            ("PTPN22", "STAT4"): {
                "synergy_type": "additive",
                "synergy_factor": 1.2,
                "confidence": 0.7,
                "mechanism": "Complementary T-cell regulation pathways",
            },
        }

        # Multi-gene combinations (3+ genes)
        self.multi_gene_synergies = {
            ("PTPN22", "HLA-DRB1", "STAT4"): {
                "synergy_factor": 2.5,
                "complexity_penalty": 0.3,
                "coordination_requirement": "high",
            }
        }

    def design_combination_strategy(
        self,
        candidate_targets: List[Dict[str, Any]],
        disease_context: str = "rheumatoid_arthritis",
        patient_profile: Optional[Dict[str, Any]] = None,
    ) -> CombinationStrategy:
        """
        Design optimal combination therapeutic strategy.

        Args:
            candidate_targets: List of potential therapeutic targets
            disease_context: Disease being treated
            patient_profile: Patient-specific information

        Returns:
            Optimized combination strategy
        """
        # Convert candidate targets to CombinationTarget objects
        targets = [
            self._create_combination_target(target) for target in candidate_targets
        ]

        # Generate all possible combinations
        combinations_to_evaluate = self._generate_target_combinations(targets)

        # Evaluate each combination
        evaluated_combinations = []
        for combo in combinations_to_evaluate:
            evaluation = self._evaluate_combination(
                combo, disease_context, patient_profile
            )
            if evaluation["viable"]:
                evaluated_combinations.append(evaluation)

        if not evaluated_combinations:
            self.logger.warning("No viable combinations found")
            return self._create_single_target_strategy(targets[0] if targets else None)

        # Select optimal combination
        optimal_combination = self._select_optimal_combination(evaluated_combinations)

        # Optimize delivery strategy
        delivery_optimization = self._optimize_combination_delivery(optimal_combination)

        # Create final strategy
        return CombinationStrategy(
            targets=optimal_combination["targets"],
            synergy_score=optimal_combination["synergy_score"],
            delivery_schedule=delivery_optimization["schedule"],
            risk_profile=optimal_combination["risk_profile"],
            expected_efficacy=optimal_combination["expected_efficacy"],
            implementation_plan=self._create_implementation_plan(
                optimal_combination, delivery_optimization
            ),
            monitoring_strategy=self._design_monitoring_strategy(optimal_combination),
        )

    def _create_combination_target(
        self, target_data: Dict[str, Any]
    ) -> CombinationTarget:
        """Create CombinationTarget from target data."""
        return CombinationTarget(
            gene_name=target_data["gene_name"],
            therapeutic_strategy=target_data.get("strategy", "base_editing"),
            priority=target_data.get("priority", 2),
            delivery_requirements=target_data.get("delivery_requirements", {}),
            expected_effect_size=target_data.get("effect_size", 0.5),
            safety_profile=target_data.get("safety_profile", {"overall_risk": 0.3}),
        )

    def _generate_target_combinations(
        self, targets: List[CombinationTarget]
    ) -> List[List[CombinationTarget]]:
        """Generate all viable target combinations."""
        combinations_list = []

        # Single targets
        for target in targets:
            combinations_list.append([target])

        # Pairs
        for combo in combinations(targets, 2):
            combinations_list.append(list(combo))

        # Triplets (limited to avoid complexity explosion)
        if len(targets) >= 3:
            for combo in combinations(targets, 3):
                combinations_list.append(list(combo))

        # Quadruplets only for very specific cases
        if len(targets) >= 4 and len(targets) <= 5:
            for combo in combinations(targets, 4):
                # Only include if all targets are low-risk
                if all(target.risk_score < 0.3 for target in combo):
                    combinations_list.append(list(combo))

        return combinations_list

    def _evaluate_combination(
        self,
        targets: List[CombinationTarget],
        disease_context: str,
        patient_profile: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Evaluate a specific combination of targets."""
        gene_names = [target.gene_name for target in targets]

        # Calculate synergy
        synergy_score = self._calculate_synergy_score(targets)

        # Calculate combined efficacy
        combined_efficacy = self._calculate_combined_efficacy(targets, synergy_score)

        # Assess combined safety
        combined_risk = self._assess_combined_risk(targets)

        # Check feasibility
        feasibility = self._assess_combination_feasibility(targets)

        # Calculate implementation complexity
        complexity = self._calculate_implementation_complexity(targets)

        # Patient-specific adjustments
        if patient_profile:
            combined_efficacy *= self._get_patient_efficacy_modifier(
                targets, patient_profile
            )
            combined_risk *= self._get_patient_risk_modifier(targets, patient_profile)

        # Overall viability assessment
        viable = (
            combined_efficacy >= 0.3
            and combined_risk <= 0.7
            and feasibility >= 0.4
            and complexity != "impossible"
        )

        return {
            "targets": targets,
            "gene_names": gene_names,
            "synergy_score": synergy_score,
            "expected_efficacy": combined_efficacy,
            "risk_profile": {
                "overall_risk": combined_risk,
                "individual_risks": [target.risk_score for target in targets],
            },
            "feasibility": feasibility,
            "complexity": complexity,
            "viable": viable,
            "benefit_risk_ratio": combined_efficacy / max(0.1, combined_risk),
        }

    def _calculate_synergy_score(self, targets: List[CombinationTarget]) -> float:
        """Calculate synergy score for combination of targets."""
        if len(targets) == 1:
            return 1.0  # No synergy for single target

        gene_names = [target.gene_name for target in targets]
        total_synergy = 1.0

        # Pairwise synergies
        for i in range(len(gene_names)):
            for j in range(i + 1, len(gene_names)):
                gene_pair = tuple(sorted([gene_names[i], gene_names[j]]))

                if gene_pair in self.synergy_models:
                    synergy_data = self.synergy_models[gene_pair]
                    synergy_factor = synergy_data["synergy_factor"]
                    confidence = synergy_data["confidence"]

                    # Apply synergy with confidence weighting
                    weighted_synergy = 1.0 + (synergy_factor - 1.0) * confidence
                    total_synergy *= weighted_synergy
                else:
                    # Default interaction based on gene network
                    default_synergy = self._estimate_default_synergy(
                        gene_names[i], gene_names[j]
                    )
                    total_synergy *= default_synergy

        # Multi-gene synergies
        if len(gene_names) >= 3:
            gene_tuple = tuple(sorted(gene_names))
            if gene_tuple in self.multi_gene_synergies:
                multi_synergy = self.multi_gene_synergies[gene_tuple]
                synergy_boost = multi_synergy["synergy_factor"]
                complexity_penalty = multi_synergy["complexity_penalty"]

                total_synergy *= synergy_boost * (1.0 - complexity_penalty)

        # Diminishing returns for large combinations
        if len(targets) > 3:
            diminishing_factor = 0.9 ** (len(targets) - 3)
            total_synergy *= diminishing_factor

        return min(3.0, max(0.5, total_synergy))  # Cap synergy between 0.5x and 3.0x

    def _estimate_default_synergy(self, gene1: str, gene2: str) -> float:
        """Estimate synergy between genes based on known interactions."""
        # Check direct interactions
        if gene1 in self.gene_interactions:
            if gene2 in self.gene_interactions[gene1]:
                interaction = self.gene_interactions[gene1][gene2]

                interaction_factors = {
                    "synergistic": 1.4,
                    "additive": 1.2,
                    "epistatic": 1.3,
                    "compensatory": 1.1,
                    "antagonistic": 0.9,
                }

                base_factor = interaction_factors.get(interaction["type"], 1.1)
                strength = interaction["strength"]

                return 1.0 + (base_factor - 1.0) * strength

        # Check pathway overlap
        shared_pathways = 0
        for pathway, genes in self.pathway_interactions.items():
            if gene1 in genes and gene2 in genes:
                shared_pathways += 1

        if shared_pathways > 0:
            return 1.0 + (shared_pathways * 0.1)  # 10% boost per shared pathway

        return 1.05  # Minimal positive interaction by default

    def _calculate_combined_efficacy(
        self, targets: List[CombinationTarget], synergy_score: float
    ) -> float:
        """Calculate combined efficacy of target combination."""
        # Base efficacy from individual targets
        individual_efficacies = [target.expected_effect_size for target in targets]

        if len(individual_efficacies) == 1:
            return individual_efficacies[0]

        # Combine efficacies (geometric mean with synergy adjustment)
        geometric_mean = np.prod(individual_efficacies) ** (
            1 / len(individual_efficacies)
        )

        # Apply synergy
        combined_efficacy = geometric_mean * synergy_score

        # Apply diminishing returns
        max_possible_efficacy = 0.95  # Never achieve 100% efficacy
        combined_efficacy = max_possible_efficacy * (1 - np.exp(-2 * combined_efficacy))

        return min(max_possible_efficacy, combined_efficacy)

    def _assess_combined_risk(self, targets: List[CombinationTarget]) -> float:
        """Assess combined risk of target combination."""
        individual_risks = [target.risk_score for target in targets]

        if len(individual_risks) == 1:
            return individual_risks[0]

        # Risk combination is more than additive due to interactions
        combined_risk = 1.0 - np.prod([1.0 - risk for risk in individual_risks])

        # Add complexity penalty
        complexity_penalty = (len(targets) - 1) * 0.05  # 5% per additional target
        combined_risk += complexity_penalty

        # Add interaction risk
        interaction_risk = self._calculate_interaction_risk(targets)
        combined_risk += interaction_risk

        return min(1.0, combined_risk)

    def _calculate_interaction_risk(self, targets: List[CombinationTarget]) -> float:
        """Calculate additional risk from target interactions."""
        gene_names = [target.gene_name for target in targets]
        interaction_risk = 0.0

        # Check for known problematic interactions
        problematic_combinations = {
            ("PTPN22", "STAT4"): 0.05,  # Both affect T-cell function
            ("HLA-DRB1", "HLA-DQB1"): 0.08,  # Both MHC class II
        }

        for combo, risk in problematic_combinations.items():
            if all(gene in gene_names for gene in combo):
                interaction_risk += risk

        # General pathway overload risk
        pathway_loads = defaultdict(int)
        for gene in gene_names:
            for pathway, pathway_genes in self.pathway_interactions.items():
                if gene in pathway_genes:
                    pathway_loads[pathway] += 1

        # Penalize pathways with multiple targets
        for pathway, load in pathway_loads.items():
            if load > 1:
                interaction_risk += (
                    load - 1
                ) * 0.03  # 3% per additional gene in pathway

        return interaction_risk

    def _assess_combination_feasibility(
        self, targets: List[CombinationTarget]
    ) -> float:
        """Assess technical feasibility of combination."""
        # Delivery compatibility
        delivery_compatibility = self._assess_delivery_compatibility(targets)

        # Manufacturing feasibility
        manufacturing_feasibility = self._assess_manufacturing_feasibility(targets)

        # Regulatory feasibility
        regulatory_feasibility = self._assess_regulatory_feasibility(targets)

        # Timing feasibility
        timing_feasibility = self._assess_timing_feasibility(targets)

        # Combined feasibility score
        feasibility_scores = [
            delivery_compatibility,
            manufacturing_feasibility,
            regulatory_feasibility,
            timing_feasibility,
        ]

        return np.mean(feasibility_scores)

    def _assess_delivery_compatibility(self, targets: List[CombinationTarget]) -> float:
        """Assess compatibility of delivery requirements."""
        delivery_methods = []
        target_tissues = []

        for target in targets:
            req = target.delivery_requirements
            delivery_methods.append(req.get("method", "unknown"))
            target_tissues.append(req.get("tissue", "systemic"))

        # Check method compatibility
        unique_methods = set(delivery_methods)
        method_compatibility = (
            1.0 - (len(unique_methods) - 1) * 0.2
        )  # 20% penalty per additional method

        # Check tissue compatibility
        unique_tissues = set(target_tissues)
        tissue_compatibility = (
            1.0 - (len(unique_tissues) - 1) * 0.15
        )  # 15% penalty per additional tissue

        return max(0.0, (method_compatibility + tissue_compatibility) / 2)

    def _optimize_combination_delivery(
        self, combination_eval: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Optimize delivery strategy for combination therapy."""
        targets = combination_eval["targets"]

        # Determine delivery approach
        if len(targets) == 1:
            approach = "single_delivery"
        elif self._can_deliver_simultaneously(targets):
            approach = "simultaneous"
        else:
            approach = "sequential"

        # Create delivery schedule
        if approach == "simultaneous":
            schedule = self._create_simultaneous_schedule(targets)
        else:
            schedule = self._create_sequential_schedule(targets)

        # Optimize timing
        optimized_schedule = self._optimize_delivery_timing(schedule, targets)

        return {
            "approach": approach,
            "schedule": optimized_schedule,
            "rationale": self._explain_delivery_rationale(approach, targets),
        }

    def _can_deliver_simultaneously(self, targets: List[CombinationTarget]) -> bool:
        """Check if targets can be delivered simultaneously."""
        # Check delivery method compatibility
        methods = [
            target.delivery_requirements.get("method", "LNP") for target in targets
        ]
        if len(set(methods)) > 1:
            return False

        # Check safety profile compatibility
        total_risk = sum(target.risk_score for target in targets)
        if total_risk > 0.8:
            return False

        # Check for known incompatibilities
        gene_names = [target.gene_name for target in targets]
        incompatible_pairs = [("HLA-DRB1", "HLA-DQB1")]  # Example

        for pair in incompatible_pairs:
            if all(gene in gene_names for gene in pair):
                return False

        return True

    def _create_sequential_schedule(
        self, targets: List[CombinationTarget]
    ) -> DeliverySchedule:
        """Create sequential delivery schedule."""
        # Sort targets by priority and safety
        sorted_targets = sorted(targets, key=lambda t: (t.priority, t.risk_score))

        schedule_items = []
        current_time = 0

        for i, target in enumerate(sorted_targets):
            # Determine delay between deliveries
            if i == 0:
                delay = 0
            else:
                # Base delay on safety and interaction concerns
                prev_target = sorted_targets[i - 1]
                delay = self._calculate_delivery_delay(prev_target, target)

            current_time += delay

            schedule_items.append(
                {
                    "target": target,
                    "delivery_time_weeks": current_time,
                    "monitoring_period_weeks": 4,  # Monitor for 4 weeks after each delivery
                    "safety_checkpoints": [1, 2, 4],  # Weeks 1, 2, and 4 post-delivery
                }
            )

        return DeliverySchedule(
            delivery_type="sequential",
            schedule_items=schedule_items,
            total_duration_weeks=current_time + 4,  # Plus final monitoring
            monitoring_requirements=self._create_sequential_monitoring(schedule_items),
        )

    def _calculate_delivery_delay(
        self, prev_target: CombinationTarget, current_target: CombinationTarget
    ) -> int:
        """Calculate optimal delay between sequential deliveries."""
        # Base delay
        base_delay = 2  # 2 weeks minimum

        # Safety-based adjustment
        safety_adjustment = int(
            prev_target.risk_score * 4
        )  # Up to 4 weeks for high-risk

        # Interaction-based adjustment
        gene_interaction = self.gene_interactions.get(prev_target.gene_name, {})
        if current_target.gene_name in gene_interaction:
            interaction_type = gene_interaction[current_target.gene_name]["type"]
            if interaction_type in ["antagonistic", "compensatory"]:
                base_delay += 2  # Extra delay for potentially interfering targets

        return base_delay + safety_adjustment

    def _select_optimal_combination(
        self, evaluated_combinations: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Select optimal combination from evaluated options."""
        # Score combinations
        for combo in evaluated_combinations:
            score = self._calculate_combination_score(combo)
            combo["overall_score"] = score

        # Sort by score
        evaluated_combinations.sort(key=lambda x: x["overall_score"], reverse=True)

        # Return best combination
        return evaluated_combinations[0]

    def _calculate_combination_score(self, combination: Dict[str, Any]) -> float:
        """Calculate overall score for combination."""
        efficacy = combination["expected_efficacy"]
        risk = combination["risk_profile"]["overall_risk"]
        feasibility = combination["feasibility"]
        synergy = combination["synergy_score"]

        # Weighted score
        score = (
            efficacy * 0.4
            + (1.0 - risk) * 0.3
            + feasibility * 0.2
            + min(synergy / 2.0, 1.0) * 0.1  # Normalize synergy contribution
        )

        return score

    def generate_combination_report(
        self, combination_strategy: CombinationStrategy
    ) -> str:
        """Generate comprehensive combination therapy report."""
        report = f"""
COMBINATION THERAPY STRATEGY REPORT
=================================

Number of Targets: {len(combination_strategy.targets)}
Targets: {', '.join([t.gene_name for t in combination_strategy.targets])}
Synergy Score: {combination_strategy.synergy_score:.2f}
Expected Efficacy: {combination_strategy.expected_efficacy:.2f}
Overall Risk: {combination_strategy.risk_profile.get('overall_risk', 'N/A')}

DELIVERY STRATEGY:
{combination_strategy.delivery_schedule.delivery_type.upper()}
Total Duration: {combination_strategy.delivery_schedule.total_duration_weeks} weeks

TARGET DETAILS:
"""

        for i, target in enumerate(combination_strategy.targets, 1):
            report += f"""
Target {i}: {target.gene_name}
- Strategy: {target.therapeutic_strategy}
- Priority: {target.priority}
- Effect Size: {target.expected_effect_size:.2f}
- Risk Score: {target.risk_score:.2f}
"""

        report += f"""
IMPLEMENTATION PLAN:
{chr(10).join('- ' + step for step in combination_strategy.implementation_plan)}

MONITORING STRATEGY:
{chr(10).join('- ' + req for req in combination_strategy.monitoring_strategy)}
        """

        return report.strip()

    # Additional helper methods would be implemented here...
    def _create_implementation_plan(
        self, combination: Dict[str, Any], delivery_opt: Dict[str, Any]
    ) -> List[str]:
        """Create implementation plan for combination therapy."""
        plan = [
            "Complete preclinical validation for all targets",
            "Establish manufacturing protocols for combination",
            "Obtain regulatory approvals for combination approach",
            f"Implement {delivery_opt['approach']} delivery strategy",
            "Execute comprehensive monitoring protocol",
        ]

        if delivery_opt["approach"] == "sequential":
            plan.append("Monitor for cumulative effects between deliveries")

        return plan

    def _design_monitoring_strategy(self, combination: Dict[str, Any]) -> List[str]:
        """Design monitoring strategy for combination therapy."""
        monitoring = [
            "Individual target confirmation for each gene",
            "Synergistic effect assessment",
            "Combined safety monitoring",
            "Efficacy endpoint evaluation",
        ]

        # Add gene-specific monitoring
        for target in combination["targets"]:
            monitoring.append(f"{target.gene_name}-specific functional assays")

        return monitoring
