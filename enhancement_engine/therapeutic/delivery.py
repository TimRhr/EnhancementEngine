"""
Therapeutic delivery system module for Enhancement Engine.

This module provides comprehensive delivery strategies for therapeutic CRISPR
interventions including tissue-specific targeting, dosing optimization,
and delivery vehicle selection.
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass
from enum import Enum

from ..models.therapeutic_data_classes import (
    DeliveryMethod, TargetTissue, TherapeuticStrategy,
    DeliveryEfficiency, TissueDistribution
)
from ..models.disease_constants import THERAPEUTIC_TARGETS


class DeliveryVehicle(Enum):
    """Types of delivery vehicles."""
    LNP = "lipid_nanoparticles"
    AAV = "adeno_associated_virus"
    LENTIVIRUS = "lentiviral_vector"
    ELECTROPORATION = "electroporation"
    MICROINJECTION = "direct_injection"
    LIPOSOMES = "liposomal_delivery"
    PROTEIN_DELIVERY = "ribonucleoprotein"
    EXOSOMES = "exosome_mediated"


@dataclass
class DeliveryParameters:
    """Parameters for therapeutic delivery."""
    vehicle: DeliveryVehicle
    target_tissue: TargetTissue
    dose_mg_kg: float
    administration_route: str
    treatment_schedule: str
    bioavailability: float
    half_life_hours: float
    
    @property
    def effective_dose(self) -> float:
        """Calculate effective dose accounting for bioavailability."""
        return self.dose_mg_kg * self.bioavailability


class TherapeuticDeliverySystem:
    """Main class for therapeutic delivery optimization."""
    
    def __init__(self):
        """Initialize delivery system."""
        self.logger = logging.getLogger(__name__)
        self._load_delivery_profiles()
    
    def _load_delivery_profiles(self) -> None:
        """Load delivery vehicle profiles."""
        self.delivery_profiles = {
            DeliveryVehicle.LNP: {
                'tissue_tropism': {
                    TargetTissue.LIVER: 0.8,
                    TargetTissue.SPLEEN: 0.6,
                    TargetTissue.BONE_MARROW: 0.4,
                    TargetTissue.SYNOVIAL: 0.2,
                    TargetTissue.BRAIN: 0.1
                },
                'efficiency': 0.7,
                'immunogenicity': 0.3,
                'stability_hours': 24,
                'cost_factor': 1.0
            },
            DeliveryVehicle.AAV: {
                'tissue_tropism': {
                    TargetTissue.MUSCLE: 0.9,
                    TargetTissue.LIVER: 0.7,
                    TargetTissue.BRAIN: 0.6,
                    TargetTissue.SYNOVIAL: 0.4,
                    TargetTissue.BONE_MARROW: 0.3
                },
                'efficiency': 0.8,
                'immunogenicity': 0.6,
                'stability_hours': 72,
                'cost_factor': 2.5
            },
            DeliveryVehicle.ELECTROPORATION: {
                'tissue_tropism': {
                    TargetTissue.MUSCLE: 0.9,
                    TargetTissue.SKIN: 0.8,
                    TargetTissue.SYNOVIAL: 0.7,
                    TargetTissue.LIVER: 0.2,
                    TargetTissue.BRAIN: 0.1
                },
                'efficiency': 0.6,
                'immunogenicity': 0.1,
                'stability_hours': 6,
                'cost_factor': 0.5
            }
        }
    
    def optimize_delivery_strategy(self, target_gene: str, therapeutic_strategy: TherapeuticStrategy,
                                 target_tissue: TargetTissue) -> DeliveryParameters:
        """
        Optimize delivery strategy for specific therapeutic application.
        
        Args:
            target_gene: Gene to be targeted
            therapeutic_strategy: Type of therapeutic intervention
            target_tissue: Primary target tissue
            
        Returns:
            Optimized delivery parameters
        """
        # Get gene-specific requirements
        gene_requirements = self._get_gene_delivery_requirements(target_gene)
        
        # Score delivery vehicles
        vehicle_scores = {}
        for vehicle, profile in self.delivery_profiles.items():
            score = self._calculate_delivery_score(
                vehicle, profile, target_tissue, therapeutic_strategy, gene_requirements
            )
            vehicle_scores[vehicle] = score
        
        # Select best vehicle
        best_vehicle = max(vehicle_scores.keys(), key=lambda x: vehicle_scores[x])
        
        # Optimize parameters
        optimized_params = self._optimize_delivery_parameters(
            best_vehicle, target_tissue, therapeutic_strategy
        )
        
        self.logger.info(f"Optimized delivery for {target_gene}: {best_vehicle.value}")
        return optimized_params
    
    def _get_gene_delivery_requirements(self, gene: str) -> Dict[str, Any]:
        """Get delivery requirements specific to gene."""
        requirements = {
            'PTPN22': {
                'primary_tissue': TargetTissue.IMMUNE_CELLS,
                'required_efficiency': 0.6,
                'systemic_delivery': True,
                'persistence_needed': False
            },
            'HLA-DRB1': {
                'primary_tissue': TargetTissue.ANTIGEN_PRESENTING,
                'required_efficiency': 0.8,
                'systemic_delivery': True,
                'persistence_needed': True
            },
            'STAT4': {
                'primary_tissue': TargetTissue.T_CELLS,
                'required_efficiency': 0.5,
                'systemic_delivery': True,
                'persistence_needed': False
            }
        }
        
        return requirements.get(gene, {
            'primary_tissue': TargetTissue.SYSTEMIC,
            'required_efficiency': 0.5,
            'systemic_delivery': True,
            'persistence_needed': False
        })
    
    def _calculate_delivery_score(self, vehicle: DeliveryVehicle, profile: Dict,
                                target_tissue: TargetTissue, strategy: TherapeuticStrategy,
                                requirements: Dict) -> float:
        """Calculate delivery score for vehicle-tissue combination."""
        # Base tropism score
        tropism_score = profile['tissue_tropism'].get(target_tissue, 0.1)
        
        # Efficiency weight
        efficiency_weight = 0.3
        efficiency_score = profile['efficiency'] * efficiency_weight
        
        # Immunogenicity penalty (lower is better)
        immunogenicity_penalty = profile['immunogenicity'] * 0.2
        
        # Strategy-specific adjustments
        strategy_bonus = 0.0
        if strategy == TherapeuticStrategy.BASE_EDITING:
            # Base editing needs sustained expression
            if profile['stability_hours'] > 24:
                strategy_bonus = 0.1
        elif strategy == TherapeuticStrategy.GENE_REPLACEMENT:
            # Gene replacement needs high efficiency
            if profile['efficiency'] > 0.7:
                strategy_bonus = 0.2
        
        # Cost consideration (lower cost is better)
        cost_penalty = (profile['cost_factor'] - 1.0) * 0.1
        
        total_score = (tropism_score + efficiency_score + strategy_bonus 
                      - immunogenicity_penalty - cost_penalty)
        
        return max(0.0, total_score)
    
    def _optimize_delivery_parameters(self, vehicle: DeliveryVehicle, 
                                    target_tissue: TargetTissue,
                                    strategy: TherapeuticStrategy) -> DeliveryParameters:
        """Optimize specific delivery parameters."""
        profile = self.delivery_profiles[vehicle]
        
        # Base dose optimization
        base_dose = self._calculate_optimal_dose(vehicle, target_tissue, strategy)
        
        # Administration route
        admin_route = self._select_administration_route(vehicle, target_tissue)
        
        # Treatment schedule
        schedule = self._optimize_treatment_schedule(vehicle, strategy)
        
        return DeliveryParameters(
            vehicle=vehicle,
            target_tissue=target_tissue,
            dose_mg_kg=base_dose,
            administration_route=admin_route,
            treatment_schedule=schedule,
            bioavailability=profile['efficiency'],
            half_life_hours=profile['stability_hours']
        )
    
    def _calculate_optimal_dose(self, vehicle: DeliveryVehicle, 
                              target_tissue: TargetTissue,
                              strategy: TherapeuticStrategy) -> float:
        """Calculate optimal dose for delivery combination."""
        # Base doses by vehicle type (mg/kg)
        base_doses = {
            DeliveryVehicle.LNP: 1.0,
            DeliveryVehicle.AAV: 1e12,  # viral genomes/kg
            DeliveryVehicle.ELECTROPORATION: 0.1,
            DeliveryVehicle.LENTIVIRUS: 1e8,
            DeliveryVehicle.PROTEIN_DELIVERY: 5.0
        }
        
        base_dose = base_doses.get(vehicle, 1.0)
        
        # Tissue-specific adjustments
        tissue_factors = {
            TargetTissue.BRAIN: 0.5,      # Lower dose for BBB
            TargetTissue.LIVER: 1.2,      # Higher clearance
            TargetTissue.SYNOVIAL: 2.0,   # Local delivery
            TargetTissue.IMMUNE_CELLS: 0.8
        }
        
        tissue_factor = tissue_factors.get(target_tissue, 1.0)
        
        # Strategy-specific adjustments
        strategy_factors = {
            TherapeuticStrategy.BASE_EDITING: 1.5,     # Needs more protein
            TherapeuticStrategy.GENE_REPLACEMENT: 2.0, # Needs high efficiency
            TherapeuticStrategy.GENE_SILENCING: 0.8    # Less material needed
        }
        
        strategy_factor = strategy_factors.get(strategy, 1.0)
        
        return base_dose * tissue_factor * strategy_factor
    
    def _select_administration_route(self, vehicle: DeliveryVehicle,
                                   target_tissue: TargetTissue) -> str:
        """Select optimal administration route."""
        routes = {
            (DeliveryVehicle.LNP, TargetTissue.LIVER): "intravenous",
            (DeliveryVehicle.LNP, TargetTissue.SYNOVIAL): "intra-articular",
            (DeliveryVehicle.AAV, TargetTissue.MUSCLE): "intramuscular",
            (DeliveryVehicle.AAV, TargetTissue.BRAIN): "intrathecal",
            (DeliveryVehicle.ELECTROPORATION, TargetTissue.MUSCLE): "intramuscular",
            (DeliveryVehicle.ELECTROPORATION, TargetTissue.SYNOVIAL): "intra-articular"
        }
        
        return routes.get((vehicle, target_tissue), "intravenous")
    
    def _optimize_treatment_schedule(self, vehicle: DeliveryVehicle,
                                   strategy: TherapeuticStrategy) -> str:
        """Optimize treatment schedule."""
        # Base schedules by vehicle
        base_schedules = {
            DeliveryVehicle.LNP: "single_dose",
            DeliveryVehicle.AAV: "single_dose",
            DeliveryVehicle.ELECTROPORATION: "3_doses_weekly",
            DeliveryVehicle.PROTEIN_DELIVERY: "weekly_x4"
        }
        
        base_schedule = base_schedules.get(vehicle, "single_dose")
        
        # Strategy adjustments
        if strategy == TherapeuticStrategy.GENE_SILENCING:
            # May need repeated dosing
            if vehicle == DeliveryVehicle.LNP:
                return "monthly_x3"
        elif strategy == TherapeuticStrategy.BASE_EDITING:
            # Single treatment often sufficient
            return "single_dose"
        
        return base_schedule
    
    def predict_tissue_distribution(self, delivery_params: DeliveryParameters) -> TissueDistribution:
        """
        Predict tissue distribution of therapeutic payload.
        
        Args:
            delivery_params: Delivery parameters
            
        Returns:
            Predicted tissue distribution
        """
        profile = self.delivery_profiles[delivery_params.vehicle]
        
        # Base distribution from vehicle profile
        base_distribution = profile['tissue_tropism'].copy()
        
        # Adjust for administration route
        route_adjustments = {
            "intra-articular": {
                TargetTissue.SYNOVIAL: 2.0,
                TargetTissue.LIVER: 0.3
            },
            "intrathecal": {
                TargetTissue.BRAIN: 3.0,
                TargetTissue.LIVER: 0.1
            },
            "intramuscular": {
                TargetTissue.MUSCLE: 2.0
            }
        }
        
        adjustments = route_adjustments.get(delivery_params.administration_route, {})
        
        # Apply adjustments
        adjusted_distribution = {}
        for tissue, base_level in base_distribution.items():
            adjustment = adjustments.get(tissue, 1.0)
            adjusted_distribution[tissue] = min(1.0, base_level * adjustment)
        
        # Calculate efficiency per tissue
        tissue_efficiencies = {}
        for tissue, distribution in adjusted_distribution.items():
            efficiency = distribution * delivery_params.bioavailability
            tissue_efficiencies[tissue] = efficiency
        
        return TissueDistribution(
            distribution_profile=adjusted_distribution,
            efficiency_per_tissue=tissue_efficiencies,
            peak_concentration_hours=delivery_params.half_life_hours / 2,
            clearance_half_life=delivery_params.half_life_hours
        )
    
    def calculate_delivery_efficiency(self, delivery_params: DeliveryParameters,
                                    target_cells: int = 1e6) -> DeliveryEfficiency:
        """
        Calculate overall delivery efficiency.
        
        Args:
            delivery_params: Delivery parameters
            target_cells: Number of target cells
            
        Returns:
            Delivery efficiency metrics
        """
        # Get tissue distribution
        tissue_dist = self.predict_tissue_distribution(delivery_params)
        
        # Calculate cells reached
        primary_tissue_efficiency = tissue_dist.efficiency_per_tissue.get(
            delivery_params.target_tissue, 0.1
        )
        
        cells_reached = int(target_cells * primary_tissue_efficiency)
        
        # Calculate successful edits (depends on CRISPR efficiency)
        editing_efficiency = 0.6  # Assume 60% editing efficiency
        successful_edits = int(cells_reached * editing_efficiency)
        
        # Calculate persistence
        vehicle_profile = self.delivery_profiles[delivery_params.vehicle]
        persistence_score = min(1.0, vehicle_profile['stability_hours'] / 48)
        
        return DeliveryEfficiency(
            cells_reached=cells_reached,
            successful_edits=successful_edits,
            editing_percentage=successful_edits / target_cells * 100,
            tissue_penetration=primary_tissue_efficiency,
            persistence_score=persistence_score,
            off_tissue_distribution=self._calculate_off_tissue_exposure(tissue_dist)
        )
    
    def _calculate_off_tissue_exposure(self, tissue_dist: TissueDistribution) -> float:
        """Calculate exposure to non-target tissues."""
        primary_tissue = max(tissue_dist.efficiency_per_tissue.values())
        total_exposure = sum(tissue_dist.efficiency_per_tissue.values())
        
        off_target_exposure = total_exposure - primary_tissue
        return off_target_exposure / total_exposure if total_exposure > 0 else 0.0
    
    def optimize_combination_delivery(self, targets: List[str],
                                    strategies: List[TherapeuticStrategy]) -> List[DeliveryParameters]:
        """
        Optimize delivery for combination therapy.
        
        Args:
            targets: List of target genes
            strategies: Corresponding therapeutic strategies
            
        Returns:
            List of optimized delivery parameters
        """
        if len(targets) != len(strategies):
            raise ValueError("Number of targets must match number of strategies")
        
        optimized_deliveries = []
        
        for target, strategy in zip(targets, strategies):
            # Get tissue requirements
            requirements = self._get_gene_delivery_requirements(target)
            primary_tissue = requirements.get('primary_tissue', TargetTissue.SYSTEMIC)
            
            # Optimize individual delivery
            delivery_params = self.optimize_delivery_strategy(target, strategy, primary_tissue)
            optimized_deliveries.append(delivery_params)
        
        # Check for compatibility and optimize timing
        optimized_deliveries = self._optimize_combination_timing(optimized_deliveries)
        
        return optimized_deliveries
    
    def _optimize_combination_timing(self, deliveries: List[DeliveryParameters]) -> List[DeliveryParameters]:
        """Optimize timing for combination therapies."""
        # Simple strategy: space out different vehicles to avoid interference
        vehicle_counts = {}
        for delivery in deliveries:
            vehicle_counts[delivery.vehicle] = vehicle_counts.get(delivery.vehicle, 0) + 1
        
        # If multiple same vehicles, add delays
        vehicle_delays = {}
        for delivery in deliveries:
            if vehicle_counts[delivery.vehicle] > 1:
                delay_count = vehicle_delays.get(delivery.vehicle, 0)
                vehicle_delays[delivery.vehicle] = delay_count + 1
                
                if delay_count > 0:
                    # Modify schedule to include delay
                    original_schedule = delivery.treatment_schedule
                    delivery.treatment_schedule = f"{original_schedule}_delay_{delay_count}weeks"
        
        return deliveries
    
    def assess_delivery_safety(self, delivery_params: DeliveryParameters) -> Dict[str, Any]:
        """
        Assess safety of delivery strategy.
        
        Args:
            delivery_params: Delivery parameters to assess
            
        Returns:
            Safety assessment results
        """
        vehicle_profile = self.delivery_profiles[delivery_params.vehicle]
        
        # Base safety scores
        immunogenicity_risk = vehicle_profile['immunogenicity']
        
        # Dose-dependent toxicity
        dose_toxicity = min(1.0, delivery_params.dose_mg_kg / 10.0)  # Simplified
        
        # Route-specific risks
        route_risks = {
            "intravenous": 0.3,
            "intra-articular": 0.2,
            "intramuscular": 0.1,
            "intrathecal": 0.8  # Higher risk for CNS
        }
        
        route_risk = route_risks.get(delivery_params.administration_route, 0.3)
        
        # Overall safety score (lower is safer)
        overall_risk = (immunogenicity_risk + dose_toxicity + route_risk) / 3
        
        return {
            'overall_risk_score': overall_risk,
            'immunogenicity_risk': immunogenicity_risk,
            'dose_toxicity': dose_toxicity,
            'route_risk': route_risk,
            'safety_level': 'High' if overall_risk < 0.3 else 
                           'Moderate' if overall_risk < 0.6 else 'Low',
            'monitoring_requirements': self._get_monitoring_requirements(delivery_params),
            'contraindications': self._get_contraindications(delivery_params)
        }
    
    def _get_monitoring_requirements(self, delivery_params: DeliveryParameters) -> List[str]:
        """Get monitoring requirements for delivery method."""
        requirements = ["Standard vital signs", "Injection site monitoring"]
        
        if delivery_params.vehicle == DeliveryVehicle.AAV:
            requirements.extend([
                "Anti-AAV antibody levels",
                "Liver function tests",
                "Complete blood count"
            ])
        
        if delivery_params.administration_route == "intravenous":
            requirements.append("Infusion reaction monitoring")
        
        if delivery_params.administration_route == "intrathecal":
            requirements.extend([
                "Neurological assessment",
                "CSF pressure monitoring",
                "MRI surveillance"
            ])
        
        return requirements
    
    def _get_contraindications(self, delivery_params: DeliveryParameters) -> List[str]:
        """Get contraindications for delivery method."""
        contraindications = []
        
        if delivery_params.vehicle == DeliveryVehicle.AAV:
            contraindications.extend([
                "Pre-existing anti-AAV antibodies",
                "Active hepatitis",
                "Severe immunodeficiency"
            ])
        
        if delivery_params.administration_route == "intrathecal":
            contraindications.extend([
                "Increased intracranial pressure",
                "CNS infection",
                "Coagulopathy"
            ])
        
        if delivery_params.dose_mg_kg > 5.0:
            contraindications.append("Pregnancy")
        
        return contraindications
