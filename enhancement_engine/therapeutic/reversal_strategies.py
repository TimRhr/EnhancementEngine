"""
Reversal strategies module for therapeutic CRISPR interventions.

This module provides comprehensive reversal and safety mechanisms including:
- Counter-editing strategies
- Emergency shutdown systems
- Temporal control mechanisms
- Risk mitigation protocols
- Long-term safety monitoring
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass
from enum import Enum
from datetime import datetime, timedelta

from ..models.therapeutic_data_classes import (
    ReversalStrategy, EmergencyProtocol, SafetyTrigger,
    ReversalEfficiency, TemporalControl
)


class ReversalType(Enum):
    """Types of reversal mechanisms."""
    COUNTER_EDITING = "counter_base_editing"
    GENE_REPLACEMENT = "replacement_therapy"
    SUICIDE_GENE = "suicide_gene_activation"
    IMMUNE_CLEARANCE = "immune_mediated_clearance"
    PHARMACOLOGICAL = "drug_mediated_reversal"
    NATURAL_DEGRADATION = "natural_turnover"


class TriggerCondition(Enum):
    """Conditions that trigger reversal."""
    ADVERSE_EVENT = "adverse_event_detected"
    PATIENT_REQUEST = "patient_withdrawal"
    EFFICACY_LOSS = "loss_of_efficacy"
    OFF_TARGET_EFFECTS = "off_target_detected"
    IMMUNE_REACTION = "immune_response"
    PROTOCOL_DEVIATION = "protocol_violation"


@dataclass
class ReversalPlan:
    """Comprehensive reversal plan for therapeutic intervention."""
    primary_strategy: ReversalStrategy
    backup_strategies: List[ReversalStrategy]
    trigger_conditions: List[TriggerCondition]
    implementation_timeline: Dict[str, int]  # Days to implement
    success_probability: float
    safety_considerations: List[str]
    monitoring_requirements: List[str]


class TherapeuticReversalSystem:
    """Main system for managing therapeutic reversals."""
    
    def __init__(self):
        """Initialize reversal system."""
        self.logger = logging.getLogger(__name__)
        self._load_reversal_protocols()
        self._initialize_emergency_systems()
    
    def _load_reversal_protocols(self) -> None:
        """Load reversal protocols for different therapeutic strategies."""
        self.reversal_protocols = {
            'base_editing': {
                'primary_method': ReversalType.COUNTER_EDITING,
                'feasibility': 0.8,
                'timeline_days': 14,
                'requirements': ['Original guide RNA', 'Counter base editor', 'Delivery system']
            },
            'gene_replacement': {
                'primary_method': ReversalType.GENE_REPLACEMENT,
                'feasibility': 0.6,
                'timeline_days': 30,
                'requirements': ['Wild-type gene copy', 'Excision system', 'New delivery']
            },
            'gene_silencing': {
                'primary_method': ReversalType.NATURAL_DEGRADATION,
                'feasibility': 0.9,
                'timeline_days': 7,
                'requirements': ['Discontinue treatment', 'RNA degradation pathway']
            },
            'knockout': {
                'primary_method': ReversalType.GENE_REPLACEMENT,
                'feasibility': 0.4,
                'timeline_days': 45,
                'requirements': ['Functional gene replacement', 'Targeted integration']
            }
        }
    
    def _initialize_emergency_systems(self) -> None:
        """Initialize emergency reversal systems."""
        self.emergency_protocols = {
            'immediate': {  # < 24 hours
                'methods': [ReversalType.PHARMACOLOGICAL, ReversalType.IMMUNE_CLEARANCE],
                'triggers': [TriggerCondition.ADVERSE_EVENT, TriggerCondition.IMMUNE_REACTION],
                'success_rate': 0.7
            },
            'short_term': {  # 1-7 days
                'methods': [ReversalType.SUICIDE_GENE, ReversalType.COUNTER_EDITING],
                'triggers': [TriggerCondition.OFF_TARGET_EFFECTS, TriggerCondition.PROTOCOL_DEVIATION],
                'success_rate': 0.8
            },
            'medium_term': {  # 1-4 weeks
                'methods': [ReversalType.GENE_REPLACEMENT, ReversalType.COUNTER_EDITING],
                'triggers': [TriggerCondition.EFFICACY_LOSS, TriggerCondition.PATIENT_REQUEST],
                'success_rate': 0.6
            }
        }
    
    def design_reversal_strategy(self, therapeutic_data: Dict[str, Any]) -> ReversalPlan:
        """
        Design comprehensive reversal strategy for therapeutic intervention.
        
        Args:
            therapeutic_data: Complete therapeutic intervention data
            
        Returns:
            Comprehensive reversal plan
        """
        strategy_type = therapeutic_data.get('strategy_type', 'base_editing')
        target_gene = therapeutic_data.get('target_gene', 'unknown')
        delivery_method = therapeutic_data.get('delivery_method', 'LNP')
        
        # Get base reversal protocol
        if strategy_type not in self.reversal_protocols:
            self.logger.warning(f"No reversal protocol for strategy: {strategy_type}")
            strategy_type = 'base_editing'  # Default fallback
        
        protocol = self.reversal_protocols[strategy_type]
        
        # Design primary reversal strategy
        primary_strategy = self._design_primary_reversal(
            therapeutic_data, protocol
        )
        
        # Design backup strategies
        backup_strategies = self._design_backup_strategies(
            therapeutic_data, primary_strategy
        )
        
        # Identify trigger conditions
        trigger_conditions = self._identify_trigger_conditions(therapeutic_data)
        
        # Create implementation timeline
        timeline = self._create_reversal_timeline(
            primary_strategy, backup_strategies
        )
        
        # Assess success probability
        success_prob = self._assess_reversal_success_probability(
            primary_strategy, backup_strategies, therapeutic_data
        )
        
        return ReversalPlan(
            primary_strategy=primary_strategy,
            backup_strategies=backup_strategies,
            trigger_conditions=trigger_conditions,
            implementation_timeline=timeline,
            success_probability=success_prob,
            safety_considerations=self._get_reversal_safety_considerations(
                primary_strategy, therapeutic_data
            ),
            monitoring_requirements=self._get_reversal_monitoring_requirements(
                primary_strategy, therapeutic_data
            )
        )
    
    def _design_primary_reversal(self, therapeutic_data: Dict[str, Any],
                               protocol: Dict[str, Any]) -> ReversalStrategy:
        """Design primary reversal strategy."""
        strategy_type = therapeutic_data.get('strategy_type', 'base_editing')
        target_gene = therapeutic_data.get('target_gene', 'unknown')
        
        if strategy_type == 'base_editing':
            return self._design_counter_editing_strategy(therapeutic_data)
        elif strategy_type == 'gene_replacement':
            return self._design_replacement_reversal_strategy(therapeutic_data)
        elif strategy_type == 'gene_silencing':
            return self._design_silencing_reversal_strategy(therapeutic_data)
        elif strategy_type == 'knockout':
            return self._design_knockout_reversal_strategy(therapeutic_data)
        else:
            return self._design_generic_reversal_strategy(therapeutic_data)
    
    def _design_counter_editing_strategy(self, therapeutic_data: Dict[str, Any]) -> ReversalStrategy:
        """Design counter base editing strategy."""
        original_edit = therapeutic_data.get('target_edit', 'C>T')
        target_position = therapeutic_data.get('target_position', 0)
        
        # Determine counter edit
        counter_edits = {
            'C>T': 'T>C',
            'A>G': 'G>A',
            'T>C': 'C>T',
            'G>A': 'A>G'
        }
        
        counter_edit = counter_edits.get(original_edit, 'unknown')
        
        # Select appropriate base editor
        if counter_edit in ['T>C', 'A>G']:
            editor_type = 'ABE'  # Adenine base editor
        else:
            editor_type = 'CBE'  # Cytosine base editor
        
        return ReversalStrategy(
            reversal_type=ReversalType.COUNTER_EDITING,
            mechanism="Counter base editing to revert original modification",
            components={
                'editor_type': editor_type,
                'counter_edit': counter_edit,
                'target_position': target_position,
                'guide_rna': f"Reverse complement of original guide",
                'delivery_method': therapeutic_data.get('delivery_method', 'LNP')
            },
            estimated_efficiency=0.75,
            implementation_time_days=14,
            safety_profile={
                'reversibility': True,
                'off_target_risk': 0.3,
                'immunogenicity': 0.2
            },
            requirements=[
                'Original guide RNA sequence',
                'Appropriate base editor',
                'Delivery vehicle',
                'Target cell access'
            ]
        )
    
    def _design_replacement_reversal_strategy(self, therapeutic_data: Dict[str, Any]) -> ReversalStrategy:
        """Design gene replacement reversal strategy."""
        target_gene = therapeutic_data.get('target_gene', 'unknown')
        
        return ReversalStrategy(
            reversal_type=ReversalType.GENE_REPLACEMENT,
            mechanism="Replace modified gene with wild-type version",
            components={
                'wild_type_sequence': f"Wild-type {target_gene} gene",
                'delivery_method': 'AAV',
                'integration_system': 'Homology-directed repair',
                'selection_marker': 'Optional safety marker'
            },
            estimated_efficiency=0.5,
            implementation_time_days=30,
            safety_profile={
                'reversibility': False,
                'off_target_risk': 0.4,
                'immunogenicity': 0.5
            },
            requirements=[
                'Wild-type gene sequence',
                'Efficient delivery system',
                'HDR enhancement',
                'Long-term monitoring'
            ]
        )
    
    def _design_silencing_reversal_strategy(self, therapeutic_data: Dict[str, Any]) -> ReversalStrategy:
        """Design gene silencing reversal strategy."""
        return ReversalStrategy(
            reversal_type=ReversalType.NATURAL_DEGRADATION,
            mechanism="Discontinue silencing treatment and allow natural RNA degradation",
            components={
                'discontinuation_protocol': 'Stop delivery of silencing agents',
                'degradation_pathway': 'Natural RNA turnover',
                'monitoring_markers': 'Gene expression recovery'
            },
            estimated_efficiency=0.9,
            implementation_time_days=7,
            safety_profile={
                'reversibility': True,
                'off_target_risk': 0.1,
                'immunogenicity': 0.1
            },
            requirements=[
                'Treatment discontinuation',
                'Expression monitoring',
                'Recovery assessment'
            ]
        )
    
    def _design_backup_strategies(self, therapeutic_data: Dict[str, Any],
                                primary_strategy: ReversalStrategy) -> List[ReversalStrategy]:
        """Design backup reversal strategies."""
        backup_strategies = []
        
        # Pharmacological reversal as first backup
        pharma_reversal = self._design_pharmacological_reversal(therapeutic_data)
        if pharma_reversal:
            backup_strategies.append(pharma_reversal)
        
        # Immune clearance as second backup
        immune_reversal = self._design_immune_clearance_strategy(therapeutic_data)
        if immune_reversal:
            backup_strategies.append(immune_reversal)
        
        # Suicide gene system if applicable
        if therapeutic_data.get('suicide_gene_included', False):
            suicide_reversal = self._design_suicide_gene_strategy(therapeutic_data)
            backup_strategies.append(suicide_reversal)
        
        return backup_strategies
    
    def _design_pharmacological_reversal(self, therapeutic_data: Dict[str, Any]) -> Optional[ReversalStrategy]:
        """Design pharmacological reversal strategy."""
        target_gene = therapeutic_data.get('target_gene', '')
        
        # Gene-specific pharmacological interventions
        pharma_interventions = {
            'PTPN22': {
                'drug': 'LYP inhibitor',
                'mechanism': 'Pharmacological inhibition of modified protein',
                'availability': 'Experimental'
            },
            'STAT4': {
                'drug': 'JAK inhibitor',
                'mechanism': 'Pathway inhibition downstream of STAT4',
                'availability': 'Approved (tofacitinib)'
            },
            'HLA-DRB1': {
                'drug': 'Immunosuppressants',
                'mechanism': 'General immune suppression',
                'availability': 'Approved'
            }
        }
        
        if target_gene in pharma_interventions:
            intervention = pharma_interventions[target_gene]
            
            return ReversalStrategy(
                reversal_type=ReversalType.PHARMACOLOGICAL,
                mechanism=intervention['mechanism'],
                components={
                    'drug': intervention['drug'],
                    'dosing': 'Standard therapeutic dosing',
                    'duration': 'Until reversal achieved',
                    'monitoring': 'Drug levels and efficacy markers'
                },
                estimated_efficiency=0.6,
                implementation_time_days=1,
                safety_profile={
                    'reversibility': True,
                    'off_target_risk': 0.3,
                    'immunogenicity': 0.1
                },
                requirements=[
                    'Drug availability',
                    'Patient tolerance',
                    'Dosing optimization'
                ]
            )
        
        return None
    
    def _design_immune_clearance_strategy(self, therapeutic_data: Dict[str, Any]) -> ReversalStrategy:
        """Design immune-mediated clearance strategy."""
        delivery_method = therapeutic_data.get('delivery_method', 'LNP')
        
        return ReversalStrategy(
            reversal_type=ReversalType.IMMUNE_CLEARANCE,
            mechanism="Enhance immune clearance of modified cells",
            components={
                'immunostimulation': 'Activate NK cells and CTLs',
                'targeting_antibodies': 'Anti-therapeutic antibodies',
                'complement_activation': 'Complement-mediated lysis',
                'delivery_method': delivery_method
            },
            estimated_efficiency=0.4,
            implementation_time_days=3,
            safety_profile={
                'reversibility': False,
                'off_target_risk': 0.6,
                'immunogenicity': 0.8
            },
            requirements=[
                'Immune system activation',
                'Targeting specificity',
                'Immune monitoring',
                'Safety management'
            ]
        )
    
    def _design_suicide_gene_strategy(self, therapeutic_data: Dict[str, Any]) -> ReversalStrategy:
        """Design suicide gene activation strategy."""
        return ReversalStrategy(
            reversal_type=ReversalType.SUICIDE_GENE,
            mechanism="Activate integrated suicide gene to eliminate modified cells",
            components={
                'suicide_gene': 'HSV-tk or iCasp9',
                'activating_drug': 'Ganciclovir or AP1903',
                'targeting_specificity': 'Modified cells only',
                'activation_kinetics': 'Rapid cell death'
            },
            estimated_efficiency=0.85,
            implementation_time_days=2,
            safety_profile={
                'reversibility': False,
                'off_target_risk': 0.2,
                'immunogenicity': 0.3
            },
            requirements=[
                'Suicide gene integration',
                'Activating drug availability',
                'Specific targeting',
                'Careful dosing'
            ]
        )
    
    def create_emergency_protocol(self, therapeutic_data: Dict[str, Any],
                                reversal_plan: ReversalPlan) -> EmergencyProtocol:
        """
        Create emergency reversal protocol.
        
        Args:
            therapeutic_data: Therapeutic intervention data
            reversal_plan: Planned reversal strategies
            
        Returns:
            Emergency protocol with immediate actions
        """
        # Define emergency triggers
        emergency_triggers = [
            SafetyTrigger(
                condition=TriggerCondition.ADVERSE_EVENT,
                severity_threshold='moderate',
                response_time_hours=6,
                automatic_activation=True
            ),
            SafetyTrigger(
                condition=TriggerCondition.IMMUNE_REACTION,
                severity_threshold='mild',
                response_time_hours=2,
                automatic_activation=True
            ),
            SafetyTrigger(
                condition=TriggerCondition.OFF_TARGET_EFFECTS,
                severity_threshold='any',
                response_time_hours=12,
                automatic_activation=False
            )
        ]
        
        # Create escalation plan
        escalation_steps = [
            {
                'level': 1,
                'timeframe': '0-6 hours',
                'actions': [
                    'Discontinue any ongoing treatment',
                    'Administer supportive care',
                    'Initiate monitoring protocols',
                    'Prepare reversal agents'
                ]
            },
            {
                'level': 2,
                'timeframe': '6-24 hours',
                'actions': [
                    'Implement primary reversal strategy',
                    'Administer immunosuppression if needed',
                    'Continuous patient monitoring',
                    'Prepare backup reversal methods'
                ]
            },
            {
                'level': 3,
                'timeframe': '24-72 hours',
                'actions': [
                    'Implement backup reversal strategies',
                    'Consider immune clearance methods',
                    'Specialized medical support',
                    'Long-term safety planning'
                ]
            }
        ]
        
        return EmergencyProtocol(
            protocol_id=f"EMERGENCY_{datetime.now().strftime('%Y%m%d_%H%M')}",
            triggers=emergency_triggers,
            escalation_steps=escalation_steps,
            contact_information={
                'primary_investigator': 'PI contact info',
                'medical_monitor': 'Medical monitor 24/7',
                'emergency_services': 'Local emergency services',
                'regulatory_authority': 'Regulatory reporting contact'
            },
            reversal_agents=self._prepare_reversal_agents(reversal_plan),
            monitoring_requirements=self._create_emergency_monitoring(),
            documentation_requirements=[
                'Time of event onset',
                'Severity assessment',
                'Actions taken',
                'Patient response',
                'Reversal efficiency'
            ]
        )
    
    def _prepare_reversal_agents(self, reversal_plan: ReversalPlan) -> Dict[str, Any]:
        """Prepare reversal agents for emergency use."""
        agents = {}
        
        # Primary reversal agents
        primary = reversal_plan.primary_strategy
        agents['primary'] = {
            'type': primary.reversal_type.value,
            'components': primary.components,
            'preparation_time': '2-4 hours',
            'storage_requirements': 'Ultra-low temperature',
            'expiration': '6 months'
        }
        
        # Backup agents
        for i, backup in enumerate(reversal_plan.backup_strategies):
            agents[f'backup_{i+1}'] = {
                'type': backup.reversal_type.value,
                'components': backup.components,
                'preparation_time': f"{backup.implementation_time_days} days",
                'storage_requirements': 'Standard conditions'
            }
        
        # Emergency drugs
        agents['emergency_drugs'] = {
            'immunosuppressants': ['Methylprednisolone', 'Cyclophosphamide'],
            'antihistamines': ['Diphenhydramine', 'Epinephrine'],
            'supportive_care': ['IV fluids', 'Oxygen', 'Monitoring equipment']
        }
        
        return agents
    
    def assess_reversal_efficiency(self, reversal_data: Dict[str, Any]) -> ReversalEfficiency:
        """
        Assess efficiency of reversal attempt.
        
        Args:
            reversal_data: Data from reversal implementation
            
        Returns:
            Efficiency assessment results
        """
        # Extract reversal metrics
        target_reduction = reversal_data.get('target_gene_expression_reduction', 0.0)
        functional_reversal = reversal_data.get('functional_parameter_reversal', 0.0)
        safety_events = reversal_data.get('safety_events', [])
        time_to_effect = reversal_data.get('time_to_effect_days', 0)
        
        # Calculate efficiency metrics
        molecular_efficiency = min(1.0, target_reduction)
        functional_efficiency = min(1.0, functional_reversal)
        
        # Safety efficiency (fewer events = higher efficiency)
        safety_efficiency = max(0.0, 1.0 - len(safety_events) * 0.2)
        
        # Temporal efficiency (faster = better)
        expected_time = reversal_data.get('expected_time_days', 14)
        temporal_efficiency = max(0.1, expected_time / max(time_to_effect, 1))
        
        # Overall efficiency
        overall_efficiency = np.mean([
            molecular_efficiency,
            functional_efficiency,
            safety_efficiency,
            min(1.0, temporal_efficiency)
        ])
        
        return ReversalEfficiency(
            overall_efficiency=overall_efficiency,
            molecular_efficiency=molecular_efficiency,
            functional_efficiency=functional_efficiency,
            safety_efficiency=safety_efficiency,
            temporal_efficiency=min(1.0, temporal_efficiency),
            success_criteria_met=overall_efficiency >= 0.7,
            recommendations=self._generate_efficiency_recommendations(
                molecular_efficiency, functional_efficiency, safety_efficiency
            )
        )
    
    def _generate_efficiency_recommendations(self, molecular: float,
                                           functional: float, safety: float) -> List[str]:
        """Generate recommendations based on efficiency metrics."""
        recommendations = []
        
        if molecular < 0.5:
            recommendations.append("Consider alternative reversal strategy")
            recommendations.append("Optimize delivery for better target engagement")
        
        if functional < 0.5:
            recommendations.append("Monitor functional recovery over longer timeframe")
            recommendations.append("Consider combination reversal approaches")
        
        if safety < 0.7:
            recommendations.append("Enhanced safety monitoring required")
            recommendations.append("Consider dose reduction for reversal agent")
        
        if all(metric >= 0.7 for metric in [molecular, functional, safety]):
            recommendations.append("Reversal successful - continue monitoring")
        
        return recommendations
    
    def generate_reversal_report(self, reversal_plan: ReversalPlan,
                               emergency_protocol: EmergencyProtocol) -> str:
        """Generate comprehensive reversal strategy report."""
        report = f"""
THERAPEUTIC REVERSAL STRATEGY REPORT
==================================

PRIMARY REVERSAL STRATEGY:
Type: {reversal_plan.primary_strategy.reversal_type.value}
Mechanism: {reversal_plan.primary_strategy.mechanism}
Estimated Efficiency: {reversal_plan.primary_strategy.estimated_efficiency:.1%}
Implementation Time: {reversal_plan.primary_strategy.implementation_time_days} days

BACKUP STRATEGIES:
{chr(10).join([f"- {strategy.reversal_type.value}: {strategy.estimated_efficiency:.1%} efficiency" 
               for strategy in reversal_plan.backup_strategies])}

TRIGGER CONDITIONS:
{chr(10).join([f"- {trigger.value}" for trigger in reversal_plan.trigger_conditions])}

SUCCESS PROBABILITY: {reversal_plan.success_probability:.1%}

EMERGENCY PROTOCOL:
Protocol ID: {emergency_protocol.protocol_id}
Emergency Triggers: {len(emergency_protocol.triggers)}
Escalation Levels: {len(emergency_protocol.escalation_steps)}

SAFETY CONSIDERATIONS:
{chr(10).join([f"- {consideration}" for consideration in reversal_plan.safety_considerations])}

MONITORING REQUIREMENTS:
{chr(10).join([f"- {requirement}" for requirement in reversal_plan.monitoring_requirements])}
        """
        
        return report.strip()
    
    # Additional helper methods...
    def _identify_trigger_conditions(self, therapeutic_data: Dict[str, Any]) -> List[TriggerCondition]:
        """Identify relevant trigger conditions for reversal."""
        triggers = [TriggerCondition.ADVERSE_EVENT, TriggerCondition.PATIENT_REQUEST]
        
        # Add strategy-specific triggers
        strategy_type = therapeutic_data.get('strategy_type', 'base_editing')
        if strategy_type in ['knockout', 'gene_replacement']:
            triggers.append(TriggerCondition.OFF_TARGET_EFFECTS)
        
        # Add delivery-specific triggers
        delivery_method = therapeutic_data.get('delivery_method', 'LNP')
        if delivery_method in ['AAV', 'lentivirus']:
            triggers.append(TriggerCondition.IMMUNE_REACTION)
        
        return triggers
