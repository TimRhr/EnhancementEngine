"""
Quality Control Framework for Therapeutic CRISPR Applications.

This module provides comprehensive quality control and validation for therapeutic
genetic interventions, including pre-clinical validation, manufacturing QC,
and clinical monitoring protocols.
"""

import logging
import re
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any, Set
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from datetime import datetime, timedelta

from ..models.data_classes import GuideRNA
from ..models.therapeutic_data_classes import (
    CorrectionStrategy,
    TherapeuticTarget,
    TherapeuticEfficacy,
)
from ..models.disease_constants import QUALITY_CONTROL_SPECS


class QCStatus(Enum):
    """Quality control status levels."""

    PASS = "pass"
    FAIL = "fail"
    WARNING = "warning"
    PENDING = "pending"
    NOT_TESTED = "not_tested"


class QCCategory(Enum):
    """Quality control categories."""

    DESIGN_VALIDATION = "design_validation"
    IN_VITRO_TESTING = "in_vitro_testing"
    MANUFACTURING = "manufacturing"
    RELEASE_TESTING = "release_testing"
    CLINICAL_MONITORING = "clinical_monitoring"
    POST_MARKET_SURVEILLANCE = "post_market_surveillance"


@dataclass
class QCTest:
    """Individual quality control test."""

    test_name: str
    test_category: QCCategory
    test_method: str
    acceptance_criteria: Dict[str, Any]
    test_result: Optional[Any] = None
    status: QCStatus = QCStatus.NOT_TESTED
    test_date: Optional[datetime] = None
    operator: Optional[str] = None
    comments: str = ""

    @property
    def is_critical(self) -> bool:
        """Check if test is critical for release."""
        critical_tests = {
            "guide_specificity",
            "editing_efficiency",
            "viability",
            "sterility",
            "endotoxin",
            "identity",
        }
        return any(critical in self.test_name.lower() for critical in critical_tests)


@dataclass
class QCBatch:
    """Quality control batch for a therapeutic product."""

    batch_id: str
    product_name: str
    manufacturing_date: datetime
    tests: List[QCTest]
    overall_status: QCStatus = QCStatus.PENDING
    release_date: Optional[datetime] = None
    expiry_date: Optional[datetime] = None

    @property
    def critical_tests_passed(self) -> bool:
        """Check if all critical tests have passed."""
        critical_tests = [test for test in self.tests if test.is_critical]
        return all(test.status == QCStatus.PASS for test in critical_tests)

    @property
    def release_eligible(self) -> bool:
        """Check if batch is eligible for release."""
        return self.critical_tests_passed and self.overall_status == QCStatus.PASS


@dataclass
class ValidationReport:
    """Comprehensive validation report."""

    validation_type: str
    target_product: str
    validation_date: datetime
    test_results: List[QCTest]
    overall_conclusion: str
    recommendations: List[str]
    next_actions: List[str]
    reviewer: str
    approval_status: str = "pending"


class DesignValidationEngine:
    """Validates therapeutic CRISPR designs before implementation."""

    def __init__(self):
        """Initialize design validation engine."""
        self.logger = logging.getLogger(__name__)

    def validate_guide_design(
        self, guide_rna: GuideRNA, target_sequence: str
    ) -> List[QCTest]:
        """
        Validate guide RNA design quality.

        Args:
            guide_rna: Guide RNA to validate
            target_sequence: Target sequence context

        Returns:
            List of QC test results
        """
        tests = []

        # GC content validation
        gc_test = self._validate_gc_content(guide_rna)
        tests.append(gc_test)

        # Specificity validation
        specificity_test = self._validate_specificity(guide_rna)
        tests.append(specificity_test)

        # Efficiency prediction validation
        efficiency_test = self._validate_efficiency_prediction(guide_rna)
        tests.append(efficiency_test)

        # Off-target validation
        off_target_test = self._validate_off_targets(guide_rna)
        tests.append(off_target_test)

        # Secondary structure validation
        structure_test = self._validate_secondary_structure(guide_rna)
        tests.append(structure_test)

        return tests

    def _validate_gc_content(self, guide_rna: GuideRNA) -> QCTest:
        """Validate GC content of guide RNA."""
        gc_content = guide_rna.gc_content

        # Acceptance criteria: 30-70% GC content
        criteria = {"min_gc": 30, "max_gc": 70, "optimal_range": (40, 60)}

        if criteria["min_gc"] <= gc_content <= criteria["max_gc"]:
            if (
                criteria["optimal_range"][0]
                <= gc_content
                <= criteria["optimal_range"][1]
            ):
                status = QCStatus.PASS
            else:
                status = QCStatus.WARNING
        else:
            status = QCStatus.FAIL

        return QCTest(
            test_name="GC_Content_Validation",
            test_category=QCCategory.DESIGN_VALIDATION,
            test_method="sequence_analysis",
            acceptance_criteria=criteria,
            test_result={"gc_content": gc_content},
            status=status,
            test_date=datetime.now(),
            comments=f"GC content: {gc_content:.1f}%",
        )

    def _validate_specificity(self, guide_rna: GuideRNA) -> QCTest:
        """Validate guide RNA specificity."""
        specificity_score = guide_rna.specificity_score
        off_target_count = len(guide_rna.off_targets)
        high_risk_off_targets = len(guide_rna.high_risk_off_targets)

        criteria = {
            "min_specificity": 0.8,
            "max_off_targets": 5,
            "max_high_risk_off_targets": 1,
        }

        if (
            specificity_score >= criteria["min_specificity"]
            and off_target_count <= criteria["max_off_targets"]
            and high_risk_off_targets <= criteria["max_high_risk_off_targets"]
        ):
            status = QCStatus.PASS
        elif specificity_score >= 0.6 and high_risk_off_targets <= 2:
            status = QCStatus.WARNING
        else:
            status = QCStatus.FAIL

        return QCTest(
            test_name="Specificity_Validation",
            test_category=QCCategory.DESIGN_VALIDATION,
            test_method="bioinformatics_analysis",
            acceptance_criteria=criteria,
            test_result={
                "specificity_score": specificity_score,
                "off_target_count": off_target_count,
                "high_risk_off_targets": high_risk_off_targets,
            },
            status=status,
            test_date=datetime.now(),
            comments=f"Specificity: {specificity_score:.3f}, Off-targets: {off_target_count}",
        )

    def _validate_efficiency_prediction(self, guide_rna: GuideRNA) -> QCTest:
        """Validate efficiency prediction for guide RNA."""
        efficiency = guide_rna.efficiency_score.overall_efficiency

        criteria = {"min_efficiency": 0.5, "optimal_efficiency": 0.7}

        if efficiency >= criteria["optimal_efficiency"]:
            status = QCStatus.PASS
        elif efficiency >= criteria["min_efficiency"]:
            status = QCStatus.WARNING
        else:
            status = QCStatus.FAIL

        return QCTest(
            test_name="Efficiency_Prediction",
            test_category=QCCategory.DESIGN_VALIDATION,
            test_method="machine_learning_prediction",
            acceptance_criteria=criteria,
            test_result={"predicted_efficiency": efficiency},
            status=status,
            test_date=datetime.now(),
            comments=f"Predicted efficiency: {efficiency:.3f}",
        )

    def _validate_off_targets(self, guide_rna: GuideRNA) -> QCTest:
        """Validate off-target analysis quality."""
        off_targets = guide_rna.off_targets

        # Check for off-targets in essential genes
        essential_hits = sum(1 for ot in off_targets if ot.essential_gene)

        # Check for high-confidence off-targets
        high_conf_hits = sum(1 for ot in off_targets if ot.score > 0.5)

        criteria = {
            "max_essential_hits": 0,
            "max_high_confidence_hits": 2,
            "analysis_depth": "genome_wide",
        }

        if (
            essential_hits <= criteria["max_essential_hits"]
            and high_conf_hits <= criteria["max_high_confidence_hits"]
        ):
            status = QCStatus.PASS
        elif essential_hits <= 1 and high_conf_hits <= 4:
            status = QCStatus.WARNING
        else:
            status = QCStatus.FAIL

        return QCTest(
            test_name="Off_Target_Analysis",
            test_category=QCCategory.DESIGN_VALIDATION,
            test_method="computational_prediction",
            acceptance_criteria=criteria,
            test_result={
                "essential_gene_hits": essential_hits,
                "high_confidence_hits": high_conf_hits,
                "total_off_targets": len(off_targets),
            },
            status=status,
            test_date=datetime.now(),
            comments=f"Essential gene hits: {essential_hits}, High-conf hits: {high_conf_hits}",
        )

    def _validate_secondary_structure(self, guide_rna: GuideRNA) -> QCTest:
        """Validate guide RNA secondary structure."""
        # Simplified secondary structure validation
        sequence = guide_rna.sequence

        # Check for problematic sequences
        poly_stretches = max(
            len(match.group())
            for match in [re.search(f"{base}+", sequence) for base in "ATCG"]
            if match
        )

        hairpin_potential = self._calculate_hairpin_potential(sequence)

        criteria = {"max_poly_stretch": 4, "max_hairpin_score": 0.3}

        if (
            poly_stretches <= criteria["max_poly_stretch"]
            and hairpin_potential <= criteria["max_hairpin_score"]
        ):
            status = QCStatus.PASS
        elif poly_stretches <= 6 and hairpin_potential <= 0.5:
            status = QCStatus.WARNING
        else:
            status = QCStatus.FAIL

        return QCTest(
            test_name="Secondary_Structure",
            test_category=QCCategory.DESIGN_VALIDATION,
            test_method="structure_prediction",
            acceptance_criteria=criteria,
            test_result={
                "max_poly_stretch": poly_stretches,
                "hairpin_potential": hairpin_potential,
            },
            status=status,
            test_date=datetime.now(),
            comments=f"Max poly stretch: {poly_stretches}, Hairpin score: {hairpin_potential:.3f}",
        )

    def _calculate_hairpin_potential(self, sequence: str) -> float:
        """Calculate potential for hairpin formation (simplified)."""
        # Simplified hairpin potential calculation
        import re

        # Look for inverted repeats
        hairpin_score = 0.0
        for i in range(len(sequence) - 6):
            for j in range(i + 6, len(sequence)):
                substr1 = sequence[i : i + 3]
                substr2 = sequence[j : j + 3]

                # Check if reverse complement
                complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
                rev_comp = "".join(
                    complement_map.get(base, base) for base in reversed(substr2)
                )

                if substr1 == rev_comp:
                    hairpin_score += 0.1

        return min(1.0, hairpin_score)


class InVitroTestingEngine:
    """Manages in vitro testing protocols for therapeutic validation."""

    def __init__(self):
        """Initialize in vitro testing engine."""
        self.logger = logging.getLogger(__name__)

    def design_in_vitro_validation(
        self,
        correction_strategy: CorrectionStrategy,
        therapeutic_target: TherapeuticTarget,
    ) -> List[QCTest]:
        """
        Design in vitro validation protocol.

        Args:
            correction_strategy: Correction strategy to validate
            therapeutic_target: Therapeutic target information

        Returns:
            List of required in vitro tests
        """
        tests = []

        # Basic functionality tests
        tests.extend(self._design_functionality_tests(correction_strategy))

        # Safety tests
        tests.extend(self._design_safety_tests(correction_strategy))

        # Efficacy tests
        tests.extend(
            self._design_efficacy_tests(correction_strategy, therapeutic_target)
        )

        # Specificity tests
        tests.extend(self._design_specificity_tests(correction_strategy))

        return tests

    def _design_functionality_tests(
        self, correction_strategy: CorrectionStrategy
    ) -> List[QCTest]:
        """Design functionality validation tests."""
        tests = []

        # Editing efficiency test
        tests.append(
            QCTest(
                test_name="Editing_Efficiency_In_Vitro",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="targeted_sequencing",
                acceptance_criteria={
                    "min_efficiency": 0.3,
                    "target_efficiency": 0.7,
                    "measurement_method": "NGS",
                    "cell_lines": ["HEK293T", "disease_relevant_cells"],
                },
                comments="Measure on-target editing efficiency in relevant cell models",
            )
        )

        # Viability test
        tests.append(
            QCTest(
                test_name="Cell_Viability",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="MTT_assay",
                acceptance_criteria={
                    "min_viability": 0.8,
                    "time_points": [24, 48, 72],
                    "controls": ["untreated", "mock_transfected"],
                },
                comments="Assess cell viability after treatment",
            )
        )

        # Delivery efficiency test
        tests.append(
            QCTest(
                test_name="Delivery_Efficiency",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="flow_cytometry",
                acceptance_criteria={
                    "min_delivery": 0.5,
                    "reporter": "fluorescent_protein",
                    "analysis_method": "FACS",
                },
                comments="Measure delivery system efficiency",
            )
        )

        return tests

    def _design_safety_tests(
        self, correction_strategy: CorrectionStrategy
    ) -> List[QCTest]:
        """Design safety validation tests."""
        tests = []

        # Off-target validation
        tests.append(
            QCTest(
                test_name="Off_Target_Validation",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="GUIDE_seq",
                acceptance_criteria={
                    "max_off_targets": 5,
                    "detection_threshold": 0.01,
                    "validation_sites": "top_10_predicted",
                },
                comments="Validate predicted off-target sites experimentally",
            )
        )

        # Genotoxicity assessment
        tests.append(
            QCTest(
                test_name="Genotoxicity_Assessment",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="comet_assay",
                acceptance_criteria={
                    "max_dna_damage": 1.5,  # Fold increase over control
                    "recovery_time": 24,  # Hours
                    "dose_response": True,
                },
                comments="Assess DNA damage and repair capacity",
            )
        )

        # Chromosomal stability
        tests.append(
            QCTest(
                test_name="Chromosomal_Stability",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="karyotype_analysis",
                acceptance_criteria={
                    "max_aberrations": 0.05,  # 5% of cells
                    "analysis_time": "2_weeks_post_treatment",
                    "cell_passages": 10,
                },
                comments="Long-term chromosomal stability assessment",
            )
        )

        return tests

    def _design_efficacy_tests(
        self,
        correction_strategy: CorrectionStrategy,
        therapeutic_target: TherapeuticTarget,
    ) -> List[QCTest]:
        """Design efficacy validation tests."""
        tests = []

        # Functional correction test
        tests.append(
            QCTest(
                test_name="Functional_Correction",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="functional_assay",
                acceptance_criteria={
                    "min_functional_improvement": 0.5,
                    "disease_relevant_readout": True,
                    "quantitative_measurement": True,
                },
                comments=f"Measure functional correction relevant to {therapeutic_target.disease_association}",
            )
        )

        # Protein expression test
        tests.append(
            QCTest(
                test_name="Protein_Expression",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="western_blot",
                acceptance_criteria={
                    "expected_expression_change": "gene_dependent",
                    "quantification_method": "densitometry",
                    "normalization": "loading_controls",
                },
                comments="Verify expected protein expression changes",
            )
        )

        # Cellular phenotype test
        tests.append(
            QCTest(
                test_name="Cellular_Phenotype",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="phenotypic_assay",
                acceptance_criteria={
                    "phenotype_rescue": True,
                    "quantitative_metrics": True,
                    "statistical_significance": 0.05,
                },
                comments="Assess rescue of disease-relevant cellular phenotypes",
            )
        )

        return tests

    def _design_specificity_tests(
        self, correction_strategy: CorrectionStrategy
    ) -> List[QCTest]:
        """Design specificity validation tests."""
        tests = []

        # On-target specificity
        tests.append(
            QCTest(
                test_name="On_Target_Specificity",
                test_category=QCCategory.IN_VITRO_TESTING,
                test_method="deep_sequencing",
                acceptance_criteria={
                    "target_site_editing": ">90%",
                    "intended_edit_percentage": ">80%",
                    "unintended_indels": "<10%",
                },
                comments="Verify specificity at intended target site",
            )
        )

        # Bystander editing assessment (for base editors)
        if correction_strategy.editing_method == "base_editing":
            tests.append(
                QCTest(
                    test_name="Bystander_Editing",
                    test_category=QCCategory.IN_VITRO_TESTING,
                    test_method="amplicon_sequencing",
                    acceptance_criteria={
                        "max_bystander_edits": 1,
                        "editing_window_analysis": True,
                        "base_editing_specificity": ">90%",
                    },
                    comments="Assess unintended bystander base edits",
                )
            )

        return tests


class ManufacturingQCEngine:
    """Manages manufacturing quality control for therapeutic products."""

    def __init__(self):
        """Initialize manufacturing QC engine."""
        self.logger = logging.getLogger(__name__)

    def create_manufacturing_batch(
        self, product_name: str, manufacturing_date: datetime
    ) -> QCBatch:
        """
        Create new manufacturing batch with required QC tests.

        Args:
            product_name: Name of therapeutic product
            manufacturing_date: Date of manufacturing

        Returns:
            QC batch with required tests
        """
        batch_id = f"{product_name}_{manufacturing_date.strftime('%Y%m%d')}_{np.random.randint(1000, 9999)}"

        # Define required manufacturing tests
        required_tests = self._define_manufacturing_tests(product_name)

        batch = QCBatch(
            batch_id=batch_id,
            product_name=product_name,
            manufacturing_date=manufacturing_date,
            tests=required_tests,
            expiry_date=manufacturing_date + timedelta(days=365),  # 1 year shelf life
        )

        return batch

    def _define_manufacturing_tests(self, product_name: str) -> List[QCTest]:
        """Define required manufacturing QC tests."""
        tests = []

        # Identity test
        tests.append(
            QCTest(
                test_name="Identity_Verification",
                test_category=QCCategory.MANUFACTURING,
                test_method="sequencing",
                acceptance_criteria={
                    "sequence_match": "100%",
                    "method": "sanger_sequencing",
                    "coverage": "full_construct",
                },
                comments="Verify product identity by sequencing",
            )
        )

        # Purity test
        tests.append(
            QCTest(
                test_name="Purity_Assessment",
                test_category=QCCategory.MANUFACTURING,
                test_method="HPLC",
                acceptance_criteria={
                    "min_purity": 0.95,
                    "max_impurities": 0.05,
                    "method_validation": True,
                },
                comments="Assess product purity",
            )
        )

        # Concentration test
        tests.append(
            QCTest(
                test_name="Concentration_Determination",
                test_category=QCCategory.MANUFACTURING,
                test_method="UV_spectroscopy",
                acceptance_criteria={
                    "target_concentration": "label_claim",
                    "tolerance": "±10%",
                    "method": "A260_measurement",
                },
                comments="Determine product concentration",
            )
        )

        # Sterility test
        tests.append(
            QCTest(
                test_name="Sterility_Test",
                test_category=QCCategory.MANUFACTURING,
                test_method="USP_71",
                acceptance_criteria={
                    "microbial_growth": "negative",
                    "incubation_time": "14_days",
                    "test_volume": "sample_dependent",
                },
                comments="Sterility testing per USP <71>",
            )
        )

        # Endotoxin test
        tests.append(
            QCTest(
                test_name="Endotoxin_Test",
                test_category=QCCategory.MANUFACTURING,
                test_method="LAL_assay",
                acceptance_criteria={
                    "max_endotoxin": "5.0_EU_per_kg",
                    "method": "kinetic_chromogenic",
                    "sensitivity": "0.1_EU_per_mL",
                },
                comments="Bacterial endotoxin testing",
            )
        )

        # pH test
        tests.append(
            QCTest(
                test_name="pH_Measurement",
                test_category=QCCategory.MANUFACTURING,
                test_method="pH_meter",
                acceptance_criteria={
                    "target_pH": "7.4",
                    "tolerance": "±0.2",
                    "temperature": "25°C",
                },
                comments="pH measurement of final product",
            )
        )

        # Osmolality test
        tests.append(
            QCTest(
                test_name="Osmolality_Test",
                test_category=QCCategory.MANUFACTURING,
                test_method="freezing_point_depression",
                acceptance_criteria={
                    "target_osmolality": "280-320_mOsm_per_kg",
                    "method_precision": "±5%",
                },
                comments="Osmolality measurement",
            )
        )

        return tests

    def evaluate_batch_release(self, batch: QCBatch) -> bool:
        """
        Evaluate if batch meets release criteria.

        Args:
            batch: QC batch to evaluate

        Returns:
            True if batch can be released
        """
        # Check if all critical tests are complete and passed
        if not batch.critical_tests_passed:
            self.logger.warning(f"Batch {batch.batch_id}: Critical tests not passed")
            return False

        # Check for any failed tests
        failed_tests = [test for test in batch.tests if test.status == QCStatus.FAIL]
        if failed_tests:
            self.logger.warning(
                f"Batch {batch.batch_id}: {len(failed_tests)} tests failed"
            )
            return False

        # Check for warnings in critical tests
        critical_warnings = [
            test
            for test in batch.tests
            if test.is_critical and test.status == QCStatus.WARNING
        ]
        if critical_warnings:
            self.logger.warning(f"Batch {batch.batch_id}: Critical tests have warnings")
            return False

        # Batch can be released
        batch.overall_status = QCStatus.PASS
        batch.release_date = datetime.now()

        self.logger.info(f"Batch {batch.batch_id}: Approved for release")
        return True


class ClinicalMonitoringEngine:
    """Manages clinical monitoring and post-market surveillance."""

    def __init__(self):
        """Initialize clinical monitoring engine."""
        self.logger = logging.getLogger(__name__)

    def design_clinical_monitoring_protocol(
        self,
        therapeutic_target: TherapeuticTarget,
        predicted_efficacy: TherapeuticEfficacy,
    ) -> List[QCTest]:
        """
        Design clinical monitoring protocol.

        Args:
            therapeutic_target: Therapeutic target information
            predicted_efficacy: Predicted efficacy

        Returns:
            List of clinical monitoring tests
        """
        tests = []

        # Safety monitoring
        tests.extend(self._design_safety_monitoring(therapeutic_target))

        # Efficacy monitoring
        tests.extend(
            self._design_efficacy_monitoring(therapeutic_target, predicted_efficacy)
        )

        # Biomarker monitoring
        tests.extend(self._design_biomarker_monitoring(therapeutic_target))

        # Long-term surveillance
        tests.extend(self._design_long_term_surveillance(therapeutic_target))

        return tests

    def _design_safety_monitoring(
        self, therapeutic_target: TherapeuticTarget
    ) -> List[QCTest]:
        """Design safety monitoring tests."""
        tests = []

        # Vital signs monitoring
        tests.append(
            QCTest(
                test_name="Vital_Signs_Monitoring",
                test_category=QCCategory.CLINICAL_MONITORING,
                test_method="standard_clinical_assessment",
                acceptance_criteria={
                    "frequency": "continuous_24h_then_q4h_x24h",
                    "parameters": ["BP", "HR", "RR", "Temp", "O2_sat"],
                    "alert_thresholds": "age_appropriate",
                },
                comments="Continuous monitoring for acute reactions",
            )
        )

        # Laboratory safety monitoring
        tests.append(
            QCTest(
                test_name="Laboratory_Safety_Panel",
                test_category=QCCategory.CLINICAL_MONITORING,
                test_method="clinical_laboratory",
                acceptance_criteria={
                    "frequency": "baseline_24h_1wk_1mo_3mo_6mo",
                    "tests": ["CBC", "CMP", "LFTs", "coag_studies"],
                    "grade_3_4_toxicity": "dose_limiting",
                },
                comments="Comprehensive laboratory safety assessment",
            )
        )

        # Immune response monitoring
        tests.append(
            QCTest(
                test_name="Immune_Response_Monitoring",
                test_category=QCCategory.CLINICAL_MONITORING,
                test_method="immunoassay",
                acceptance_criteria={
                    "anti_cas_antibodies": "baseline_1mo_3mo_6mo",
                    "T_cell_response": "baseline_1mo_6mo",
                    "inflammatory_markers": "baseline_24h_1wk_1mo",
                },
                comments="Monitor for immune reactions to CRISPR components",
            )
        )

        return tests

    def _design_efficacy_monitoring(
        self,
        therapeutic_target: TherapeuticTarget,
        predicted_efficacy: TherapeuticEfficacy,
    ) -> List[QCTest]:
        """Design efficacy monitoring tests."""
        tests = []

        # Disease-specific efficacy measures
        if therapeutic_target.disease_association == "rheumatoid_arthritis":
            tests.append(
                QCTest(
                    test_name="RA_Disease_Activity",
                    test_category=QCCategory.CLINICAL_MONITORING,
                    test_method="clinical_assessment",
                    acceptance_criteria={
                        "measures": ["DAS28", "CDAI", "SDAI", "HAQ"],
                        "frequency": "baseline_2wk_1mo_3mo_6mo_12mo",
                        "primary_endpoint": "DAS28_remission_6mo",
                    },
                    comments="Standard RA disease activity measures",
                )
            )

        # Functional assessment
        tests.append(
            QCTest(
                test_name="Functional_Assessment",
                test_category=QCCategory.CLINICAL_MONITORING,
                test_method="standardized_questionnaire",
                acceptance_criteria={
                    "instruments": ["disease_specific", "generic_QoL"],
                    "frequency": "baseline_1mo_3mo_6mo_12mo",
                    "validated_instruments": True,
                },
                comments="Patient-reported functional outcomes",
            )
        )

        return tests

    def _design_biomarker_monitoring(
        self, therapeutic_target: TherapeuticTarget
    ) -> List[QCTest]:
        """Design biomarker monitoring tests."""
        tests = []

        # Editing confirmation
        tests.append(
            QCTest(
                test_name="Editing_Confirmation",
                test_category=QCCategory.CLINICAL_MONITORING,
                test_method="NGS",
                acceptance_criteria={
                    "sample_types": ["blood", "target_tissue_if_accessible"],
                    "frequency": "1wk_1mo_3mo_6mo",
                    "detection_limit": "1%",
                    "quantitative": True,
                },
                comments="Confirm and quantify genetic editing",
            )
        )

        # Off-target monitoring
        tests.append(
            QCTest(
                test_name="Off_Target_Monitoring",
                test_category=QCCategory.CLINICAL_MONITORING,
                test_method="targeted_sequencing",
                acceptance_criteria={
                    "sites_monitored": "top_10_predicted",
                    "frequency": "1mo_6mo_12mo",
                    "detection_threshold": "0.1%",
                },
                comments="Monitor predicted off-target sites",
            )
        )

        # Protein biomarkers
        tests.append(
            QCTest(
                test_name="Protein_Biomarkers",
                test_category=QCCategory.CLINICAL_MONITORING,
                test_method="immunoassay",
                acceptance_criteria={
                    "biomarkers": "disease_and_target_specific",
                    "frequency": "baseline_1mo_3mo_6mo",
                    "quantitative_measurement": True,
                },
                comments="Disease and intervention-specific protein biomarkers",
            )
        )

        return tests

    def _design_long_term_surveillance(
        self, therapeutic_target: TherapeuticTarget
    ) -> List[QCTest]:
        """Design long-term surveillance tests."""
        tests = []

        # Cancer surveillance
        tests.append(
            QCTest(
                test_name="Cancer_Surveillance",
                test_category=QCCategory.POST_MARKET_SURVEILLANCE,
                test_method="clinical_examination_imaging",
                acceptance_criteria={
                    "frequency": "annually_for_15_years",
                    "screening": "age_appropriate_cancer_screening",
                    "enhanced_surveillance": "if_indicated",
                },
                comments="Long-term cancer risk monitoring",
            )
        )

        # Reproductive health (if applicable)
        tests.append(
            QCTest(
                test_name="Reproductive_Health_Monitoring",
                test_category=QCCategory.POST_MARKET_SURVEILLANCE,
                test_method="clinical_assessment",
                acceptance_criteria={
                    "frequency": "if_pregnancy_desired",
                    "genetic_counseling": "mandatory",
                    "partner_testing": "recommended",
                },
                comments="Reproductive health and genetic counseling",
            )
        )

        # Neurological assessment (for systemic interventions)
        tests.append(
            QCTest(
                test_name="Neurological_Assessment",
                test_category=QCCategory.POST_MARKET_SURVEILLANCE,
                test_method="clinical_neurological_exam",
                acceptance_criteria={
                    "frequency": "annually_for_5_years",
                    "cognitive_testing": "if_indicated",
                    "detailed_assessment": "if_symptoms",
                },
                comments="Monitor for delayed neurological effects",
            )
        )

        return tests


class QualityControlOrchestrator:
    """Main quality control orchestrator for therapeutic applications."""

    def __init__(self):
        """Initialize quality control orchestrator."""
        self.design_validator = DesignValidationEngine()
        self.in_vitro_tester = InVitroTestingEngine()
        self.manufacturing_qc = ManufacturingQCEngine()
        self.clinical_monitor = ClinicalMonitoringEngine()
        self.logger = logging.getLogger(__name__)

    def create_comprehensive_qc_plan(
        self,
        correction_strategy: CorrectionStrategy,
        therapeutic_target: TherapeuticTarget,
        predicted_efficacy: TherapeuticEfficacy,
    ) -> ValidationReport:
        """
        Create comprehensive quality control plan.

        Args:
            correction_strategy: Correction strategy to validate
            therapeutic_target: Therapeutic target
            predicted_efficacy: Predicted efficacy

        Returns:
            Comprehensive validation report
        """
        try:
            all_tests = []

            # Design validation
            if correction_strategy.guide_rnas:
                design_tests = self.design_validator.validate_guide_design(
                    correction_strategy.guide_rnas[0], "target_sequence"
                )
                all_tests.extend(design_tests)

            # In vitro validation
            in_vitro_tests = self.in_vitro_tester.design_in_vitro_validation(
                correction_strategy, therapeutic_target
            )
            all_tests.extend(in_vitro_tests)

            # Manufacturing QC
            manufacturing_batch = self.manufacturing_qc.create_manufacturing_batch(
                f"{therapeutic_target.gene_symbol}_therapeutic", datetime.now()
            )
            all_tests.extend(manufacturing_batch.tests)

            # Clinical monitoring
            clinical_tests = self.clinical_monitor.design_clinical_monitoring_protocol(
                therapeutic_target, predicted_efficacy
            )
            all_tests.extend(clinical_tests)

            # Generate recommendations
            recommendations = self._generate_qc_recommendations(
                all_tests, therapeutic_target
            )

            # Generate next actions
            next_actions = self._generate_next_actions(all_tests)

            validation_report = ValidationReport(
                validation_type="comprehensive_therapeutic_validation",
                target_product=f"{therapeutic_target.gene_symbol}_{therapeutic_target.target_variant}_therapy",
                validation_date=datetime.now(),
                test_results=all_tests,
                overall_conclusion=self._generate_overall_conclusion(all_tests),
                recommendations=recommendations,
                next_actions=next_actions,
                reviewer="QC_System",
            )

            self.logger.info(
                f"Created comprehensive QC plan with {len(all_tests)} tests"
            )
            return validation_report

        except Exception as e:
            self.logger.error(f"QC plan creation failed: {e}")
            raise Exception(f"Failed to create QC plan: {e}")

    def _generate_qc_recommendations(
        self, tests: List[QCTest], therapeutic_target: TherapeuticTarget
    ) -> List[str]:
        """Generate QC recommendations based on test requirements."""
        recommendations = []

        # Count tests by category
        category_counts = defaultdict(int)
        for test in tests:
            category_counts[test.test_category] += 1

        recommendations.append(
            f"Comprehensive validation plan includes {len(tests)} tests across {len(category_counts)} categories"
        )

        # Category-specific recommendations
        if category_counts[QCCategory.DESIGN_VALIDATION] > 0:
            recommendations.append(
                "Complete design validation before proceeding to in vitro testing"
            )

        if category_counts[QCCategory.IN_VITRO_TESTING] > 0:
            recommendations.append(
                "Establish validated in vitro models representative of target disease"
            )

        if category_counts[QCCategory.MANUFACTURING] > 0:
            recommendations.append("Implement GMP-compliant manufacturing processes")

        if category_counts[QCCategory.CLINICAL_MONITORING] > 0:
            recommendations.append(
                "Establish robust clinical monitoring infrastructure before patient treatment"
            )

        # Target-specific recommendations
        if therapeutic_target.target_tissue == "hematopoietic_system":
            recommendations.append(
                "Special attention to hematopoietic toxicity monitoring required"
            )

        if therapeutic_target.therapeutic_approach.value == "replacement":
            recommendations.append(
                "Enhanced identity and purity testing for replacement constructs"
            )

        return recommendations

    def _generate_next_actions(self, tests: List[QCTest]) -> List[str]:
        """Generate next actions based on test plan."""
        next_actions = []

        # Immediate actions
        not_tested = [test for test in tests if test.status == QCStatus.NOT_TESTED]
        if not_tested:
            next_actions.append(f"Execute {len(not_tested)} pending validation tests")

        # Critical path actions
        design_tests = [
            test for test in tests if test.test_category == QCCategory.DESIGN_VALIDATION
        ]
        if design_tests:
            next_actions.append(
                "Complete design validation as prerequisite for subsequent testing"
            )

        # Regulatory actions
        next_actions.extend(
            [
                "Submit validation plan to regulatory authorities if required",
                "Establish quality agreements with testing laboratories",
                "Train personnel on QC procedures and acceptance criteria",
                "Implement electronic QC data management system",
            ]
        )

        return next_actions

    def _generate_overall_conclusion(self, tests: List[QCTest]) -> str:
        """Generate overall conclusion from test plan."""
        total_tests = len(tests)

        # Count by category
        category_counts = defaultdict(int)
        for test in tests:
            category_counts[test.test_category] += 1

        critical_tests = sum(1 for test in tests if test.is_critical)

        conclusion = f"""
Comprehensive quality control plan established for therapeutic genetic intervention.

Total Tests: {total_tests}
Critical Tests: {critical_tests}

Test Distribution:
- Design Validation: {category_counts[QCCategory.DESIGN_VALIDATION]}
- In Vitro Testing: {category_counts[QCCategory.IN_VITRO_TESTING]}
- Manufacturing QC: {category_counts[QCCategory.MANUFACTURING]}
- Clinical Monitoring: {category_counts[QCCategory.CLINICAL_MONITORING]}
- Post-Market Surveillance: {category_counts[QCCategory.POST_MARKET_SURVEILLANCE]}

This plan provides comprehensive coverage of quality, safety, and efficacy requirements
for therapeutic genetic intervention development and clinical application.
        """

        return conclusion.strip()

    def generate_regulatory_submission_package(
        self, validation_report: ValidationReport
    ) -> Dict[str, Any]:
        """
        Generate regulatory submission package.

        Args:
            validation_report: Comprehensive validation report

        Returns:
            Regulatory submission package
        """
        submission_package = {
            "executive_summary": {
                "product_name": validation_report.target_product,
                "indication": "therapeutic_genetic_correction",
                "validation_status": validation_report.approval_status,
                "submission_date": datetime.now().isoformat(),
            },
            "quality_control_summary": {
                "total_tests": len(validation_report.test_results),
                "critical_tests": sum(
                    1 for test in validation_report.test_results if test.is_critical
                ),
                "test_categories": list(
                    set(
                        test.test_category.value
                        for test in validation_report.test_results
                    )
                ),
            },
            "manufacturing_information": {
                "gmp_compliance": "yes",
                "quality_system": "ISO_13485",
                "manufacturing_sites": "to_be_determined",
            },
            "clinical_plan": {
                "phase_1_safety": "dose_escalation_design",
                "primary_endpoints": ["safety", "preliminary_efficacy"],
                "monitoring_plan": "comprehensive_as_per_protocol",
            },
            "risk_management": {
                "risk_evaluation": "comprehensive_rems_if_required",
                "mitigation_strategies": validation_report.recommendations,
                "monitoring_strategy": "long_term_surveillance_plan",
            },
            "regulatory_pathway": {
                "fda_pathway": "BLA_biologics_license_application",
                "ema_pathway": "CAT_committee_advanced_therapies",
                "breakthrough_designation": "to_be_requested",
            },
        }

        return submission_package
