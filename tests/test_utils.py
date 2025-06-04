import pytest

pytest.importorskip("Bio")
pytest.importorskip("numpy")
pytest.importorskip("pandas")

from enhancement_engine.models.data_classes import (
    EfficiencyScore,
    GuideRNA,
    PAMSite,
    OffTarget,
    SafetyScore,
)


def test_efficiency_score_average():
    score = EfficiencyScore(on_target_score=0.5, doench_score=0.7)
    assert abs(score.overall_efficiency - 0.6) < 1e-6


def test_guide_rna_safety_score():
    pam = PAMSite(sequence="NGG", position=1, strand="+", cas_type=None)
    off1 = OffTarget(sequence="A" * 20, chromosome="1", position=1, strand="+", mismatches=1)
    off2 = OffTarget(sequence="A" * 20, chromosome="1", position=2, strand="+", mismatches=3, essential_gene=True)
    guide = GuideRNA(
        sequence="A" * 20,
        pam_site=pam,
        target_position=1,
        efficiency_score=EfficiencyScore(on_target_score=0.9),
        gc_content=50.0,
        off_targets=[off1, off2],
    )
    expected = 100 - 20 - 5 - 30
    assert guide.overall_safety_score == expected


def test_safety_score_risk_level():
    score = SafetyScore(75, 70, 70, 70, 70)
    assert score.risk_level.value == "low"
