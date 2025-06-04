import pytest


def test_simple_gene_analysis(monkeypatch):
    pytest.importorskip("Bio")
    pytest.importorskip("numpy")
    pytest.importorskip("pandas")
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from enhancement_engine.core.engine import EnhancementEngine
    from enhancement_engine.models.data_classes import (
        GeneInfo,
        GuideRNA,
        PAMSite,
        EfficiencyScore,
        OffTarget,
        SafetyScore,
        VariantInfo,
        EnhancementGain,
        VariantEffect,
        ProteinEffect,
        EnhancementCategory,
    )

    engine = EnhancementEngine("test@example.com")

    def fake_get_gene_information(name):
        return GeneInfo(
            name="Gene",
            symbol=name,
            gene_id="1",
            chromosome="1",
            start_pos=1,
            end_pos=1000,
            description="d",
            ncbi_id="1",
            refseq_id="NM_000000",
        )

    def fake_analyze_gene_sequence(info):
        seq = SeqRecord(Seq("ATG" + "A" * 297 + "TGA"), id="x")
        return {
            "sequence_available": True,
            "sequence_record": seq,
            "sequence_length": len(seq.seq),
            "coding_region": None,
            "composition": {},
            "regulatory_elements": [],
        }

    def fake_identify_enhancement_target(info, variant):
        return 10

    def fake_design_crispr_guides(info, position, seq_analysis):
        pam = PAMSite(sequence="NGG", position=position, strand="+", cas_type=engine.config.cas_type)
        guide = GuideRNA(
            sequence="A" * 20,
            pam_site=pam,
            target_position=position,
            efficiency_score=EfficiencyScore(on_target_score=0.9),
            gc_content=50.0,
            off_targets=[OffTarget(sequence="A" * 20, chromosome="1", position=1, strand="+", mismatches=0)],
        )
        return {"best_guide": guide, "alternative_guides": [], "design_summary": {"total_guides": 1}}

    def fake_assess_safety(guide, gene, variant):
        return SafetyScore(90, 90, 90, 90, 90)

    def fake_simulate_effects(info, variant, tissue):
        return VariantEffect(
            variant=VariantInfo(name=variant),
            protein_effect=ProteinEffect(stability_change=0.0, activity_change=1.0),
            enhancement_gain=EnhancementGain(
                category=EnhancementCategory.COGNITIVE,
                primary_metric="m",
                baseline_value=1.0,
                enhanced_value=2.0,
                improvement_factor=2.0,
            ),
            side_effects=[],
        )

    monkeypatch.setattr(engine, "_get_gene_information", fake_get_gene_information)
    monkeypatch.setattr(engine, "_analyze_gene_sequence", fake_analyze_gene_sequence)
    monkeypatch.setattr(engine, "_identify_enhancement_target", fake_identify_enhancement_target)
    monkeypatch.setattr(engine, "_design_crispr_guides", fake_design_crispr_guides)
    monkeypatch.setattr(engine, "_assess_safety", fake_assess_safety)
    monkeypatch.setattr(engine, "_simulate_enhancement_effects", fake_simulate_effects)
    monkeypatch.setattr(engine, "_generate_recommendations", lambda *a, **k: ["ok"])
    monkeypatch.setattr(engine, "_calculate_confidence_score", lambda *a, **k: 0.9)
    monkeypatch.setattr(engine, "_generate_warnings", lambda *a, **k: [])

    report = engine.analyze_gene("COMT")
    assert report.gene_name == "COMT"
    assert report.best_guide.sequence == "A" * 20
    assert report.feasibility_score > 0
    assert report.confidence_score == 0.9
