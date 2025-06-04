from enhancement_engine import EnhancementEngine


def main() -> None:
    """Run a very small demonstration analysis."""
    engine = EnhancementEngine(email="demo@example.com")
    report = engine.analyze_gene("COMT", "Val158Met")

    print(f"Gene: {report.gene_name}")
    print(f"Feasibility score: {report.feasibility_score:.1f}")
    print(f"Safety score: {report.safety_assessment.overall_score:.1f}")
    gain = report.predicted_effect.enhancement_gain.improvement_factor
    print(f"Predicted improvement: {gain:.1f}x")


if __name__ == "__main__":
    main()
