"""Example of running a COMT gene analysis."""
from enhancement_engine import EnhancementEngine


def main() -> None:
    engine = EnhancementEngine(email="demo@example.com")
    report = engine.analyze_gene("COMT", "Val158Met")

    print(report.summary)

    guide = report.best_guide
    print(f"Best guide sequence: {guide.sequence}")
    print(f"Guide efficiency: {guide.efficiency_score.overall_efficiency:.2f}")
    print(f"Off-target sites: {len(guide.off_targets)}")

    safety = report.safety_assessment
    print(f"Overall safety score: {safety.overall_score:.1f}")

    gain = report.predicted_effect.enhancement_gain.improvement_factor
    print(f"Expected improvement factor: {gain:.1f}x")


if __name__ == "__main__":
    main()
