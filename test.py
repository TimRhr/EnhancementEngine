from enhancement_engine import EnhancementEngine

# Initialize the engine
engine = EnhancementEngine(email="your.email@domain.com")

# Analyze an enhancement gene
report = engine.analyze_gene("COMT", "Val158Met")

print(f"Enhancement gain: {report.enhancement_gain}")
print(f"Safety score: {report.safety_assessment.safety_score}/100")