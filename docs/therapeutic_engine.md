# Therapeutic Engine

The therapeutic pipeline extends the core engine to evaluate CRISPR-based disease correction strategies.
It performs disease risk analysis, designs correction guides and estimates clinical safety.

## Quick Example

```python
from enhancement_engine.core.therapeutic_engine import TherapeuticEnhancementEngine

engine = TherapeuticEnhancementEngine(email="you@example.com")
report = engine.analyze_disease_gene(
    "PTPN22", "R620W", {"age": 45}, disease="rheumatoid_arthritis"
)
print(report.summary)
```

### Command Line

```bash
enhancement-engine therapeutic PTPN22 --variant R620W \
    --disease rheumatoid_arthritis -e you@example.com
```

### Web Interface

Launch the app with `enhancement-engine-web` and open `/therapeutic` or POST to
`/api/therapeutic` with `gene`, `variant` and `disease` fields.
