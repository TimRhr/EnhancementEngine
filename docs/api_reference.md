# API Reference

This page provides a high level overview of the main classes and functions exposed by the Enhancement Engine.

## `EnhancementEngine`

Central class used to run analyses. Typical workflow:

```python
from enhancement_engine import EnhancementEngine
engine = EnhancementEngine(email="you@example.com")
report = engine.analyze_gene("COMT", "Val158Met")
```

Important methods include:

- `analyze_gene(gene_name, variant, target_tissue="general")` – run a single gene analysis.
- `batch_analysis(genes, variants=None)` – analyze multiple genes at once.
- `clear_cache()` – remove cached sequence data.
- `get_database_stats()` – show information about cached data.

## Data Classes

Several convenience data classes are re-exported from `enhancement_engine` for easy access:

- `GeneInfo` – basic information about a gene.
- `GuideRNA` – CRISPR guide with efficiency and specificity scores.
- `EnhancementReport` – combined results from a gene analysis.
- `SafetyScore` – evaluation of potential off‑targets.
- `VariantEffect` – predicted improvement factor and other effects.

## Command Line Interface

The `enhancement-engine` command exposes the same functionality as the Python API. Run `enhancement-engine --help` to see all available commands.
