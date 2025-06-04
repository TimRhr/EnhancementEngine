# User Guide

This guide covers the basic setup steps for the Enhancement Engine and shows how to run your first analyses.

## Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/TimRhr/EnhancementEngine.git
   cd EnhancementEngine
   ```
2. **Install the package**
   Install the library in editable mode so the command line interface is available:
   ```bash
   pip install -e .
   ```
3. **Install dependencies**
   All required third‑party packages can be installed with:
   ```bash
   pip install -r requirements.txt
   ```

## Basic Usage

Create an instance of `EnhancementEngine` and run an analysis:

```python
from enhancement_engine import EnhancementEngine

engine = EnhancementEngine(email="your.email@domain.com")
report = engine.analyze_gene("COMT", "Val158Met")

print(f"Enhancement gain: {report.enhancement_gain}")
print(f"Safety score: {report.safety_assessment.safety_score}/100")
```

The package also provides a command line entry point:

```bash
# Show available commands
enhancement-engine --help

# Analyze a gene
enhancement-engine analyze COMT -e your.email@domain.com
```
