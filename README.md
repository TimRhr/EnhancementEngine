# 🧬 Enhancement Engine

**Comprehensive simulation and analysis tool for human genetic enhancement**

Enhancement Engine is a Python framework for simulating, analyzing, and evaluating genetic modifications for human enhancement. It combines CRISPR guide design, off-target prediction, safety assessment, and phenotype simulation in a unified platform.

## ⚡ Quick Start

```python
from enhancement_engine import EnhancementEngine

# Initialize the engine
engine = EnhancementEngine(email="your.email@domain.com")

# Analyze an enhancement gene
report = engine.analyze_gene("COMT", "Val158Met")

print(f"Enhancement gain: {report.enhancement_gain}")
print(f"Safety score: {report.safety_assessment.safety_score}/100")
```

## 🎯 Key Features

- **🔬 CRISPR Guide Design**: Optimal guide RNA selection with efficiency scoring
- **⚠️ Safety Assessment**: Comprehensive off-target and risk analysis
- **📊 Effect Simulation**: Predict enhancement outcomes and side effects
- **🧬 Gene Database**: Curated database of enhancement-relevant genes
- **📈 Batch Analysis**: Process multiple genes and variants simultaneously
- **📋 Detailed Reporting**: Publication-ready analysis reports

## 🚀 Installation

### From GitHub
```bash
git clone https://github.com/TimRhr/EnhancementEngine.git
cd EnhancementEngine
pip install -e .
```

### Dependencies
```bash
pip install -r requirements.txt
```

## 📖 Usage Examples

### Single Gene Analysis
```python
from enhancement_engine import EnhancementEngine

engine = EnhancementEngine(email="researcher@university.edu")

# Analyze COMT gene for cognitive enhancement
result = engine.analyze_gene(
    gene_name="COMT",
    variant="Val158Met",
    target_tissue="brain"
)

print(f"Predicted cognitive improvement: {result.enhancement_gain.cognitive}")
print(f"Off-target sites found: {len(result.safety_assessment.off_targets)}")
```

### Batch Analysis
```python
# Analyze multiple enhancement genes
enhancement_genes = ["COMT", "BDNF", "ACTN3", "FOXO3"]
batch_report = engine.batch_analysis(enhancement_genes)

# Compare enhancement potentials
for gene, report in batch_report.results.items():
    print(f"{gene}: Safety={report.safety_score}, Effect={report.enhancement_gain}")
```

### Custom CRISPR Design
```python
from enhancement_engine.core import CRISPRDesigner

designer = CRISPRDesigner(cas_type="cas9")
guides = designer.design_guides(
    sequence=gene_sequence,
    target_region=(100, 200)
)

# Get best guide with safety metrics
best_guide = guides[0]
print(f"Guide: {best_guide.sequence}")
print(f"Efficiency: {best_guide.efficiency_score}")
print(f"Safety: {best_guide.specificity_score}")
```

### Command-Line Interface

Enhancement Engine also provides a console entry point named
`enhancement-engine`. Some example commands are shown below:

```bash
# Show available commands
enhancement-engine --help

# Analyze a single gene
enhancement-engine analyze COMT -e your.email@domain.com

# Analyze multiple genes at once
enhancement-engine batch COMT,BDNF,ACTN3 -e your.email@domain.com

# Inspect cached data statistics
enhancement-engine cache stats -e your.email@domain.com
```

### Web Interface

An optional web-based UI built with [Flask](https://flask.palletsprojects.com/) is included.
Install the extra dependency and launch the app:

```bash
pip install Flask
python examples/webapp.py
```

The server runs on `http://127.0.0.1:5000/`. Open this URL in a browser to use
the HTML form or POST a JSON payload with `gene` and `variant` to `/analyze` to
retrieve an analysis report. A brief walkthrough is available in
[docs/webapp_guide.md](docs/webapp_guide.md).

## 🩺 Therapeutic Engine

The project also includes a therapeutic pipeline aimed at disease correction.
The `TherapeuticEnhancementEngine` integrates risk assessment, CRISPR design and
clinical safety modules.

### CLI Usage

Run a therapeutic analysis from the command line with:

```bash
enhancement-engine therapeutic PTPN22 --variant R620W \
    --disease rheumatoid_arthritis -e you@example.com
```

### Webapp

Start the Flask app and navigate to `/therapeutic` for an HTML form or POST JSON
to `/api/therapeutic`.

## 🧬 Supported Enhancement Genes

| Gene | Function | Enhancement Type |
|------|----------|------------------|
| **COMT** | Dopamine metabolism | Cognitive performance |
| **BDNF** | Neuroplasticity | Learning & memory |
| **ACTN3** | Muscle fiber type | Athletic performance |
| **FOXO3** | Cellular longevity | Lifespan extension |
| **EPO** | Oxygen transport | Endurance |
| **MYOSTATIN** | Muscle growth inhibitor | Strength & muscle mass |

## 📊 Analysis Pipeline

```
Gene Input → Sequence Analysis → CRISPR Design → Safety Assessment → Effect Simulation → Report Generation
```

## ⚠️ Safety & Ethics

This tool is designed for **research and educational purposes**. Any practical application of genetic enhancement technologies should:

- Follow all applicable laws and regulations
- Undergo proper ethical review
- Include comprehensive safety testing
- Consider long-term consequences
- Respect individual autonomy and consent

## 🤝 Contributing

We welcome contributions! Please read our [CONTRIBUTING.md](CONTRIBUTING.md) for code style, testing instructions and the pull request workflow.

### Development Setup
```bash
git clone https://github.com/TimRhr/EnhancementEngine.git
cd EnhancementEngine
pip install -e ".[dev]"
pre-commit install
```

### Running Tests
```bash
pytest tests/ --cov=enhancement_engine
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Documentation

- [User Guide](docs/user_guide.md)
- [API Reference](docs/api_reference.md)
- [Therapeutic Engine](docs/therapeutic_engine.md)
- [Examples](examples/)
- [FAQ](docs/faq.md)

## 🔗 Related Tools

- [CRISPOR](http://crispor.org/) - CRISPR guide design
- [Cas-OFFinder](http://www.rgenome.net/cas-offinder/) - Off-target prediction
- [AlphaFold](https://alphafold.ebi.ac.uk/) - Protein structure prediction

## 📧 Contact

- **Author**: Tim Röhr
- **GitHub**: [@TimRhr](https://github.com/TimRhr)
- **Issues**: [GitHub Issues](https://github.com/TimRhr/EnhancementEngine/issues)

## ⭐ Acknowledgments

Built with:
- [BioPython](https://biopython.org/) for biological data processing
- [scikit-learn](https://scikit-learn.org/) for machine learning
- [Plotly](https://plotly.com/) for interactive visualizations

---

**Disclaimer**: This software is for research purposes only. The authors are not responsible for any misuse of this technology.