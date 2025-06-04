from flask import Flask, request, jsonify

from .core.engine import EnhancementEngine

app = Flask(__name__)

# Create a single engine instance for simplicity
_engine = EnhancementEngine(email="webapp@example.com")


@app.route("/")
def index():
    """Simple index route."""
    return "Enhancement Engine", 200


@app.route("/analyze")
def analyze():
    """Analyze a gene via query parameters."""
    gene = request.args.get("gene")
    variant = request.args.get("variant", "enhancement_variant")
    target = request.args.get("target", "general")

    if not gene:
        return jsonify({"error": "gene parameter required"}), 400

    report = _engine.analyze_gene(gene, variant, target)
    return jsonify({"gene_name": report.gene_name})
