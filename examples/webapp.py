"""Minimal Flask web interface for Enhancement Engine."""

from flask import Flask, jsonify, request

from enhancement_engine import EnhancementEngine

app = Flask(__name__)
engine = EnhancementEngine(email="demo@example.com")


@app.route("/")
def index() -> str:
    """Simple index route."""
    return "Enhancement Engine Web API"


@app.route("/analyze", methods=["POST"])
def analyze() -> tuple:
    """Run a gene analysis and return JSON."""
    data = request.get_json() or {}
    gene = data.get("gene")
    variant = data.get("variant")
    if not gene or not variant:
        return jsonify({"error": "gene and variant required"}), 400

    report = engine.analyze_gene(gene, variant)
    return (
        jsonify(
            gene=gene,
            variant=variant,
            enhancement_gain=report.enhancement_gain,
            safety_score=report.safety_assessment.safety_score,
        ),
        200,
    )


if __name__ == "__main__":
    app.run(debug=True)

