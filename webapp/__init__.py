from flask import Flask, render_template, request
from enhancement_engine import EnhancementEngine

app = Flask(__name__)

# For demo purposes we create the engine once. In a real app this might be
# configured differently or use dependency injection.
engine = EnhancementEngine(email="demo@example.com")


@app.route("/", methods=["GET"])
def index() -> str:
    """Render the input form."""
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze() -> str:
    """Run a simple gene analysis and display the result."""
    gene = request.form.get("gene", "").strip()
    variant = request.form.get("variant", "enhancement_variant").strip()
    if not gene:
        return render_template("index.html", error="Gene name is required")

    report = engine.analyze_gene(gene, variant)
    return render_template("result.html", report=report)


__all__ = ["app"]
