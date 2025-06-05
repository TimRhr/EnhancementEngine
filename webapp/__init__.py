"""Flask web interface for the Enhancement Engine."""

import os
from flask import Flask, render_template, request, jsonify
from enhancement_engine import EnhancementEngine

# Template directory configuration - use the main templates folder
TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "templates")
app = Flask(__name__, template_folder=TEMPLATE_DIR)


@app.route("/", methods=["GET"])
def index() -> str:
    """Render the input form."""
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze():
    """Run a gene analysis and display the result."""
    try:
        # Get form data
        gene = request.form.get("gene", "").strip()
        variant = request.form.get("variant", "enhancement_variant").strip()
        email = request.form.get("email", "").strip()
        
        # Validation
        if not gene:
            return render_template("index.html", error="Gene name is required")
        
        if not email:
            return render_template("index.html", error="Email is required for NCBI access")
        
        # Validate email format
        import re
        if not re.match(r'^[^@\s]+@[^@\s]+\.[^@\s]+$', email):
            return render_template("index.html", error="Please enter a valid email address")
        
        # Create engine with user's email
        engine = EnhancementEngine(email=email)
        
        # Run analysis
        report = engine.analyze_gene(gene, variant)
        
        # Return results
        return render_template(
            "result.html",
            gene=gene,
            variant=variant,
            report=report,
        )
        
    except Exception as e:
        error_msg = f"Analysis failed: {str(e)}"
        return render_template("index.html", error=error_msg)


@app.route("/api/analyze", methods=["POST"])
def api_analyze():
    """API endpoint for JSON analysis requests."""
    try:
        data = request.get_json()
        
        if not data:
            return jsonify({"error": "No JSON data provided"}), 400
        
        gene = data.get("gene", "").strip()
        variant = data.get("variant", "enhancement_variant").strip()
        email = data.get("email", "").strip()
        
        if not gene:
            return jsonify({"error": "Gene name is required"}), 400
        
        if not email:
            return jsonify({"error": "Email is required"}), 400
        
        # Create engine and run analysis
        engine = EnhancementEngine(email=email)
        report = engine.analyze_gene(gene, variant)
        
        # Convert report to JSON-serializable format
        result = {
            "gene_name": report.gene_name,
            "target_variant": report.target_variant,
            "feasibility_score": round(report.feasibility_score, 2),
            "safety_score": round(report.safety_assessment.overall_score, 2),
            "enhancement_factor": round(report.predicted_effect.enhancement_gain.improvement_factor, 2),
            "confidence": round(report.confidence_score, 2),
            "summary": report.summary,
            "recommendations": report.recommendations[:3],  # Top 3 recommendations
            "warnings": report.warnings
        }
        
        return jsonify(result)
        
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/health", methods=["GET"])
def health_check():
    """Health check endpoint."""
    return jsonify({"status": "healthy", "service": "Enhancement Engine Web API"})


if __name__ == "__main__":
    app.run(debug=True)