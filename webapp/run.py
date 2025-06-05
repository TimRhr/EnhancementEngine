#!/usr/bin/env python3
"""
Enhanced Flask web interface for the Enhancement Engine.
Supports both local development and containerized deployment.
"""

import os
import sys
import logging
from pathlib import Path
from typing import Optional

# Add the parent directory to Python path for local development
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

from flask import Flask, render_template, request, jsonify, flash, redirect, url_for
from enhancement_engine import EnhancementEngine
from enhancement_engine.core.engine import EnhancementEngineError


def setup_logging() -> None:
    """Setup logging with platform-independent paths."""
    # Determine log directory based on environment
    if os.getenv('CONTAINER_ENV', 'false').lower() == 'true':
        # Container environment
        log_dir = Path('/app/logs')
    else:
        # Local development
        log_dir = Path(__file__).parent.parent / 'logs'
    
    # Create log directory if it doesn't exist
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging configuration
    log_file = log_dir / 'webapp.log'
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )


def create_app(config: Optional[dict] = None) -> Flask:
    """Create and configure Flask application."""
    # Determine template directory
    if os.getenv('CONTAINER_ENV', 'false').lower() == 'true':
        template_dir = '/app/templates'
    else:
        # Local development - check both possible locations
        local_templates = Path(__file__).parent.parent / 'templates'
        webapp_templates = Path(__file__).parent / 'templates'
        
        if local_templates.exists():
            template_dir = str(local_templates)
        elif webapp_templates.exists():
            template_dir = str(webapp_templates)
        else:
            # Create templates directory if it doesn't exist
            local_templates.mkdir(parents=True, exist_ok=True)
            template_dir = str(local_templates)
    
    app = Flask(__name__, template_folder=template_dir)
    
    # Configuration
    app.config.update({
        'SECRET_KEY': os.getenv('SECRET_KEY', 'development-key-change-in-production'),
        'DEBUG': os.getenv('FLASK_DEBUG', 'true').lower() == 'true',
        'HOST': os.getenv('FLASK_HOST', '127.0.0.1'),
        'PORT': int(os.getenv('FLASK_PORT', '5000')),
    })
    
    if config:
        app.config.update(config)
    
    # Demo email for development (should be configurable in production)
    demo_email = os.getenv('DEMO_EMAIL', 'demo@example.com')
    
    # Initialize enhancement engine (with error handling)
    try:
        engine = EnhancementEngine(email=demo_email)
        app.logger.info(f"Enhancement Engine initialized with email: {demo_email}")
    except Exception as e:
        app.logger.error(f"Failed to initialize Enhancement Engine: {e}")
        engine = None
    
    @app.route("/", methods=["GET"])
    def index():
        """Render the input form."""
        return render_template("index.html")
    
    @app.route("/analyze", methods=["POST"])
    def analyze():
        """Run gene analysis and display results."""
        if not engine:
            flash("Enhancement Engine not available. Please check configuration.", "error")
            return redirect(url_for('index'))
        
        # Get form data
        gene = request.form.get("gene", "").strip()
        variant = request.form.get("variant", "enhancement_variant").strip()
        email = request.form.get("email", "").strip()
        
        # Validation
        if not gene:
            flash("Gene name is required", "error")
            return redirect(url_for('index'))
        
        if not variant:
            variant = "enhancement_variant"
        
        try:
            app.logger.info(f"Analyzing gene: {gene}, variant: {variant}")
            
            # Run analysis
            report = engine.analyze_gene(gene, variant)
            
            app.logger.info(f"Analysis completed for {gene}")
            
            return render_template(
                "result.html",
                gene=gene,
                variant=variant,
                report=report,
            )
            
        except EnhancementEngineError as e:
            app.logger.error(f"Enhancement analysis error: {e}")
            flash(f"Analysis failed: {str(e)}", "error")
            return redirect(url_for('index'))
        
        except Exception as e:
            app.logger.error(f"Unexpected error during analysis: {e}")
            flash("An unexpected error occurred. Please try again.", "error")
            return redirect(url_for('index'))
    
    @app.route("/api/analyze", methods=["POST"])
    def api_analyze():
        """JSON API endpoint for analysis."""
        if not engine:
            return jsonify({"error": "Enhancement Engine not available"}), 500
        
        try:
            data = request.get_json()
            if not data:
                return jsonify({"error": "No JSON data provided"}), 400
            
            gene = data.get("gene", "").strip()
            variant = data.get("variant", "enhancement_variant").strip()
            
            if not gene:
                return jsonify({"error": "Gene name is required"}), 400
            
            # Run analysis
            report = engine.analyze_gene(gene, variant)
            
            # Convert report to JSON-serializable format
            result = {
                "gene_name": report.gene_name,
                "target_variant": report.target_variant,
                "feasibility_score": report.feasibility_score,
                "safety_score": report.safety_assessment.overall_score,
                "enhancement_factor": report.predicted_effect.enhancement_gain.improvement_factor,
                "confidence": report.confidence_score,
                "summary": report.summary,
                "recommendations": report.recommendations[:3]  # Top 3 recommendations
            }
            
            return jsonify(result)
            
        except EnhancementEngineError as e:
            return jsonify({"error": str(e)}), 400
        except Exception as e:
            app.logger.error(f"API error: {e}")
            return jsonify({"error": "Internal server error"}), 500
    
    @app.route("/health", methods=["GET"])
    def health_check():
        """Health check endpoint."""
        status = {
            "status": "healthy" if engine else "degraded",
            "engine_available": engine is not None,
            "timestamp": "2025-01-03T10:00:00Z"
        }
        return jsonify(status)
    
    @app.errorhandler(404)
    def not_found(error):
        """Handle 404 errors."""
        return render_template("error.html", 
                             error_code=404, 
                             error_message="Page not found"), 404
    
    @app.errorhandler(500)
    def internal_error(error):
        """Handle 500 errors."""
        return render_template("error.html", 
                             error_code=500, 
                             error_message="Internal server error"), 500
    
    return app


def main():
    """Main application entry point."""
    # Setup logging
    setup_logging()
    
    # Create Flask app
    app = create_app()
    
    # Get configuration from environment
    host = app.config['HOST']
    port = app.config['PORT']
    debug = app.config['DEBUG']
    
    app.logger.info(f"Starting Enhancement Engine Web App")
    app.logger.info(f"Debug mode: {debug}")
    app.logger.info(f"Listening on: http://{host}:{port}")
    
    try:
        app.run(host=host, port=port, debug=debug)
    except Exception as e:
        app.logger.error(f"Failed to start application: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()