#!/usr/bin/env python3
"""
Enhanced Flask application runner for Enhancement Engine.
Optimized for Docker deployment with proper configuration management.
"""

import os
import sys
import logging
from pathlib import Path

# Add the project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from webapp import app

def setup_logging():
    """Configure logging for the application."""
    log_level = os.getenv('LOG_LEVEL', 'INFO').upper()
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Create logs directory if it doesn't exist
    log_dir = Path('/app/logs')
    log_dir.mkdir(exist_ok=True)
    
    logging.basicConfig(
        level=getattr(logging, log_level),
        format=log_format,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_dir / 'enhancement_engine.log')
        ]
    )

def configure_app():
    """Configure Flask application settings."""
    # Basic Flask configuration
    app.config['DEBUG'] = os.getenv('FLASK_DEBUG', 'false').lower() == 'true'
    app.config['TESTING'] = False
    
    # Security configuration
    app.config['SECRET_KEY'] = os.getenv('SECRET_KEY', 'dev-key-change-in-production')
    
    # Enhancement Engine specific configuration
    app.config['DEMO_EMAIL'] = os.getenv('DEMO_EMAIL', 'demo@example.com')
    app.config['CACHE_ENABLED'] = os.getenv('CACHE_ENABLED', 'true').lower() == 'true'
    app.config['CACHE_DIR'] = os.getenv('CACHE_DIR', '/app/data/cache')
    
    # Ensure cache directory exists
    Path(app.config['CACHE_DIR']).mkdir(parents=True, exist_ok=True)
    
    return app

def main():
    """Main function to run the Flask application."""
    setup_logging()
    
    logger = logging.getLogger(__name__)
    logger.info("Starting Enhancement Engine Web Application")
    
    # Configure the application
    app = configure_app()
    
    # Get host and port from environment variables
    host = os.getenv('FLASK_HOST', '0.0.0.0')
    port = int(os.getenv('FLASK_PORT', 5000))
    debug = os.getenv('FLASK_DEBUG', 'false').lower() == 'true'
    
    logger.info(f"Starting server on {host}:{port}")
    logger.info(f"Debug mode: {debug}")
    logger.info(f"Demo email: {app.config['DEMO_EMAIL']}")
    logger.info(f"Cache enabled: {app.config['CACHE_ENABLED']}")
    
    try:
        # Run the Flask application
        app.run(
            host=host,
            port=port,
            debug=debug,
            threaded=True,
            use_reloader=False  # Disable reloader in production
        )
    except Exception as e:
        logger.error(f"Failed to start application: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()