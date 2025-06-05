#!/usr/bin/env python3
"""
Run script for Enhancement Engine Web Application.

This script starts the Flask development server for the Enhancement Engine
web interface.
"""

import os
import sys
from webapp import app

def main():
    """Start the Flask development server."""
    print("🧬 Starting Enhancement Engine Web Application...")
    print("=" * 50)
    print("Application will be available at:")
    print("📍 http://127.0.0.1:5000/")
    print("📍 http://localhost:5000/")
    print("=" * 50)
    print("💡 Available endpoints:")
    print("   GET  /           - Analysis form")
    print("   POST /analyze    - Run analysis")
    print("   POST /api/analyze - JSON API endpoint")
    print("   GET  /health     - Health check")
    print("=" * 50)
    print("⚠️  Safety Notice:")
    print("   This tool is for research and educational purposes only.")
    print("   Any practical application requires proper ethical review,")
    print("   safety testing, and regulatory approval.")
    print("=" * 50)
    
    try:
        # Set debug mode based on environment
        debug_mode = os.getenv('FLASK_DEBUG', 'True').lower() == 'true'
        
        # Start the Flask app
        app.run(
            host='127.0.0.1',
            port=5000,
            debug=debug_mode,
            use_reloader=debug_mode,
            threaded=True
        )
        
    except KeyboardInterrupt:
        print("\n🛑 Application stopped by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Error starting application: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()