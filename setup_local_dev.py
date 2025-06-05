#!/usr/bin/env python3
"""
Setup script for local Enhancement Engine development.
"""

import os
import sys
from pathlib import Path

def create_directories():
    """Create necessary directories for local development."""
    directories = [
        'logs',
        'data/cache',
        'data/cached_sequences',
        'templates',
        'static'
    ]
    
    base_dir = Path(__file__).parent
    
    for dir_path in directories:
        full_path = base_dir / dir_path
        full_path.mkdir(parents=True, exist_ok=True)
        print(f"✓ Created directory: {full_path}")

def create_env_file():
    """Create .env file if it doesn't exist."""
    env_path = Path(__file__).parent / '.env'
    env_example_path = Path(__file__).parent / '.env.example'
    
    if not env_path.exists():
        if env_example_path.exists():
            # Copy from .env.example
            with open(env_example_path, 'r') as src:
                content = src.read()
            with open(env_path, 'w') as dst:
                dst.write(content)
            print(f"✓ Created .env file from .env.example")
        else:
            # Create basic .env file
            env_content = """# Enhancement Engine Local Development
FLASK_HOST=127.0.0.1
FLASK_PORT=5000
FLASK_DEBUG=true
SECRET_KEY=local-development-key
DEMO_EMAIL=demo@example.com
CONTAINER_ENV=false
LOG_LEVEL=INFO
"""
            with open(env_path, 'w') as f:
                f.write(env_content)
            print(f"✓ Created basic .env file")
    else:
        print(f"ℹ .env file already exists")

def check_dependencies():
    """Check if required dependencies are installed."""
    required_packages = [
        'flask',
        'biopython',
        'pandas', 
        'numpy',
        'scipy',
        'requests'
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"✓ {package} is installed")
        except ImportError:
            missing_packages.append(package)
            print(f"✗ {package} is missing")
    
    if missing_packages:
        print(f"\n⚠️  Missing packages: {', '.join(missing_packages)}")
        print(f"Install with: pip install {' '.join(missing_packages)}")
        return False
    
    return True

def create_test_templates():
    """Create basic templates if they don't exist."""
    templates_dir = Path(__file__).parent / 'templates'
    
    # Check if templates exist
    required_templates = ['index.html', 'result.html', 'error.html']
    missing_templates = []
    
    for template in required_templates:
        template_path = templates_dir / template
        if not template_path.exists():
            missing_templates.append(template)
    
    if missing_templates:
        print(f"⚠️  Missing templates: {', '.join(missing_templates)}")
        print("Please make sure the template files exist in the templates/ directory")
        return False
    
    print("✓ All required templates found")
    return True

def main():
    """Main setup function."""
    print("🧬 Enhancement Engine - Local Development Setup")
    print("=" * 50)
    
    try:
        # Create directories
        print("\n📁 Creating directories...")
        create_directories()
        
        # Create .env file
        print("\n⚙️  Setting up configuration...")
        create_env_file()
        
        # Check dependencies
        print("\n📦 Checking dependencies...")
        deps_ok = check_dependencies()
        
        # Check templates
        print("\n📄 Checking templates...")
        templates_ok = create_test_templates()
        
        print("\n" + "=" * 50)
        
        if deps_ok and templates_ok:
            print("✅ Setup complete! You can now run:")
            print("   python webapp/run.py")
            print("\n🌐 The app will be available at: http://127.0.0.1:5000")
        else:
            print("❌ Setup incomplete. Please resolve the issues above.")
            sys.exit(1)
            
    except Exception as e:
        print(f"❌ Setup failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()