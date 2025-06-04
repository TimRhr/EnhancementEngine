from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="enhancement-engine",
    version="0.1.0",
    author="Tim Röhr",
    author_email="tim.roehr@outlook.com",
    description="Comprehensive simulation and analysis tool for human genetic enhancement",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TimRhr/EnhancementEngine",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
    python_requires=">=3.8",
    install_requires=requirements + ["Flask"],
    extras_require={
        "dev": ["pytest", "pytest-cov", "black", "flake8"],
        "ml": ["torch", "tensorflow"],
        "md": ["mdanalysis"],
    },
    entry_points={
        "console_scripts": [
            "enhancement-engine=enhancement_engine.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "enhancement_engine": ["data/*.json", "data/*.csv"],
    },
)