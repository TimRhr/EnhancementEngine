# Contributing to Enhancement Engine

Thank you for considering a contribution! We welcome improvements, bug fixes and new features. This document explains the basic workflow.

## Code Style

- Follow **PEP 8** and keep code formatted with [`black`](https://black.readthedocs.io/).
- Run [`flake8`](https://flake8.pycqa.org/) to check for style issues.
- Use 4 spaces per indentation level and include type hints and docstrings for public functions and classes.

## Development Setup

```bash
# clone your fork and install dependencies
pip install -e ".[dev]"
pre-commit install
```

Running `pre-commit` locally will format your code and run linters before each commit.

## Testing

All functionality must be covered by tests. Execute the full suite with:

```bash
pytest tests/ --cov=enhancement_engine
```

Ensure `flake8` and `black` run without errors:

```bash
flake8
black --check enhancement_engine tests
```

## Submitting Pull Requests

1. Fork the repository and create a branch based on `main`.
2. Make your changes and add tests.
3. Verify that `black`, `flake8` and `pytest` succeed.
4. Push your branch and open a pull request describing the change and linking any related issues.
5. One of the maintainers will review and merge after checks pass.

Thank you for helping to improve the project!
