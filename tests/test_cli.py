import subprocess
import sys
import pytest


def test_cli_help():
    pytest.importorskip("Bio")
    pytest.importorskip("pandas")
    pytest.importorskip("numpy")
    result = subprocess.run([
        sys.executable,
        "-m",
        "enhancement_engine.cli",
        "--help",
    ], capture_output=True, text=True)
    assert result.returncode == 0
    assert "usage" in result.stdout.lower()
