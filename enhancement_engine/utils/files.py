import json
import csv
from pathlib import Path
from typing import Any, Dict, List


def load_json(path: str) -> Any:
    """Load a JSON file and return the parsed data."""
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def save_json(data: Any, path: str) -> None:
    """Save data as JSON to the provided path."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


def load_text(path: str) -> str:
    """Read a text file and return its content."""
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def save_text(text: str, path: str) -> None:
    """Write text content to a file."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def load_lines(path: str) -> List[str]:
    """Return a list of stripped lines from a text file."""
    with open(path, "r", encoding="utf-8") as fh:
        return [line.rstrip("\n") for line in fh]


def load_csv(path: str) -> List[Dict[str, str]]:
    """Load a CSV file into a list of dictionaries."""
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return list(reader)


def save_csv(rows: List[Dict[str, Any]], path: str, fieldnames: List[str]) -> None:
    """Save a list of dictionaries to a CSV file."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
