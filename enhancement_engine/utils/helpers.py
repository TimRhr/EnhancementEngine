import re
from datetime import datetime


def slugify(value: str) -> str:
    """Simple slugify helper for filenames or identifiers."""
    value = value.lower()
    value = re.sub(r"[^a-z0-9]+", "-", value)
    return value.strip("-")


def current_timestamp() -> str:
    """Return a timestamp string suitable for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")
