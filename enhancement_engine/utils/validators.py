import re
from pathlib import Path


def is_valid_email(email: str) -> bool:
    """Basic email validation."""
    if not email:
        return False
    pattern = r"^[^@\s]+@[^@\s]+\.[^@\s]+$"
    return re.match(pattern, email) is not None


def validate_email(email: str) -> None:
    """Raise ValueError if the email address is invalid."""
    if not is_valid_email(email):
        raise ValueError(f"Invalid email address: {email}")


def validate_file(path: str) -> None:
    """Raise FileNotFoundError if the given path does not exist."""
    if not Path(path).is_file():
        raise FileNotFoundError(f"File not found: {path}")
