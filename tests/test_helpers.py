import pytest
pytest.importorskip("Bio")

from enhancement_engine.utils.helpers import slugify, current_timestamp
from enhancement_engine.utils.validators import is_valid_email
import re


def test_slugify_basic():
    assert slugify('Hello World!') == 'hello-world'
    assert slugify('Python--Rocks') == 'python-rocks'


def test_current_timestamp_format():
    ts = current_timestamp()
    assert re.match(r"\d{8}_\d{6}", ts)


def test_is_valid_email():
    assert is_valid_email('user@example.com')
    assert not is_valid_email('invalid-email')
