import os
import pytest


@pytest.fixture
def test_data_dir():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "fixtures"))
