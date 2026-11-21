"""Shared fixtures; tests exercise installed/public APIs, not notebook globals."""

from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")


@pytest.fixture
def reference():
    def read(name):
        return pd.read_csv(Path(__file__).parent / "fixtures" / name)
    return read
