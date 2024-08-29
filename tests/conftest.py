import pandas as pd
import pytest

from environmental_biophysics.environmental_constants import (
    environmental_biophysics_proj_path,
)


@pytest.fixture()
def soil_df():
    return pd.read_csv(environmental_biophysics_proj_path / "tests/data/soil_data.csv")


@pytest.fixture()
def crop_df():
    return pd.read_csv(environmental_biophysics_proj_path / "tests/data/crop_data.csv")
