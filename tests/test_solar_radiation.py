from environmental_biophysics.solar_radiation import get_height_weight_factor


def test_get_height_weight_factor():
    assert get_height_weight_factor(0.75, 1.52, 0.56, 3) == 0.89
