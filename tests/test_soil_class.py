from environmental_biophysics.crop_model.soil_class import Soil


def test_soil_class(soil_df):
    soil = Soil(soil_df=soil_df)
    # TODO(gustavo): test update water content in the profile,
    #  once crop list is available.
    assert soil is not None
