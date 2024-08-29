from environmental_biophysics.crop_model.crop_class import Crop
from environmental_biophysics.crop_model.soil_class import Soil


def test_crop_class(crop_df, soil_df):
    soil = Soil(soil_df=soil_df)
    crop = Crop(crop_no=0, sim_length=100, crop_df=crop_df, soil=soil)
    assert crop is not None
