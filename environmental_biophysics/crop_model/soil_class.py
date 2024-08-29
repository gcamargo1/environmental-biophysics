import numpy as np
import pandas as pd

from environmental_biophysics.soil import (
    FIELD_CAPACITY_WATER_POT,
    PERMNT_WILT_POINT_WATER_POT,
    WATER_DENSITY,
    get_air_entry_pot,
    get_b_value,
    get_bulk_density,
    get_organic_matter,
    get_sat_water_content,
    get_vol_water_content_33_j_kg,
    get_vol_water_content_1500_jkg,
    get_water_content,
    get_water_potential,
)


class Soil:
    """Soil instance."""

    def __init__(
        self,
        soil_df: pd.DataFrame,  # TODO (Gustavo): validate soil_df dataframe
        daily_ref_evap_transp: float = 5,
        *,
        manual_organic_matter: bool = False,
        manual_bulk_density: bool = False,
        manual_field_capacity: bool = False,
        manual_permanent_wilt_point: bool = False,
    ) -> None:
        self.total_layers = soil_df.shape[0]
        self.layers = range(self.total_layers)
        self.daily_ref_evap_transp = daily_ref_evap_transp  # weather dependent
        self.manual_organic_matter = manual_organic_matter
        self.manual_bulk_density = manual_bulk_density
        self.manual_field_capacity = manual_field_capacity
        self.manual_permanent_wilt_point = manual_permanent_wilt_point

        # Arrays declaration
        self.layer_thickness = np.zeros(self.total_layers)
        self.cum_depth = np.zeros(self.total_layers)
        self.clay = np.zeros(self.total_layers)
        self.sand = np.zeros(self.total_layers)
        self.organic_matter = np.zeros(self.total_layers)
        self.bulk_density = np.zeros(self.total_layers)  # Mg/m3
        self.field_capacity = np.zeros(self.total_layers)  # m3/m3
        self.mean_field_capacity = 0
        self.perm_wilt_point = np.zeros(self.total_layers)  # m3/m3
        self.porosity = np.zeros(self.total_layers)  # m3/m3
        self.b_value = np.zeros(self.total_layers)
        self.kl = np.zeros(self.total_layers)  # cm2 / day
        self.air_entry_potential = np.zeros(self.total_layers)
        self.field_capacity_water_potential = np.zeros(self.total_layers)
        self.plant_avail_water = np.zeros(self.total_layers)
        self.init_plant_avail_water = np.zeros(self.total_layers)
        self.water_content = np.zeros(self.total_layers)
        self.init_water_avail = np.zeros(self.total_layers)
        self.water_potential = np.zeros(self.total_layers)
        self.perm_wilt_point_pot = np.zeros(self.total_layers)
        self.sat_water_potential = np.zeros(self.total_layers)

        for lyr in self.layers:
            self.layer_thickness[lyr] = soil_df["soil_layer_thickness_m"][lyr]
            self.cum_depth[lyr] = soil_df["cumulative_depth_m"][lyr]
            self.clay[lyr] = soil_df["clay_percentage"][lyr] / 100
            self.sand[lyr] = soil_df["sand_percentage"][lyr] / 100
            if self.manual_organic_matter:
                self.organic_matter[lyr] = soil_df["organic_matter_percentage"][lyr]
            else:
                self.organic_matter[lyr] = get_organic_matter(self.clay[lyr])
            if self.manual_bulk_density:
                self.bulk_density[lyr] = soil_df["bulk_density_Mg_m3"][lyr]
            else:
                self.bulk_density[lyr] = get_bulk_density(
                    self.clay[lyr], self.sand[lyr], self.organic_matter[lyr]
                )
            self.porosity[lyr] = get_sat_water_content(self.bulk_density[lyr])
            if self.manual_field_capacity:
                self.field_capacity[lyr] = soil_df["field_capacity_m3_m3"][lyr]
            else:
                self.field_capacity[lyr] = get_vol_water_content_33_j_kg(
                    self.clay[lyr], self.sand[lyr], self.organic_matter[lyr]
                )
            if self.manual_permanent_wilt_point:
                self.perm_wilt_point[lyr] = soil_df["perm_wilt_pt_m3_m3"][lyr]
            else:
                self.perm_wilt_point[lyr] = get_vol_water_content_1500_jkg(
                    self.clay[lyr], self.sand[lyr], self.organic_matter[lyr]
                )
            self.b_value[lyr] = get_b_value(
                self.field_capacity[lyr], self.perm_wilt_point[lyr]
            )
            self.air_entry_potential[lyr] = get_air_entry_pot(
                self.field_capacity[lyr], self.porosity[lyr], self.b_value[lyr]
            )
            if self.manual_field_capacity:  # water content
                self.field_capacity_water_potential[lyr] = get_water_potential(
                    sat_water_content=self.porosity[lyr],
                    air_entry_potential=self.air_entry_potential[lyr],
                    campbell_b=self.b_value[lyr],
                    water_content=self.field_capacity[lyr],
                )
            else:
                # calculated again using more accurate formula
                self.field_capacity[lyr] = get_water_content(
                    self.porosity[lyr],
                    self.air_entry_potential[lyr],
                    self.b_value[lyr],
                    water_potential=FIELD_CAPACITY_WATER_POT,
                )
            if not self.manual_permanent_wilt_point:
                # calculated again using more accurate formula
                self.perm_wilt_point[lyr] = get_water_content(
                    self.porosity[lyr],
                    self.air_entry_potential[lyr],
                    self.b_value[lyr],
                    PERMNT_WILT_POINT_WATER_POT,
                )
            assert (
                self.perm_wilt_point[lyr] < self.field_capacity[lyr]
            ), "Permanent wilting point must be less than field capacity"
            self.plant_avail_water[lyr] = (
                self.field_capacity[lyr] - self.perm_wilt_point[lyr]
            )
            self.init_plant_avail_water[lyr] = soil_df[
                "initial_plant_avail_water_fraction"
            ][lyr]
            self.water_content[lyr] = (
                self.init_plant_avail_water[lyr] * self.plant_avail_water[lyr]
                + self.perm_wilt_point[lyr]
            )
            self.init_water_avail[lyr] = (
                (self.water_content[lyr] - self.perm_wilt_point[lyr])
                * self.layer_thickness[lyr]
                * WATER_DENSITY
            )
            self.water_potential[lyr] = get_water_potential(
                self.porosity[lyr],
                self.air_entry_potential[lyr],
                self.b_value[lyr],
                self.water_content[lyr],
            )
            self.perm_wilt_point_pot[lyr] = get_water_potential(
                self.porosity[lyr],
                self.air_entry_potential[lyr],
                self.b_value[lyr],
                self.perm_wilt_point[lyr],
            )
            self.sat_water_potential[lyr] = get_water_potential(
                self.porosity[lyr],
                self.air_entry_potential[lyr],
                self.b_value[lyr],
                self.porosity[lyr],
            )

        for lyr in self.layers:
            self.mean_field_capacity += (
                self.field_capacity[lyr] * self.layer_thickness[lyr]
            )
        self.mean_field_capacity = self.mean_field_capacity / self.layer_thickness.sum()

    def update_water_content(self, crop_list) -> None:  # noqa: ANN001
        """Updates soil water content based on each crop water uptake."""
        for crop in crop_list:
            if True:  # crop.planted_flag:
                for lyr in self.layers:
                    self.water_content[lyr] = self.water_content[
                        lyr
                    ] - crop.water_uptake[lyr] / (
                        self.layer_thickness[lyr] * WATER_DENSITY
                    )
                    self.water_potential[lyr] = get_water_potential(
                        self.porosity[lyr],
                        self.air_entry_potential[lyr],
                        self.b_value[lyr],
                        self.water_content[lyr],
                    )
