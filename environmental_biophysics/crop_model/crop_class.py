import numpy as np

from environmental_biophysics.soil import FIELD_CAPACITY_WATER_POT


class Crop:
    def __init__(self, crop_no, sim_length, crop_df, soil) -> None:  # noqa: ANN001
        """Crop instance."""
        crop_inputs = crop_df.loc[crop_no]
        # Campbell max canopy transpiration, mm/d:
        self.campbell_max_daily_transp = crop_inputs[
            "campbell_max_full_cover_transpiration_mm_day"
        ]
        # Campbell leaf water potential at onset of stomatal closure [J/kg]:
        self.leaf_water_pot_stress_onset = crop_inputs[
            "campbell_leaf_water_potential_at_onset_of_stomatal closure_j_kg"
        ]
        # Campbell leaf water potential at wilting point [J/kg]:
        self.leaf_water_pot_wilt_point = crop_inputs[
            "campbell_leaf_water_potential_at_wilting_point_j_kg"
        ]
        self.root_depth = crop_inputs["crop_root_depth_m"]
        self.leaf_water_pot = 0  # J/kg
        self.sim_length = sim_length  # d
        self.conductance = np.ones(soil.total_layers)
        self.water_uptake = np.zeros(soil.total_layers)
        self.leaf_water_potential = np.zeros(soil.total_layers)
        self.soil_water_pot_avg = 0
        self.transp_ratio = 0  # to quantify crop water stress
        self.crop_transp = 0
        self.pot_transp = 0
        self.att_transp = 0
        self.expect_transp = 0
        self.cum_transp = 0
        self.cum_pot_transp = 0
        self.light_intercpt = 0
        self.max_pot_transp = 0
        # TODO (Gustavo): putting crop root information with soil for now.

    def get_water_uptake(self, soil) -> None:  # noqa: ANN001
        """CropSyst/Campbell model daily water uptake.

        Args:
            soil: soil class with properties

        References:
            Campbell, G. S. 1985. Soil physics with BASIC: Transport models for
             soil-plant systems. Developments in soil science. Elsevier, Amsterdam
            Campbell, G. S. 1991. Simulation of water uptake by plant roots. In
             Modeling plant and soil systems, 273-285. J. Hanks, and J. T. Ritchie,
             eds. Madison, WI: ASA/CSSA/SSSA.
        """
        daily_ref_evap_transp = soil.daily_ref_evap_transp
        root_hydr_cond = np.zeros(soil.total_layers)
        shoot_hydr_cond = np.zeros(soil.total_layers)
        plant_hydr_cond = np.zeros(soil.total_layers)
        root_activity = np.zeros(soil.total_layers)
        root_cond_adj = np.zeros(soil.total_layers)
        tot_root_cond_adj = 0
        salinity_factor = np.zeros(soil.total_layers)
        soil_water_pot_avg = 0

        # Transpiration
        self.pot_transp = daily_ref_evap_transp * self.light_intercpt
        self.max_pot_transp = self.campbell_max_daily_transp * self.light_intercpt
        self.expect_transp = min(self.pot_transp, self.max_pot_transp)  # mm/day

        # Plant hydraulic conductance (kg s m-4)
        tot_plant_hydr_cond = self.max_pot_transp / (
            FIELD_CAPACITY_WATER_POT - self.leaf_water_pot_stress_onset
        )
        # assumption of 2/3 of plant hydraulic conductance is from roots
        tot_root_hydr_cond = tot_plant_hydr_cond / 0.65

        # assumption of 1/3 of plant hydraulic conductivity is from shoots
        tot_shoot_hydr_cond = tot_plant_hydr_cond / 0.35

        for lyr in soil.layers:
            root_activity[lyr] = 1
            salinity_factor[lyr] = 1
            root_cond_adj[lyr] = (
                root_activity[lyr] * soil.root_fraction[lyr] * salinity_factor[lyr]
            )
            root_hydr_cond[lyr] = tot_root_hydr_cond * root_cond_adj[lyr]
            tot_root_cond_adj += root_cond_adj[lyr]

        # Root, shoot and plant hydraulic conductance(kg s m-4)
        for lyr in soil.layers:
            if root_cond_adj[lyr] > 0:
                shoot_hydr_cond[lyr] = (
                    tot_shoot_hydr_cond * root_cond_adj[lyr] / tot_root_cond_adj
                )
                plant_hydr_cond[lyr] = (
                    root_hydr_cond[lyr]
                    * shoot_hydr_cond[lyr]
                    / (root_hydr_cond[lyr] + shoot_hydr_cond[lyr])
                )
            else:
                plant_hydr_cond[lyr] = 0

        tot_root_hydr_cond *= tot_root_cond_adj
        tot_plant_hydr_cond = (tot_root_hydr_cond * tot_shoot_hydr_cond) / (
            tot_root_hydr_cond + tot_shoot_hydr_cond
        )

        if tot_plant_hydr_cond > 0:
            for lyr in soil.layers:
                soil_water_pot_avg += soil.water_potential[lyr] * root_cond_adj[lyr]
            leaf_water_pot = (
                soil_water_pot_avg - self.expect_transp / tot_plant_hydr_cond
            )
            if leaf_water_pot < self.leaf_water_pot_stress_onset:
                leaf_water_pot = (
                    tot_plant_hydr_cond
                    * soil_water_pot_avg
                    * (
                        self.leaf_water_pot_stress_onset
                        - self.leaf_water_pot_wilt_point
                    )
                    + self.leaf_water_pot_wilt_point * self.expect_transp
                ) / (
                    tot_plant_hydr_cond
                    * (
                        self.leaf_water_pot_stress_onset
                        - self.leaf_water_pot_wilt_point
                    )
                    + self.expect_transp
                )
            if leaf_water_pot < self.leaf_water_pot_wilt_point:
                leaf_water_pot = self.leaf_water_pot_wilt_point
                self.att_transp = 0
                transp_ratio = self.att_transp / self.expect_transp

            elif leaf_water_pot < self.leaf_water_pot_stress_onset:
                self.att_transp = (
                    self.expect_transp
                    * (leaf_water_pot - self.leaf_water_pot_wilt_point)
                    / (
                        self.leaf_water_pot_stress_onset
                        - self.leaf_water_pot_wilt_point
                    )
                )
                transp_ratio = self.att_transp / self.expect_transp

            else:
                self.att_transp = self.expect_transp
                transp_ratio = 1
            # crop water uptake (kg/m2/d = mm/d)
            for lyr in soil.layers:
                self.water_uptake[lyr] = (
                    plant_hydr_cond[lyr]
                    * (soil.water_potential[lyr] - leaf_water_pot)
                    * transp_ratio
                )
                if self.water_uptake[lyr] < 0:
                    self.water_uptake[lyr] = 0
        self.crop_transp = self.water_uptake.sum()  # mm/day
        self.cum_transp += self.crop_transp
        self.cum_pot_transp += self.expect_transp
        self.transp_ratio = self.crop_transp / self.expect_transp
