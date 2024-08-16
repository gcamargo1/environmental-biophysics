"""Solar radiation functions."""

import math

import numpy as np
from numpy.typing import NDArray

from environmental_biophysics.environmental_constants import (
    SUNSET_OR_SUNRISE_SOLAR_ZENITH_ANGLE,
)


def get_height_weight_factor(
    height_dom: float | NDArray,
    dominant_fact: float | NDArray,
    suppressed_fact: float | NDArray,
    number_species: int,
) -> float | NDArray:
    """Height domimance weight factor.

    Height dominance weighing factor calculated as a linear interpolation
    between dominant and suppressed species for radiation interception between
    species.

    Args:
        height_dom: Species dominance factor based on height
        dominant_fact: Dominant weighing factor
        suppressed_fact: Suppressed weighing factor
        number_species: number of species

    References:
        Camargo, G.G.T. 2015. Ph.D. Dissertation. Penn State University.

    Examples:
        >>> res = get_height_weight_factor(height_dom=0.75, dominant_fact=1.52,
        ...                                suppressed_fact=0.56, number_species=3)
        >>> assert res == 0.89

        >>> get_height_weight_factor(1.5, 1.52, 0.56, 3)
        1.13
    """
    assert (height_dom and dominant_fact and suppressed_fact and number_species) > 0

    # One species or species of same height
    if height_dom == 1:
        return 1

    # If species in shorter than canopy average
    if height_dom < 1:
        return ((height_dom - 1) * (suppressed_fact - 1)) / (0 - 1) + 1

    # If species in taller than canopy average
    return 1 + (height_dom - 1) / (number_species - 1) * (dominant_fact - 1)


def get_optical_air_mass(atm_press: float, solar_zenith_angle: float) -> float:
    """Optical air mass.

    Return optical air mass, or the ratio of slant path length through the
     atmosphere to zenith path length.

    Args:
        atm_press: [kPa]
        solar_zenith_angle: [deg]

    References:
        Campbell, G. S., and J. M. Norman. 1998. Introduction to environmental
         biophysics. Springer, New York. Eq. 11.12

    Examples:
        >>> get_optical_air_mass(100, 50)
        1.5357589603755304

        >>> get_optical_air_mass(91.6, 30)
        1.0441319774485631
    """
    sea_level_atm_prssr = 101.3
    seven_thousant_meters_atm_prssr = 41
    assert seven_thousant_meters_atm_prssr < atm_press <= sea_level_atm_prssr
    assert solar_zenith_angle >= 0
    deg_to_rad = math.pi / 180.0
    solar_zenith_angle *= deg_to_rad
    return atm_press / (sea_level_atm_prssr * math.cos(solar_zenith_angle))


def get_solar_radiation_extinction_coeff_black_beam(
    solar_zenith_angle: float, x_area_ratio: float
) -> float:
    """Solar radiation extinction coefficient of a canopy of black leaves.

    Return solar radiation extinction coefficient of a canopy of black leaves
     with an ellipsoidal leaf area distribution for beam radiation.

    Args:
        solar_zenith_angle: rad
        x_area_ratio: average area of canopy elements projected on to the horizontal
         plane divided by the average area projected on to a vertical plane

    References:
         Campbell, G.S., Norman, J.M., 1998. Introduction to environmental biophysics.
          Springer, New York. page 251. Eq. 15.4

    Examples:
        >>> get_solar_radiation_extinction_coeff_black_beam(0.087, 0)
        0.055576591963547875

        >>> get_solar_radiation_extinction_coeff_black_beam(0.087, 2)
        0.7254823957447912
    """
    assert (solar_zenith_angle >= 0 and x_area_ratio) >= 0
    assert x_area_ratio >= 0
    numerator = (x_area_ratio**2 + math.tan(solar_zenith_angle) ** 2) ** 0.5
    denominator = x_area_ratio + 1.774 * (x_area_ratio + 1.182) ** -0.733
    return numerator / denominator


def get_solar_radiation_extinction_coeff_black_diff(
    x_area_ratio: float, leaf_area_index: float
) -> float:
    """Solar radiation extinction.

    Return solar radiation extinction coefficient of a canopy of black leaves for
     diffuse radiation.

    Args:
        x_area_ratio: average area of canopy elements projected on to the horizontal
         plane divided by the average area projected on to a vertical plane
        leaf_area_index: leaf area index [m3/m3]

    References:
        Campbell, G.S., Norman, J.M., 1998. Introduction to environmental biophysics.
         Springer, New York. Eq. 15.5

    Examples:
        >>> get_solar_radiation_extinction_coeff_black_diff(2, 0.1)
        0.9710451784887358

        >>> get_solar_radiation_extinction_coeff_black_diff(0, 0.1)
        0.9099461266386164
    """
    assert (x_area_ratio and leaf_area_index) >= 0
    step_size = 90
    max_angle = math.pi / 2
    transm_diff = 0
    diff_ang = max_angle / step_size
    diff_ang_center = 0.5 * diff_ang
    angle = diff_ang
    # Integration loop
    while True:
        angle = angle - diff_ang_center
        transm_beam = math.exp(
            -get_solar_radiation_extinction_coeff_black_beam(angle, x_area_ratio)
            * leaf_area_index
        )
        transm_diff += 2 * transm_beam * math.sin(angle) * math.cos(angle) * diff_ang
        angle = angle + diff_ang_center + diff_ang
        if angle > (max_angle + diff_ang):
            break
    return -math.log(transm_diff) / leaf_area_index


def get_solar_beam_fraction(
    atm_press: float, solar_zenith_angle: float, atm_transmittance: float
) -> float:
    """Return radiation beam fraction.

    Args:
        atm_press: [kPa]
        solar_zenith_angle: [deg]
        atm_transmittance: atmospheric transmittance - 0.75 for clear sky

    Examples:
        >>> get_solar_beam_fraction(101.3, 0, 0.75)
        0.75

        >>> get_solar_beam_fraction(101.3, 50, 0.45)
        0.2892276326469122

    References:
        Campbell, G. S., and J. M. Norman. 1998. Introduction to environmental
         biophysics. Springer, New York. Ch. 11
    """
    assert 0 <= solar_zenith_angle <= SUNSET_OR_SUNRISE_SOLAR_ZENITH_ANGLE
    assert 0 <= atm_transmittance <= 1
    deg_to_rad = math.pi / 180.0
    solar_zenith_angle *= deg_to_rad
    optical_air_mass = get_optical_air_mass(atm_press, solar_zenith_angle)  # 11.12
    solar_perpend_frac = atm_transmittance**optical_air_mass  # Eq. 11.11
    return solar_perpend_frac * math.cos(solar_zenith_angle)  # 11.8


def get_solar_diffuse_fraction(
    atm_press: float, solar_zenith_angle: float, atm_transmittance: float
) -> float:
    """Return radiation diffuse fraction.

    Args:
        atm_press: [kPa]
        solar_zenith_angle: [deg]
        atm_transmittance: 0.75 for clear sky

    References:
        Campbell, G. S., and J. M. Norman. 1998. Introduction to environmental
         biophysics. Springer, New York. Ch. 11

    Examples:
        >>> get_solar_diffuse_fraction(101.3, 0, 0.75)
        0.075

        >>> get_solar_diffuse_fraction(101.3, 50, 0.45)
        0.10606799311188815
    """
    assert 0 <= solar_zenith_angle <= SUNSET_OR_SUNRISE_SOLAR_ZENITH_ANGLE
    assert 0 <= atm_transmittance <= 1

    deg_to_rad = math.pi / 180.0
    solar_zenith_angle *= deg_to_rad
    optical_air_mass = get_optical_air_mass(atm_press, solar_zenith_angle)  # 11.12
    return (
        0.3 * (1 - atm_transmittance**optical_air_mass) * math.cos(solar_zenith_angle)
    )  # 11.13


def get_solar_radiation_interception_sub_daily(
    atm_transm: float,
    atm_press: float,
    leaf_transm: float,
    leaf_area_index: NDArray,
    x_sp1: float,
    x_sp2: float,
    angles_deg: NDArray,
) -> tuple[NDArray, NDArray]:
    """Return sub daily radiation interception for two species.

    Args:
        atm_transm: atmospheric transmission [0-1]
        atm_press: atmospheric pressure
        leaf_transm: leaf transmission [0-1]
        leaf_area_index: leaf area index array
        x_sp1: average area of canopy elements projected on to the horizontal plane
         divided by the average area projected on to a vertical plane for species 1
        x_sp2: average area of canopy elements projected on to the horizontal plane
         divided by the average area projected on to a vertical plane for species 2
        angles_deg: angles range in degrees

    References:
         Campbell, G. S., and J. M. Norman. 1998. Introduction to environmental
          biophysics. Springer, New York.
    """
    deg_to_rad = math.pi / 180
    angles = angles_deg * deg_to_rad
    beam_frac = np.zeros(len(angles))
    diff_frac = np.zeros(len(angles))
    total_intercpt = np.zeros(len(angles))
    sp1_intercpt_alone = np.zeros([len(angles), len(leaf_area_index)])
    sp2_intercpt_alone = np.zeros([len(angles), len(leaf_area_index)])
    canopy_intercpt = np.zeros([len(angles), len(leaf_area_index)])
    sp1_intercpt = np.zeros([len(angles), len(leaf_area_index)])
    sp2_intercpt = np.zeros([len(angles), len(leaf_area_index)])
    sp1_intercpt_daily = np.zeros(len(leaf_area_index))
    sp2_intercpt_daily = np.zeros(len(leaf_area_index))

    for i, angle in enumerate(angles_deg):
        beam_frac[i] = get_solar_beam_fraction(atm_press, angle, atm_transm)
        diff_frac[i] = get_solar_diffuse_fraction(atm_press, angle, atm_transm)
        total_intercpt[i] = beam_frac[i] + diff_frac[i]

    for i, angle in enumerate(angles):
        for j, lai in enumerate(leaf_area_index):
            sp1_intercpt_alone[i, j] = (
                1
                - (
                    beam_frac[i]
                    / total_intercpt[i]
                    * math.exp(
                        -(leaf_transm**0.5)
                        * lai
                        * (
                            get_solar_radiation_extinction_coeff_black_beam(
                                angle, x_sp1
                            )
                        )
                    )
                    + diff_frac[i]
                    / total_intercpt[i]
                    * math.exp(
                        -(leaf_transm**0.5)
                        * lai
                        * get_solar_radiation_extinction_coeff_black_diff(x_sp1, lai)
                    )
                )
            ) * total_intercpt[i]
            sp2_intercpt_alone[i, j] = (
                1
                - (
                    beam_frac[i]
                    / total_intercpt[i]
                    * math.exp(
                        -(leaf_transm**0.5)
                        * lai
                        * (
                            get_solar_radiation_extinction_coeff_black_beam(
                                angle, x_sp2
                            )
                        )
                    )
                    + diff_frac[i]
                    / total_intercpt[i]
                    * math.exp(
                        -(leaf_transm**0.5)
                        * lai
                        * get_solar_radiation_extinction_coeff_black_diff(x_sp2, lai)
                    )
                )
            ) * total_intercpt[i]
            canopy_intercpt[i, j] = (
                1
                - (1 - sp1_intercpt_alone[i, j] / total_intercpt[i])
                * (1 - sp2_intercpt_alone[i, j] / total_intercpt[i])
            ) * total_intercpt[i]
            sp1_intercpt[i, j] = (
                -math.log(1 - sp1_intercpt_alone[i, j] / total_intercpt[i])
                / (
                    -math.log(1 - sp1_intercpt_alone[i, j] / total_intercpt[i])
                    + (-math.log(1 - sp2_intercpt_alone[i, j] / total_intercpt[i]))
                )
            ) * canopy_intercpt[i, j]
            sp2_intercpt[i, j] = (
                -math.log(1 - sp2_intercpt_alone[i, j] / total_intercpt[i])
                / (
                    -math.log(1 - sp1_intercpt_alone[i, j] / total_intercpt[i])
                    + (-math.log(1 - sp2_intercpt_alone[i, j] / total_intercpt[i]))
                )
            ) * canopy_intercpt[i, j]
    for i in range(len(leaf_area_index)):
        sp1_intercpt_daily[i] = sp1_intercpt[:, i].sum() / total_intercpt.sum()
        sp2_intercpt_daily[i] = sp2_intercpt[:, i].sum() / total_intercpt.sum()
    return sp1_intercpt_daily, sp2_intercpt_daily


def rad_intercpt_cycles(
    crop_list: tuple[list[float | int], list[float | int], list[float]],
) -> NDArray:
    """Returns solar radiation intercepted on each species.

    Args:
        crop_list: list of species characteristics with three items:
            radiation extinction coefficients
            leaf area index [m2/m2]
            plant height [m]

    References:
         Camargo, G.G.T. 2015. Ph.D. Dissertation. Penn State University.
          https://etda.libraries.psu.edu/files/final_submissions/10226

    Examples:
        >>> rad_intercpt_cycles(([0.5,1,1],[0.5,1,2],[0.5,1,1]))
        array([0.23758411, 0.30170163, 0.23758411])

        >>> rad_intercpt_cycles(([0.5,1,0.5],[0.6,1.2,1],[0.7,1.4,1.5]))
        array([0.13822214, 0.29363746, 0.45733724])
    """
    # Variables init
    number_species = len(crop_list)
    extinction_coeff = np.zeros(number_species)
    leaf_area_index = np.zeros(number_species)
    height = np.zeros(number_species)
    transm_rad = np.zeros(number_species)  # if plant was alone
    rad_intercpt_dom = np.zeros(number_species)  # dominant species
    rad_intercpt_suppr = np.zeros(number_species)  # supressed species
    rad_intercpt = np.zeros(number_species)  # species rad interception
    transm_rad_others = np.ones(number_species)  # non-domnt sp rad transm
    height_dom = np.zeros(number_species)  # species canopy dominance factor
    dominant_fact = np.zeros(number_species)  # Dominant weight factor
    suppressed_fact = np.zeros(number_species)  # Suppressed weight factor
    hght_wght_fct = np.zeros(number_species)
    hght_wght_fct_adj = np.zeros(number_species)
    k_lai_prod = np.zeros(number_species)
    k_lai_prod_adj = np.zeros(number_species)
    total_transm = 1

    # Read and store inputs
    for i in range(number_species):
        extinction_coeff[i] = crop_list[i][0]
        leaf_area_index[i] = crop_list[i][1]
        height[i] = crop_list[i][2]

        # Transmitted radiation if all species had same height
        transm_rad[i] = math.exp(-extinction_coeff[i] * leaf_area_index[i])

        # Intercepted radiation if species was dominant
        rad_intercpt_dom[i] = 1 - transm_rad[i]
        k_lai_prod[i] = extinction_coeff[i] * leaf_area_index[i]

    # Calculate total transmitance, interception and height dominance
    for i in range(number_species):
        # Total transmitance if all species had the same height
        total_transm *= transm_rad[i]

        # Height dominance factor
        height_dom[i] = number_species * height[i] / height.sum()

    # Total radiation interception if species had the same height
    total_interception = 1 - total_transm

    # All species but ith species rad transmission
    for i in range(number_species):
        all_species = list(range(number_species))  # list with all species
        all_species.pop(i)  # remove ith species from list
        all_species_but_ith = all_species  # list of non-dominant species
        total_transm_but_ith_sp = 1  # sum of non-dominant species transmission
        for j in all_species_but_ith:
            total_transm_but_ith_sp *= transm_rad[j]

        # Total transmitted radiation from all species but ith
        transm_rad_others[i] = total_transm_but_ith_sp

    # Radiation interception by suppressed species once all other species
    # intercepts the radiation first
    for i in range(number_species):
        rad_intercpt_suppr[i] = rad_intercpt_dom[i] * transm_rad_others[i]

    # Determine two extremes weighing factors: in dominant_fact species will
    # intercept all radiation than it can based on k and LAI, in
    # suppressed_factor, species will only intercept radiation after all the
    # other species intercepted all rad that was possible
    for i in range(number_species):
        dominant_fact[i] = (
            rad_intercpt_dom[i] / total_interception * k_lai_prod.sum() / k_lai_prod[i]
        )
        suppressed_fact[i] = (
            rad_intercpt_suppr[i]
            / total_interception
            * k_lai_prod.sum()
            / k_lai_prod[i]
        )

        # Based on species height determine a height weight factor in between
        # dominant_fact and suppressed_fact values usin linear interpolation
        hght_wght_fct[i] = get_height_weight_factor(
            height_dom[i], dominant_fact[i], suppressed_fact[i], number_species
        )
        # Adjust extinction coefficient and leaf area index product
        k_lai_prod_adj[i] = extinction_coeff[i] * leaf_area_index[i] * hght_wght_fct[i]

    for i in range(number_species):
        # Adjust height weighting factor
        hght_wght_fct_adj[i] = (
            hght_wght_fct[i] / k_lai_prod_adj.sum() * k_lai_prod.sum()
        )

        # Radiation interception for each species
        rad_intercpt[i] = (
            total_interception * hght_wght_fct_adj[i] * k_lai_prod[i] / k_lai_prod.sum()
        )
    return rad_intercpt
