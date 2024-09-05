"""Soil physics related functions."""

import math

from numpy.typing import NDArray

# Constants
MIN_SOIL_PARTICLE_DENS = 2.65  # Mg/m3
FIELD_CAPACITY_WATER_POT = -33  # J/kg
PERMNT_WILT_POINT_WATER_POT = -1500  # J/kg
WATER_DENSITY = 1000  # kg/m3


def get_bulk_density(
    clay: float | NDArray, sand: float | NDArray, organic_matter: float | NDArray
) -> float:
    """Returns bulk density (Mg/m3 or g/cm3)

    Args:
        clay: clay content (fraction)
        sand: sand content (fraction)
        organic_matter: organic matter (%)

    References:
        Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic
        estimates by texture and organic matter for hydrologic solutions. Eq. 5,6
        Soil Sci. Soc. Am. J. 70, 1569-1578.

    Examples:
        >>> get_bulk_density(clay=0.03, sand=0.92, organic_matter=1.906)
        1.43
        >>> get_bulk_density(clay=0.15, sand=0.2, organic_matter=2.29)
        1.39
    """
    x1 = (
        0.078
        + 0.278 * sand
        + 0.034 * clay
        + 0.022 * organic_matter
        - 0.018 * sand * organic_matter
        - 0.027 * clay * organic_matter
        - 0.584 * sand * clay
    )
    x2 = -0.107 + 1.636 * x1
    field_capacity = get_vol_water_content_33_j_kg(clay, sand, organic_matter)  # m3/m3
    sat_water_content = 0.043 + field_capacity + x2 - 0.097 * sand
    result = (1 - sat_water_content) * MIN_SOIL_PARTICLE_DENS
    return round(result, 2)  # type: ignore


def get_sat_water_content(bulk_density: float | NDArray) -> float:
    """Returns the saturated water content (m3/m3).

    Args:
        bulk_density: bulk density (Mg / m3)

    References:
        Campbell, G.S., 1985. Soil physics with BASIC: Transport models for soil-plant
         systems. Elsevier, Amsterdam.

    Examples:
        >>> get_sat_water_content(bulk_density=1.3)
        0.509
    """
    result = 1 - bulk_density / MIN_SOIL_PARTICLE_DENS
    return round(result, 3)  # type: ignore


def get_vol_water_content_33_j_kg(
    clay: float | NDArray, sand: float | NDArray, organic_matter: float | NDArray
) -> float:
    """Returns the volumetric water content at field capacity (33 J/kg) (m3/m3).

    Args:
        clay: clay content (fraction)
        sand: sand content (fraction)
        organic_matter: organic matter (%)

    References:
        Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic estimates by texture
         and organic matter for hydrologic solutions. Soil Sci. Soc. Am. J. 70,
         1569-1578. eq.2 R2=0.63.

    Examples:
        >>> get_vol_water_content_33_j_kg (clay=0.03, sand=0.92, organic_matter=1.906)
        0.0764
        >>> get_vol_water_content_33_j_kg (clay=0.33, sand=0.09, organic_matter=2.866)
        0.3829
    """
    x1 = (
        0.299
        - 0.251 * sand
        + 0.195 * clay
        + 0.011 * organic_matter
        + 0.006 * sand * organic_matter
        - 0.027 * clay * organic_matter
        + 0.452 * sand * clay
    )
    result = -0.015 + 0.636 * x1 + 1.283 * x1**2
    return round(result, 4)  # type: ignore


def get_vol_water_content_1500_jkg(
    clay: float | NDArray, sand: float | NDArray, organic_matter: float | NDArray
) -> float:
    """Returns the volumetric water content at permanent wilting point (1500 J/kg).

    Args:
        clay: clay content (fraction)
        sand: sand content (fraction)
        organic_matter: organic matter (%)

    References:
        Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic estimates by texture
         and organic matter for hydrologic solutions. Soil Sci. Soc. Am. J. 70,
         1569-1578. eq.1 R2=0.86

    Returns:
        Wilting point volumetric water content (m3/m3).

    Examples:
        >>> get_vol_water_content_1500_jkg (0.03,0.92,1.906)
        0.03

        >>> get_vol_water_content_1500_jkg (0.33,0.09,2.866)
        0.21
    """
    x1 = (
        0.031
        - 0.024 * sand
        + 0.487 * clay
        + 0.006 * organic_matter
        + 0.005 * sand * organic_matter
        - 0.013 * clay * organic_matter
        + 0.068 * sand * clay
    )
    result = -0.02 + 1.14 * x1
    return round(result, 2)  # type: ignore


def get_b_value(
    water_content_33_j_kg: float | NDArray, water_content_1500_j_kg: float | NDArray
) -> float:
    """Return b soil parameter.

    Args:
        water_content_33_j_kg: water content at -33 J/kg
        water_content_1500_j_kg: water content at -1500 J/kg

    References:
        Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic estimates by texture
         and organic matter for hydrologic solutions. Soil Sci. Soc. Am. J. 70,
         1569-1578.

    Examples:
        >>> get_b_value(0.08, 0.03)
        3.89
    """
    result = (math.log(1500) - math.log(33)) / (
        math.log(water_content_33_j_kg) - math.log(water_content_1500_j_kg)
    )
    return round(result, 2)


def get_air_entry_pot(
    field_capacity: float | NDArray,
    sat_water_content: float | NDArray,
    b_value: float | NDArray,
) -> float:
    """Return air entry potential.

    Args:
        field_capacity: water content at field capacity (m3/m3)
        sat_water_content: saturated water content
        b_value: soil parameter

    References:
        Kemanian, A.R., Stockle, C.O., 2010. C-Farm: A simple model to evaluate the
         carbon balance of soil profiles. Eur. J. Agron. 32, 22-29.

    Examples:
        >>> get_air_entry_pot(.08,0.5,4.33)
        -0.0118
    """
    result = -33 * (field_capacity / sat_water_content) ** b_value
    return round(result, 4)  # type: ignore


def get_water_potential(
    sat_water_content: float | NDArray,
    air_entry_potential: float | NDArray,
    campbell_b: float | NDArray,
    water_content: float | NDArray,
) -> float:
    """Returns Soil Water Potential (J/kg).

    Args:
        sat_water_content: saturation water content (m3/m3)
        air_entry_potential: air entry water potential (J/kg)
        campbell_b: soil moisture release curve parameter
        water_content: water content (m3/m3)

    References:
        Campbell, G.S., 1985. Soil physics with BASIC: Transport models for soil-plant
         systems. Elsevier, Amsterdam. Eq. 5.9

    Examples:
        >>> get_water_potential (0.5, -1.5, 5, 0.25)
        -48.0

        >>> get_water_potential (0.20, -1.0, 4, 0.25)
        -0.4096
    """
    assert sat_water_content > 0, "sat water content must be positive"
    assert water_content > 0, "water content must be positive"

    return round(
        air_entry_potential * (sat_water_content / water_content) ** campbell_b,  # type: ignore
        4,
    )


def get_water_content(
    sat_water_content: float | NDArray,
    air_entry_potential: float | NDArray,
    campbell_b: float | NDArray,
    water_potential: float | NDArray,
) -> float:
    """Returns Soil water content (m3/m3).

    Args:
        sat_water_content: saturation water content (m3/m3)
        air_entry_potential: air entry water potential (J/kg)
        campbell_b: soil moisture release curve parameter
        water_potential: matric water potential (J/kg)

    References:
        Campbell, G.S., 1985. Soil physics with BASIC: Transport models for soil-plant
         systems. Elsevier, Amsterdam. pp.80

    Examples:
        >>> get_water_content(0.5,-1.5,5,-52.7)
        0.2454
    """
    result = sat_water_content * (water_potential / air_entry_potential) ** (
        -1 / campbell_b
    )
    return round(result, 4)  # type: ignore


def get_organic_matter(clay: float | NDArray) -> float:
    """Half of carbon saturation given by the original author, and converted to organic
     matter carbon = 0.58 * organic_matter.

    Args:
        clay: clay content (0-1)

    References:
        Hassink, J., and A. P. Whitmore. 1997. A Model of the Physical Protection of
         Organic Matter in Soils. Soil Sci. Soc. Am. J. 61(1):131-139.

    Examples:
        >>> get_organic_matter(0.5)
        3.41

        >>> get_organic_matter(0.03)
        1.91
    """
    result = 1.81 + 0.032 * clay * 100
    return round(result, 2)  # type: ignore


def get_evap_lim_water_content(
    field_capacity_water_content: float | NDArray,
    air_dry_water_content: float | NDArray,
    water_content: float | NDArray,
) -> float:
    """Returns evaporation limitation based on water content [m3/m3].

    Args:
        field_capacity_water_content: field_capacity water content [m3/m3]
        air_dry_water_content: air dry water content [m3/m3]
        water_content: water content [m3/m3]

    References:
        Kemanian, A. R., and C. O. Stockle. 2010. C-Farm: A simple model to evaluate
         the carbon balance of soil profiles. Eur. J. Agron. 32(1):22-29

    Examples:
        >>> get_evap_lim_water_content(0.35, 0.05, 0.25)
        0.296

        >>> get_evap_lim_water_content(0.25, 0.05, 0.15)
        0.125
    """
    assert (
        field_capacity_water_content and air_dry_water_content and water_content
    ) >= 0
    assert (
        field_capacity_water_content and air_dry_water_content and water_content
    ) <= 1
    assert field_capacity_water_content > air_dry_water_content
    assert field_capacity_water_content >= water_content
    result = (
        (water_content - air_dry_water_content)
        / (field_capacity_water_content - air_dry_water_content)
    ) ** 3
    return round(result, 3)  # type: ignore
