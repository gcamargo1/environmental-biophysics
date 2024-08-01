import math

# Constants
MIN_SOIL_PARTICLE_DENS = 2.65  # Mg/m3


def get_bulk_density(clay: float, sand: float, organic_matter: float) -> float:
    """Returns Bulk density (Mg/m3 or g/cm3)

    clay: clay content (fraction)
    sand: sand content (fraction)
    organic_matter: organic matter (%)

    Reference: Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic
     estimates by texture and organic matter for hydrologic solutions. Eq. 5,6
     Soil Sci. Soc. Am. J. 70, 1569-1578.

    >>> get_bulk_density(0.03,0.92,1.906)
    1.43
    >>> get_bulk_density(0.15,0.2,2.29)
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
    field_capacity = vol_water_content_33_j_kg(clay, sand, organic_matter)  # m3/m3
    sat_water_content = 0.043 + field_capacity + x2 - 0.097 * sand
    return (1 - sat_water_content) * MIN_SOIL_PARTICLE_DENS


def sat_water_content(bulk_density: float) -> float:
    """Returns the saturated water content (m3/m3).

    bulk_density: bulk density (Mg / m3)

    Reference: Campbell, G.S., 1985. Soil physics with BASIC: Transport models
     for soil-plant systems. Elsevier, Amsterdam.

     >>> sat_water_content(1.3)
    0.5094339622641508
    """
    return 1 - bulk_density / MIN_SOIL_PARTICLE_DENS


def vol_water_content_33_j_kg(clay: float, sand: float, organic_matter: float) -> float:
    """Returns the volumetric water content at field capacity (33 J/kg) (m3/m3).

    clay: clay content (fraction)
    sand: sand content (fraction)
    organic_matter: organic matter (%)

    Reference: Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic
     estimates by texture and organic matter for hydrologic solutions.
     Soil Sci. Soc. Am. J. 70, 1569-1578. eq.2 R2=0.63

    >>> vol_water_content_33_j_kg (0.03,0.92,1.906)
    0.08
    >>> vol_water_content_33_j_kg (0.33,0.09,2.866)
    0.38
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
    return -0.015 + 0.636 * x1 + 1.283 * x1**2


def vol_water_content_1500_jkg(
    clay: float, sand: float, organic_matter: float
) -> float:
    """Returns the volumetric water content at field capacity (33 J/kg) (m3/m3)

    clay: clay content (fraction)
    sand: sand content (fraction)
    organic_matter: organic matter (%)

    Reference: Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic
     estimates by texture and organic matter for hydrologic solutions.
     Soil Sci. Soc. Am. J. 70, 1569-1578. eq.1 R2=0.86

    >>> vol_water_content_1500_jkg (0.03,0.92,1.906)
    0.03
    >>> vol_water_content_1500_jkg (0.33,0.09,2.866)
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
    return -0.02 + 1.14 * x1


def b_value(water_content_33_j_kg: float, water_content_1500_j_kg: float) -> float:
    """Return b soil parameter.

    water_content_33_j_kg: water content at -33 J/kg
    water_content_1500_j_kg: water content at -1500 J/kg

    Reference: Saxton, K.E., Rawls, W.J., 2006. Soil water characteristic
     estimates by texture and organic matter for hydrologic solutions. Soil Sci.
     Soc. Am. J. 70, 1569-1578.

    >>> b_value(0.08, 0.03)
    3.89
    """
    return (math.log(1500) - math.log(33)) / (
        math.log(water_content_33_j_kg) - math.log(water_content_1500_j_kg)
    )


def air_entry_pot(
    field_capacity: float, sat_water_content: float, b_value: float
) -> float:
    """Return air entry potential.

    field_capacity: water content at field capacity
    sat_water_content: saturated water content
    b_value: soil parameter

    Reference: Kemanian, A.R., Stockle, C.O., 2010. C-Farm: A simple model to
    evaluate the carbon balance of soil profiles. Eur. J. Agron. 32, 22-29.
    >>> air_entry_pot(.08,0.5,4.33)
    -0.0118
    """
    return -33 * (field_capacity / sat_water_content) ** b_value


def water_potential(
    sat_water_content: float,
    air_entry_potential: float,
    campbell_b: float,
    water_content: float,
) -> float:
    """Returns Soil Water Potential (J/kg)

    sat_water_content : saturation water content (m3/m3)
    air_entry_potential: air entry water potential (J/kg)
    campbell_b: soil moisture release curve parameter
    water_content: water content (m3/m3)

    Reference: Campbell, G.S., 1985. Soil physics with BASIC: Transport models
     for soil-plant systems. Elsevier, Amsterdam. Eq. 5.9

    >>> water_potential (0.5, -1.5, 5, 0.25)
    -48.0
    >>> water_potential (0.20, -1.0, 4, 0.25)
    -0.4096"""
    assert sat_water_content > 0, "sat water content must be positive"
    assert water_content > 0, "water content must be positive"

    return air_entry_potential * (sat_water_content / water_content) ** campbell_b


def water_content(
    sat_water_content: float,
    air_entry_potential: float,
    campbell_b: float,
    water_potential: float,
) -> float:
    """Returns Soil water content (m3/m3)

    sat_water_content : saturation water content (m3/m3)
    air_entry_potential: air entry water potential (J/kg)
    campbell_b: soil moisture release curve parameter
    water_potential: matric water potential (J/kg)

    Reference: Campbell, G.S., 1985. Soil physics with BASIC: Transport models
     for soil-plant systems. Elsevier, Amsterdam. pp.80

    >>> water_content(0.5,-1.5,5,-52.7)
    0.24
    """
    return sat_water_content * (water_potential / air_entry_potential) ** (
        -1 / campbell_b
    )


def organic_m(clay: float) -> float:
    """Half of carbon saturation given by the original author, and converted to organic
     matter carbon = 0.58 * organic_matter

    clay: clay content (0-1)

    Hassink, J., and A. P. Whitmore. 1997. A Model of the Physical Protection of
     Organic Matter in Soils. Soil Sci. Soc. Am. J. 61(1):131-139.

    >>> organic_m(0.5)
    3.41
    >>> organic_m(0.03)
    1.91
    """
    return 1.81 + 0.032 * clay * 100
