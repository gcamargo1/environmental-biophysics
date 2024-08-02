"""Water vapor and other gases functions."""


def get_vapor_pressure_air(
    vapor_pressure_temp_min: float,
    vapor_pressure_temp_max: float,
    rh_max: float,
    rh_min: float,
) -> float:
    """Return the vapor pressure of air (kPa).

    Args:
        vapor_pressure_temp_min: saturated vapor pressure of min temperature (kPa)
        vapor_pressure_temp_max: saturated vapor pressure of max temperature (kPa)
        rh_max: max relative humidity (0 - 100)
        rh_min: min relative humidity (0 - 100)

    References:
        Campbell, G.S., Norman, J.M., 1998. Introduction to environmental biophysics.
         Springer, New York.

    Returns:
        Vapor pressure of air (kPa).

    Examples:
        >>> get_vapor_pressure_air(1.817, 5.320, 87, 25)
        1.455
    """
    result = 0.5 * (
        vapor_pressure_temp_min * rh_max / 100.0
        + vapor_pressure_temp_max * rh_min / 100.0
    )
    return round(result, 3)


def get_vapor_press_defct_ave(
    max_sat_vap_press: float, min_sat_vap_press: float, air_vap_press: float
) -> float:
    """Returns Average vapor pressure deficit.

    Args:
        max_sat_vap_press:
        min_sat_vap_press:
        air_vap_press:
    """
    return (max_sat_vap_press + min_sat_vap_press) / 2.0 - air_vap_press


def get_vapor_press_defct_max(
    max_sat_vap_press: float, min_relative_humidity: float
) -> float:
    """Max vapor presure deficit."""
    return 0.67 * max_sat_vap_press * (1 - min_relative_humidity / 100.0)
