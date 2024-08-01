def vapor_pressure_air(
    vapor_pressure_temp_min: float,
    vapor_pressure_temp_max: float,
    rh_max: float,
    rh_min: float,
) -> float:
    """Return the vapor pressure of air (kPa).

    vapor_pressure_temp_min: saturated vapor pressure of min temperature (kPa)
    vapor_pressure_temp_max: saturated vapor pressure of max temperature (kPa)
    rh_max: max relative humidity (0 - 100)
    rh_min: min relative humidity (0 - 100)

    Reference: Campbell, G.S., Norman, J.M., 1998. Introduction to environmental
     biophysics. Springer, New York.

    >>> vapor_pressure_air(1.817,5.320,87,25)
    1.455
    """
    return 0.5 * (
        vapor_pressure_temp_min * rh_max / 100.0
        + vapor_pressure_temp_max * rh_min / 100.0
    )


def vapor_press_defct_ave(
    max_sat_vap_press: float, min_sat_vap_press: float, air_vap_press: float
) -> float:
    """max_sat_vap_press: max sat vapor pressure - kPa
    min_sat_vap_press: min sat vapor pressure - kPa
    air_vap_press: air vapor pressure - kPa
    """
    return (max_sat_vap_press + min_sat_vap_press) / 2.0 - air_vap_press


def vapor_press_defct_max(
    max_sat_vap_press: float, min_relative_humidity: float
) -> float:
    """Max vapor presure deficit."""
    return 0.67 * max_sat_vap_press * (1 - min_relative_humidity / 100.0)
