# eos80.py

"""
EOS-80 sigma-theta (UNESCO 1981 equation of state at zero pressure, Bryden
1973 / UNESCO 1983 potential temperature), the definition behind the CalCOFI
bottle database's STheta, so CTD and bottle densities are computed alike.

    sigma_theta(t, s, p) -> kg m-3 minus 1000
    t: in-situ temperature (degC, IPTS-68 as in EOS-80; ITS-90 values are
       converted with t68 = 1.00024 t90), s: practical salinity, p: dbar.
"""

import numpy as np


def _rho0(t, s):
    """Density at zero pressure (UNESCO 1981), kg m-3."""
    rw = (999.842594 + 6.793952e-2 * t - 9.095290e-3 * t**2 + 1.001685e-4 * t**3
          - 1.120083e-6 * t**4 + 6.536332e-9 * t**5)
    a = 8.24493e-1 - 4.0899e-3 * t + 7.6438e-5 * t**2 - 8.2467e-7 * t**3 + 5.3875e-9 * t**4
    b = -5.72466e-3 + 1.0227e-4 * t - 1.6546e-6 * t**2
    return rw + a * s + b * s**1.5 + 4.8314e-4 * s**2


def _adtg(s, t, p):
    """Adiabatic lapse rate (Bryden 1973), degC per dbar."""
    ds = s - 35.0
    return ((((-2.1687e-16 * t + 1.8676e-14) * t - 4.6206e-13) * p
             + ((2.7759e-12 * t - 1.1351e-10) * ds + ((-5.4481e-14 * t + 8.733e-12) * t - 6.7795e-10) * t + 1.8741e-8)) * p
            + (-4.2393e-8 * t + 1.8932e-6) * ds + ((6.6228e-10 * t - 6.836e-8) * t + 8.5258e-6) * t + 3.5803e-5)


def potential_temperature(s, t, p, pr=0.0):
    """Potential temperature (UNESCO 1983, 4th-order Runge-Kutta), degC."""
    h = pr - p
    xk = h * _adtg(s, t, p)
    t = t + 0.5 * xk
    q = xk
    p = p + 0.5 * h
    xk = h * _adtg(s, t, p)
    t = t + 0.29289322 * (xk - q)
    q = 0.58578644 * xk + 0.121320344 * q
    xk = h * _adtg(s, t, p)
    t = t + 1.707106781 * (xk - q)
    q = 3.414213562 * xk - 4.121320344 * q
    p = p + 0.5 * h
    xk = h * _adtg(s, t, p)
    return t + (xk - 2.0 * q) / 6.0


def sigma_theta(t, s, p, its90=True):
    t, s, p = (np.asarray(x, dtype=float) for x in (t, s, p))
    t68 = t * 1.00024 if its90 else t
    return _rho0(potential_temperature(s, t68, p), s) - 1000.0


def pressure_from_depth(z, lat):
    """Saunders (1981) pressure (dbar) from depth (m)."""
    x = np.sin(np.radians(np.abs(np.asarray(lat, dtype=float))))
    c1 = 5.92e-3 + x**2 * 5.25e-3
    return ((1 - c1) - np.sqrt((1 - c1) ** 2 - 8.84e-6 * np.asarray(z, dtype=float))) / 4.42e-6
