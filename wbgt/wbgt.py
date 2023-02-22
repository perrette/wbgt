from pathlib import Path
import ctypes
from ctypes import c_int, c_double, byref, POINTER
import numpy as np
from numpy.ctypeslib import ndpointer

import logging
logging.basicConfig()
logger = logging.getLogger('wbgt')
logger.setLevel(logging.INFO)


# import shared lib with ctypes
ext = ctypes.CDLL(list(Path(__file__).parent.glob('*.so'))[0])

ext.Twb.restype = c_double
ext.Twb.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int]
def Twb(tk, rh, pres, speed, solar, fdir, cza, rad):
    return ext.Twb(tk, rh, pres, speed, solar, fdir, cza, rad)

ext.Tglobe.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double]
ext.Tglobe.restype = c_double
def Tglobe(tk, rh, Pair, speed, solar, fdir, cza):
    return ext.Tglobe(tk, rh, Pair, speed, solar, fdir, cza)


ext.stab_srdt.argypes = [c_int, c_double, c_double, c_double]
ext.stab_srdt.restype = c_int
def stab_srdt(daytime, speed, solar, dT):
    return ext.stab_srdt(daytime, speed, solar, dT)


ext.est_wind_speed.argtypes = [c_double, c_double, c_int, c_int]
ext.est_wind_speed.restype = c_double
def est_wind_speed(speed, zspeed, stability_class, urban):
    return ext.est_wind_speed(speed, zspeed, stability_class, urban)


def estimate_windspeed(zspeed, speed, solar, dT, cza, urban):
    REF_HEIGHT=2

    if ( zspeed == REF_HEIGHT ):
        return speed

    if ( cza > 0.):
        daytime = 1
    else:
        daytime = 0

    stability_class = stab_srdt(daytime, speed, solar, dT);
    return est_wind_speed(speed, zspeed, stability_class, urban);


ext.calc_solar_parameters.argtypes = [c_int, c_int, c_double, c_double, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
def calc_solar_parameters(year, month, dday, lat, lon, solar):
    solar = c_double(solar)
    cza = c_double()
    fdir = c_double()
    ext.calc_solar_parameters(year, month, dday, lat, lon, byref(solar), byref(cza), byref(fdir))
    return solar.value, cza.value, fdir.value


ext.solarposition.argtypes = [c_int, c_int, c_double, c_double, c_double,
                  c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double),
                  POINTER(c_double), POINTER(c_double), POINTER(c_double)]
def solarposition(year, month, day, days_1900, latitude, longitude):
    ap_ra = c_double()
    ap_dec = c_double()
    altitude = c_double()
    refraction = c_double()
    azimuth = c_double()
    distance = c_double()
    ext.solarposition(year, month, day, days_1900, latitude, longitude,
        byref(ap_ra), byref(ap_dec), byref(altitude), byref(refraction), byref(azimuth), byref(distance))
    return ap_ra.value, ap_dec.value, altitude.value, refraction.value, azimuth.value, distance.value


def wbgt(tk, rh, pres, speed, solar, fdir, cza):
    """
    tk: temperature in Kelvin
    rh: relative humidity as fraction
    pres: pressure in mbar (typical value ~ 1000 mb)
    solar: sun input power at a given hour, lon, lat... (as returned by calc_solar_parameters in the original code)
    fdir: fraction of direct sun input
    cza : cosinus azimuh angle

    returns Tg, Tnwb, WBGT
    """
    # calculate the globe, natural wet bulb, psychrometric wet bulb, and
    # outdoor wet bulb globe temperatures
    Tg   = Tglobe(tk, rh, pres, speed, solar, fdir, cza)
    Tnwb = Twb(tk, rh, pres, speed, solar, fdir, cza, 1)
    # Tpsy = Twb(tk, rh, pres, speed, solar, fdir, cza, 0)
    # return Tg, Tnwb, Tpsy, 0.1 * (tk-273.15) + 0.2 * Tg + 0.7 * Tnwb
    Twbg = 0.1 * (tk-273.15) + 0.2 * Tg + 0.7 * Tnwb
    return Tg, Tnwb, Twbg


# main function from original code
ext.calc_wbgt.argtypes = [c_int, c_int, c_int, c_int, c_int, c_int, c_int,
              c_double, c_double, c_double, c_double, c_double, c_double,
              c_double, c_double, c_double, c_int, POINTER(c_double),
                  POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
def calc_wbgt(year, month, day, hour, minute, gmt, avg,
              lat, lon, solar, pres, Tair, relhum,
              speed, zspeed, dT, urban):
    """
    returns:
        est_speed
        Tg
        Tnwb
        Tpsy
        Twbg
    """
    est_speed = c_double()
    Tg = c_double()
    Tnwb = c_double()
    Tpsy = c_double()
    Twbg = c_double()

    return_value = ext.calc_wbgt(year, month, day, hour, minute, gmt, avg,
              lat, lon, solar, pres, Tair, relhum,
              speed, zspeed, dT, urban, byref(est_speed),
                byref(Tg), byref(Tnwb), byref(Tpsy), byref(Twbg))

    return est_speed.value, Tg.value, Tnwb.value, Tpsy.value, Twbg.value


ext.calc_wbgt_vector.argtypes = [ctypes.c_size_t,
              ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"),
              ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"),
              ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"),
                  ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_double, flags="C_CONTIGUOUS"), ndpointer(c_int, flags="C_CONTIGUOUS")]
def calc_wbgt_vector(year, month, day, hour, minute, gmt, avg,
              lat, lon, solar, pres, Tair, relhum,
              speed, zspeed, dT, urban):
    """
    returns:
        est_speed
        Tg
        Tnwb
        Tpsy
        Twbg
    """
    # Broadcast all inputs to the same dimensions
    year, month, day, hour, minute, gmt, avg, lat, lon, solar, pres, Tair, relhum, speed, zspeed, dT, urban = (np.array(a) for a in np.broadcast_arrays(*locals().values()))
    size = year.size
    shape = year.shape
    # year, month, day, hour, minute, gmt, avg, lat, lon, solar, pres, Tair, relhum, speed, zspeed, dT, urban = [np.array(a) for a in np.broadcast_arrays(2, np.ones(3))]

    if ((year > 2049) | (year < 1950)).any():
        raise ValueError("Solar calculations only defined between 1950 and 2049.")


    est_speed = np.empty(shape)
    Tg = np.empty(shape)
    Tnwb = np.empty(shape)
    Tpsy = np.empty(shape)
    Twbg = np.empty(shape)
    status = np.zeros(shape, dtype='int32')

    return_value = ext.calc_wbgt_vector(year.size, year.astype('int32'), month.astype('int32'), day.astype('int32'),
        hour.astype('int32'), minute.astype('int32'), gmt.astype('int32'), avg.astype('int32'),
        lat.astype('float64'), lon.astype('float64'), solar.astype('float64'), pres.astype('float64'), Tair.astype('float64'), relhum.astype('float64'),
        speed.astype('float64'), zspeed.astype('float64'), dT.astype('float64'), urban.astype('int32'), est_speed, Tg, Tnwb, Tpsy, Twbg, status)

    if (status != 0).any():
        logger.warning("Some values did not converge or are incorrect. Set with nan.")
        bad = status != 0
        est_speed[bad] = np.nan
        Tg[bad] = np.nan
        Tnwb[bad] = np.nan
        Tpsy[bad] = np.nan
        Twbg[bad] = np.nan

    return est_speed, Tg, Tnwb, Tpsy, Twbg