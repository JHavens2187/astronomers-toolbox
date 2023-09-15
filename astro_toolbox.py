# a toolbox for astronomers
# j. Havens, 2023

import datetime
import math
from datetime import datetime, timedelta

import astroplan
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astroplan import AltitudeConstraint, AirmassConstraint, MoonSeparationConstraint
from astroplan import Observer, FixedTarget
from astroplan import plots
from astropy.coordinates import SkyCoord, FK5
from astropy.time import Time
from astroquery.sdss import SDSS
from astroquery.simbad import Simbad

import LST_calculator as lst


# Utility Functions
def deg2rad(deg):
    return deg * np.pi / 180


def rad2deg(rad):
    return rad * 180 / np.pi


def deg2hms(deg):
    hr = deg / 15
    min = (hr - math.floor(hr)) * 60
    sec = (min - math.floor(min)) * 60
    return math.floor(hr), math.floor(min), sec


def hms2deg(hr, min, sec):
    return (hr + (min / 60) + (sec / 3600)) * 15


def dms2deg(deg, min, sec):
    return deg + (min / 60) + (sec / 3600)


def deg2dms(decdeg):
    deg = math.floor(decdeg)
    decmin = (decdeg - deg) * 60
    decsec = (decmin - math.floor(decmin)) * 60
    return deg, math.floor(decmin), decsec


def dms2hms(deg, min, sec):
    decdeg = degminsec2deg(deg, min, sec)
    dechr = 15 * decdeg
    decmin = (dechr - math.floor(dechr)) * 60
    decsec = (decmin - math.floor(decmin)) * 60
    return math.floor(dechr), math.floor(decmin), decsec


def hms2dms(hr, min, sec):
    decdeg = hms2deg(hr, min, sec)
    deg = math.floor(decdeg)
    decmin = (decdeg - deg) * 60
    decsec = (decmin - math.floor(decmin)) * 60
    return deg, math.floor(decmin), decsec


def LST(lat, lon):
    """
    calculates the local sidereal time
    :return: LST_hours, LST_minutes, LST_seconds
    """
    # define your location
    if lat and lon is None:  # assume Lawrence, KS
        latitude = 38.9717  # Your latitude in decimal degrees (positive for North, negative for South)
        longitude = -95.2353  # Your longitude in decimal degrees (positive for East, negative for West)
    else:
        latitude = lat
        longitude = lon
    return lst.LST(latitude, longitude)


# Coordinate Transformations
def equatorial_to_galactic(ra, dec):
    """
    converts equatorial coordinates to galactic coordinates
    :param ra: ra in degrees
    :param dec: dec in degrees
    :return: l and b IN DEGREES
    """
    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    galactic = c.galactic
    return galactic.l.degree, galactic.b.degree


def galactic_to_equatorial(l, b):
    """
    converts galactic coordinates to equatorial coordinates
    :param l: l in degrees
    :param b: b in degrees
    :return: ra and dec IN DEGREES
    """
    c = SkyCoord(l=l * u.degree, b=b * u.degree, frame='galactic')
    equatorial = c.icrs
    return equatorial.ra.degree, equatorial.dec.degree


# Time and Motion
def proper_motion(ra, dec, pm_ra, pm_dec, epoch='2000-01-01'):
    """
    calculates the proper motion of an object
    :param ra: ra in degrees
    :param dec: dec in degrees
    :param pm_ra: proper motion in ra in deg/yr
    :param pm_dec: proper motion in dec in deg/yr
    :param epoch: epoch of reference, default to J2000
    :return: current position as a result of proper motion IN DEGREES
    """
    current_time = Time(datetime.now())
    delta_time = current_time - Time(epoch)
    new_ra = ra + pm_ra * delta_time.to(u.yr).value
    new_dec = dec + pm_dec * delta_time.to(u.yr).value
    return new_ra, new_dec


def precession(ra, dec, epoch="2000-01-01"):
    """
    calculates the precession of an object
    :param ra: ra in degrees
    :param dec: dec in degrees
    :param epoch: epoch of reference, default to J2000
    :return: updated ra and dec IN DEGREES
    """
    epoch_time = Time(epoch)
    current_time = Time(datetime.now())
    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame=FK5, equinox=epoch_time)
    c_prec = c.transform_to(FK5(equinox=current_time))
    return c_prec.ra.degree, c_prec.dec.degree


# Celestial Object Information from SIMBAD
def get_object_coordinates(object_name):
    """
    gets the coordinates of an object from SIMBAD
    :param object_name: the name of the object
    :return: ra and dec of object IN HMS/DMS
    """

    result_table = Simbad.query_object(object_name)
    coord = SkyCoord(ra=result_table['RA'].data[0], dec=result_table['DEC'].data[0], unit=(u.deg, u.deg))
    ra = deg2hms(coord.ra.degree)
    dec = deg2dms(coord.dec.degree)
    return ra, dec


# Observational Planning
def observable_time(observer_location, object_name):
    """
    calculates if an object is observable in the next 24 hours
    :param observer_location: location of observer
    :param object_name: name of the object
    :return: boolean if object is observable
    """
    observer = Observer(location=observer_location)
    target = FixedTarget.from_name(object_name)
    # Get the current date and time in UTC format
    now = datetime.utcnow()
    # Set the time range for the next 24 hours
    start_time = Time(now)
    end_time = Time(now + timedelta(days=1))
    time_range = Time([start_time, end_time])
    constraints = [AltitudeConstraint(min=30 * u.deg), AirmassConstraint(max=2),
                   MoonSeparationConstraint(min=30 * u.deg)]
    observable = astroplan.is_observable(constraints, observer, target, time_range=time_range)
    return observable


def airmass_plot(observer_location, object_name):
    """
    plots the airmass of an object over the next 24 hours
    :param observer_location: observer location
    :param object_name: object name
    :return: plot of the airmass of the object over the next 24 hours
    """
    observer = Observer(location=observer_location)
    target = FixedTarget.from_name(object_name)
    # Get the current date and time
    now = datetime.utcnow()
    # Set the time range for the next 24 hours
    start_time = Time(now)
    end_time = Time(now + timedelta(days=1))
    time_range = Time([start_time, end_time])
    plots.plot_airmass(target, observer, time_range)
    plt.show()


# Data Handling and Visualization
def read_fits(file_path):
    # Placeholder: Implement FITS file reading
    pass


# Get spectrum data from SDSS
def get_spectrum(object_name):
    """
    Gets the spectrum of an object from SDSS.

    :param object_name: Name of the celestial object.
    :return: Spectrum of the object.
    """
    # Convert object_name to SkyCoord object
    co = SkyCoord.from_name(object_name)

    # Query SDSS
    xid = SDSS.query_region(co, spectro=True, radius=10 * u.arcsec)

    # Check if any matches are found
    if xid is None:
        return "No spectrum found for the object."

    # Extract plate, mjd, fiberID from the first match
    plate = xid['plate'][0]
    mjd = xid['mjd'][0]
    fiberID = xid['fiberID'][0]

    # Get spectra
    sp = SDSS.get_spectra(plate=plate, mjd=mjd, fiberID=fiberID)

    return sp[0][1].data['flux']


def load_catalog(file_path):
    return pd.read_csv(file_path)


# Additional Utility Functions
def angular_distance(ra1, dec1, ra2, dec2):
    """ Calculates the angular distance between two points on the celestial sphere
    :param ra1: ra in degrees of obj 1
    :param dec1: dec1 in degrees of obj
    :param ra2: ra in degrees of obj 2
    :param dec2: dec2 in degrees of obj 2
    :return: angular distance in degrees
    """

    c1 = SkyCoord(ra=ra1 * u.degree, dec=dec1 * u.degree, frame='icrs')
    c2 = SkyCoord(ra=ra2 * u.degree, dec=dec2 * u.degree, frame='icrs')
    return c1.separation(c2).degree


def velocity_redshift(online_file_url):
    # Placeholder: Fetch and calculate redshift from velocity data in an online file
    pass


def redshift_velocity(online_file_url):
    # Placeholder: Fetch and calculate velocity from redshift data in an online file
    pass

# Add more functions as needed
