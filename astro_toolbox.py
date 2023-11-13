# a toolbox for astronomers
# J. Havens, 2023

from datetime import datetime, timedelta
import datetime
import math

from astroplan import AltitudeConstraint, AirmassConstraint, MoonSeparationConstraint
from astroplan import FixedTarget
from astroplan import Observer
from astroplan.plots import plot_airmass, plot_altitude
from astropy.coordinates import SkyCoord, solar_system_ephemeris, FK5, EarthLocation, AltAz
from astropy.time import Time
from astroquery.sdss import SDSS
from astroquery.simbad import Simbad
from astroquery.skyview import SkyView
from skyfield.api import load, Topos
import LST_calculator as lst
import astroplan
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import requests






# TODO: Implement error handling for functions to make the toolbox more robust.


class AstroToolbox:

    def __init__(self, ra, dec, object_name, observer_location=None, date=None, epoch='2000-01-01'):
        """
        Initialize the AstroToolbox with common parameters.

        :param (float) ra: Right Ascension in degrees.
        :param (float) dec: Declination in degrees.
        :param (str) object_name: Name of the object
        :param (list) observer_location:  observer's location in [lat, lon] format. Default is Lawrence.
        :param (str) date: Date for observations in 'YYYY-MM-DD' format. Default is now.
        :param (str) epoch: Epoch of reference, default to J2000.
        """
        self.ra = ra  # Right Ascension in degrees
        self.dec = dec  # Declination in degrees
        self.object_name = object_name
        self.observer_location = observer_location if observer_location else [38.9717,
                                                                              -95.2353]  # Default to Lawrence, KS
        self.date = date if date else datetime.utcnow().strftime('%Y-%m-%d')  # Default to current date
        self.epoch = epoch  # Default to J2000

    # Utility Functions
    def deg2rad(self, deg):
        return deg * np.pi / 180

    def rad2deg(self, rad):
        return rad * 180 / np.pi

    def deg2hms(self, deg):
        hr = deg / 15
        min = (hr - math.floor(hr)) * 60
        sec = (min - math.floor(min)) * 60
        return math.floor(hr), math.floor(min), sec

    def hms2deg(self, hr, min, sec):
        return (hr + (min / 60) + (sec / 3600)) * 15

    def dms2deg(self, deg, min, sec):
        return deg + (min / 60) + (sec / 3600)

    def deg2dms(self, decdeg):
        deg = math.floor(decdeg)
        decmin = (decdeg - deg) * 60
        decsec = (decmin - math.floor(decmin)) * 60
        return deg, math.floor(decmin), decsec

    def dms2hms(self, deg, min, sec):
        decdeg = self.dms2deg(deg, min, sec)
        dechr = 15 * decdeg
        decmin = (dechr - math.floor(dechr)) * 60
        decsec = (decmin - math.floor(decmin)) * 60
        return math.floor(dechr), math.floor(decmin), decsec

    def hms2dms(self, hr, min, sec):
        decdeg = self.hms2deg(hr, min, sec)
        deg = math.floor(decdeg)
        decmin = (decdeg - deg) * 60
        decsec = (decmin - math.floor(decmin)) * 60
        return deg, math.floor(decmin), decsec

    def cgs_angs2jy(self, cgs_angs, cen_wave):
        """
        converts flux in units of erg/s/cm^2/angstrom to janskys
        :param cgs_angs: flux in cgs units
        :param cen_wave: reference wavelength in angstroms
        :return: flux in Jy units
        """
        c = 2.998e18  # angstroms/s
        fhz = (cen_wave ** 2 * cgs_angs) / c
        return fhz * 1e23

    def jy2cgs_angs(self, jy, cen_wave):
        """
        converts flux in units of janskys to erg/s/cm^2/angstrom
        :param jy: flux in Jy units
        :param cen_wave: reference wavelength in angstroms
        :return: flux in units of erg/s/cm^2/angstrom
        """
        c = 2.998e18  # angstroms/s
        fhz = jy / 1e23
        return fhz * (c / (cen_wave ** 2))

    def LST(self, observer_location=None):
        """
        calculates the local sidereal time
        :param (list) observer_location: observer's location in lat lon format
        :return: LST_hours, LST_minutes, LST_seconds
        """
        # define your location
        if observer_location is None:
            observer_location = self.observer_location
        lat, lon = observer_location[0], observer_location[1]
        if lat and lon is None:  # assume Lawrence, KS
            latitude = 38.9717  # Your latitude in decimal degrees (positive for North, negative for South)
            longitude = -95.2353  # Your longitude in decimal degrees (positive for East, negative for West)
        else:
            latitude = lat
            longitude = lon
        return lst.LST(latitude, longitude)

    # Coordinate Transformations
    def equatorial_to_galactic(self, ra=None, dec=None):
        """
        converts equatorial coordinates to galactic coordinates
        :param ra: ra in degrees
        :param dec: dec in degrees
        :return: l and b IN DEGREES
        """
        if ra is None:
            ra = self.ra
        if dec is None:
            dec = self.dec
        c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
        galactic = c.galactic
        return galactic.l.degree, galactic.b.degree

    def galactic_to_equatorial(self, l, b):
        """
        converts galactic coordinates to equatorial coordinates
        :param l: l in degrees
        :param b: b in degrees
        :return: ra and dec IN DEGREES
        """
        c = SkyCoord(l=l * u.degree, b=b * u.degree, frame='galactic')
        equatorial = c.icrs
        return equatorial.ra.degree, equatorial.dec.degree

    # Function to convert between equatorial coordinates and altitude and azimuth
    def radec2altaz(self, ra=None, dec=None, observer_location=None):
        """
        Converts equatorial coordinates to altitude and azimuth using both
        manual calculation and SkyCoord transformation.
        :param ra: right ascension of object in degrees
        :param dec: declination of object in degrees
        :param observer_location: observer's location in lat lon format
        :return: altitude and azimuth of object in degrees
        """
        if ra is None:
            ra = self.ra
        if dec is None:
            dec = self.dec
        if observer_location is None:
            observer_location = self.observer_location

        # Manual calculation
        ra_rad = math.radians(ra)
        dec_rad = math.radians(dec)
        lat_rad = math.radians(observer_location[0])
        lst = self.LST()
        lst_rad = math.radians(lst[0] * 15)

        sin_alt = math.sin(dec_rad) * math.sin(lat_rad) + math.cos(dec_rad) * math.cos(lat_rad) * math.cos(
            lst_rad - ra_rad)
        alt = math.degrees(math.asin(sin_alt))

        sin_az = math.sin(lst_rad - ra_rad)
        cos_az = (math.sin(dec_rad) - math.sin(alt) * math.sin(lat_rad)) / (math.cos(alt) * math.cos(lat_rad))
        az = math.degrees(math.atan2(sin_az, cos_az))
        az = az % 360  # Ensure azimuth is within 0-360 degrees

        # Ensure the RA and Dec are in Astropy units
        ra = ra * u.deg
        dec = dec * u.deg

        # Get the observer's Earth location
        location = EarthLocation(lat=observer_location[0] * u.deg, lon=observer_location[1] * u.deg)

        # Use the current time for the observation
        current_time = Time.now()

        # Create the SkyCoord object for the given RA and Dec
        sky_coord = SkyCoord(ra=ra, dec=dec, frame='icrs')

        # Create the AltAz frame with the current time and observer location
        altaz_frame = AltAz(obstime=current_time, location=location)

        # Transform the SkyCoord object to the AltAz frame
        altaz_coord = sky_coord.transform_to(altaz_frame)

        # Extract the altitude and azimuth in degrees
        alt_sky = altaz_coord.alt.deg
        az_sky = altaz_coord.az.deg

        # Check if values are close enough; if not, issue a warning
        if not (math.isclose(alt, alt_sky, abs_tol=0.1) and math.isclose(az, az_sky, abs_tol=0.1)):
            Warning("Calculated AltAz values differ from SkyCoord transformation by more than 0.1 degrees.")

        # Return the SkyCoord values regardless of the check
        return alt_sky, az_sky

    def altaz2radec(self, alt, az, lst=None, observer_location=None):
        """
        Converts horizontal coordinates (Altitude, Azimuth) to equatorial coordinates (RA, DEC).

            :param (float) alt: Altitude in degrees
            :param (float) az: Azimuth in degrees
            :param (list) observer_location:  observer's location in lat lon format
            :param (float) lst: Local Sidereal Time in hours

            :return tuple: Right Ascension and Declination in degrees
        """
        if lst is None:
            lst = self.LST()
        if observer_location is None:
            observer_location = self.observer_location
        # Convert angles to radians
        alt_rad = math.radians(alt)
        az_rad = math.radians(az)
        lat_rad = math.radians(observer_location[0])
        lst_rad = math.radians(lst[0] * 15)  # Convert LST to degrees and then to radians

        # Calculate Declination
        sin_dec = math.sin(alt_rad) * math.sin(lat_rad) + math.cos(alt_rad) * math.cos(lat_rad) * math.cos(az_rad)
        dec = math.degrees(math.asin(sin_dec))

        # Calculate Right Ascension
        sin_ra = -math.sin(az_rad) * math.cos(alt_rad)
        cos_ra = (math.sin(alt_rad) - sin_dec * math.sin(lat_rad)) / (math.cos(math.asin(sin_dec)) * math.cos(lat_rad))
        ra = lst[0] * 15 - math.degrees(math.atan2(sin_ra, cos_ra))  # Convert lst from hours to degrees

        # Normalize RA to be between 0 and 360
        ra = ra % 360

        # Ensure the altitude and azimuth are in Astropy units
        alt = alt * u.deg
        az = az * u.deg

        # Get the observer's Earth location
        location = EarthLocation(lat=observer_location[0]*u.deg, lon=observer_location[1]*u.deg)

        # Use the current time for the observation
        obstime = Time.now()

        # Create the AltAz frame and transform it to ICRS (RA-Dec)
        altaz_frame = AltAz(obstime=obstime, location=location)
        altaz_coord = SkyCoord(alt=alt, az=az, frame=altaz_frame)
        radec = altaz_coord.transform_to('icrs')

        ra_sky = radec.ra.deg
        dec_sky = radec.dec.deg

        # Check if values are close enough; if not, issue a warning
        if not (math.isclose(ra, ra_sky, abs_tol=0.1) and math.isclose(dec, dec_sky, abs_tol=0.1)):
            Warning("Calculated RA-Dec values differ from SkyCoord transformation by more than 0.1 degrees.")

        # Return the SkyCoord values regardless of the check
        return ra_sky, dec_sky

        return ra / 15, dec

    # Time and Motion
    def proper_motion(self, pm_ra, pm_dec, ra=None, dec=None, epoch=None):
        """
        calculates the proper motion of an object
        :param ra: ra in degrees
        :param dec: dec in degrees
        :param pm_ra: proper motion in ra in deg/yr
        :param pm_dec: proper motion in dec in deg/yr
        :param epoch: epoch of reference, default to J2000
        :return: current position as a result of proper motion IN DEGREES
        """
        if ra is None:
            ra = self.ra
        if dec is None:
            dec = self.dec
        if epoch is None:
            epoch = self.epoch
        current_time = Time(datetime.datetime.now())
        delta_time = current_time - Time(epoch)
        new_ra = ra + pm_ra * delta_time.to(u.yr).value
        new_dec = dec + pm_dec * delta_time.to(u.yr).value
        return new_ra, new_dec

    def precession(self, ra=None, dec=None, epoch=None):
        """
        calculates the precession of an object
        :param ra: ra in degrees
        :param dec: dec in degrees
        :param epoch: epoch of reference, default to J2000
        :return: updated ra and dec IN DEGREES
        """
        if ra is None:
            ra = self.ra
        if dec is None:
            dec = self.dec
        if epoch is None:
            epoch = self.epoch
        epoch_time = Time(epoch)
        current_time = Time(datetime.datetime.now())
        c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame=FK5, equinox=epoch_time)
        c_prec = c.transform_to(FK5(equinox=current_time))
        return c_prec.ra.degree, c_prec.dec.degree

    def planet_positions(self, date=None, observer_location=None):
        """
        calculates the positions of the planets in the solar system
        :param (list) observer_location:  observer's location in lat lon format
        :param (list) date: date of observation in YYYY-MM-DD format
        :return: (list) list of lists containing positions of barycenter of planets in the solar system
        """
        if date is None:
            date = self.date
        if observer_location is None:
            observer_location = self.observer_location
        ts = load.timescale()
        t = ts.utc(int(date[0]), int(date[1]), int(date[2]))
        eph = load('de421.bsp')
        earth = eph['earth']
        planets = ['mercury barycenter', 'venus barycenter', 'mars barycenter', 'jupiter barycenter',
                   'saturn barycenter']
        positions = {}
        positions_list = []

        for planet in planets:
            astrometric = (earth + Topos(latitude_degrees=observer_location[0],
                                         longitude_degrees=observer_location[1])).at(
                t).observe(eph[planet])
            alt, az, _ = astrometric.apparent().altaz()
            positions[planet] = {'Altitude': alt.degrees, 'Azimuth': az.degrees}
            positions_list.append([alt.degrees, az.degrees])

        return positions_list

    def planetary_phase(self, planet, date=None, observer_location=None):
        """
        Calculate the phase angle of a given planet as observed from Earth. \n
        0˚ - full phase \n
        90˚ - first quarter phase \n
        180˚ - new phase \n
        availible objects are:
        0 SOLAR SYSTEM BARYCENTER, 1 MERCURY BARYCENTER, 2 VENUS BARYCENTER, 3 EARTH BARYCENTER, 4 MARS BARYCENTER,
        5 JUPITER BARYCENTER, 6 SATURN BARYCENTER, 7 URANUS BARYCENTER, 8 NEPTUNE BARYCENTER, 9 PLUTO BARYCENTER,
        10 SUN, 199 MERCURY, 299 VENUS, 301 MOON, 399 EARTH
        :param (str) planet: The name of the planet for which to calculate the phase angle. Must be a valid name recognized by the Skyfield library
        (e.g., 'mars barycenter').
        :param (list) date: (optional) date of observation. If not provided, self.date will be used.
        :param (list) observer_location: (optional) location of observer. If not provided, self.observer_location will be used.
        :return: (float) The phase angle of the planet in degrees. The phase angle is the angle between the Sun and the planet as observed from Earth.
        It gives an idea of how much of the planet's disk is illuminated from the observer's perspective.\n
        """
        if date is None:
            date = self.date
        if observer_location is None:
            observer_location = self.observer_location
        ts = load.timescale()
        t = ts.utc(date[0], month=date[1], day=date[2])
        eph = load('de440.bsp')
        sun = eph['sun']
        earth = eph['earth']

        # Use actual observer's latitude and longitude
        observer_latitude = observer_location[0]
        observer_longitude = observer_location[1]

        astrometric_planet = (
                earth + Topos(latitude_degrees=observer_latitude, longitude_degrees=observer_longitude)).at(
            t).observe(eph[planet]).apparent()
        astrometric_sun = (earth + Topos(latitude_degrees=observer_latitude, longitude_degrees=observer_longitude)).at(
            t).observe(sun).apparent()

        # Calculate the phase angle using separation_from method
        phase_angle = astrometric_planet.separation_from(astrometric_sun)

        return phase_angle.degrees

    # Celestial Object Information from SIMBAD
    def get_object_coordinates(self, object_name=None):
        """
        gets the coordinates of an object from SIMBAD
        :param object_name: the name of the object
        :return: ra and dec of object IN HMS/DMS
        """
        if object_name is None:
            object_name = self.object_name
        result_table = Simbad.query_object(object_name)
        coord = SkyCoord(ra=result_table['RA'].data[0], dec=result_table['DEC'].data[0], unit=(u.deg, u.deg))
        ra = self.deg2hms(coord.ra.degree)
        dec = self.deg2dms(coord.dec.degree)
        return ra, dec

    def get_light_pollution_value(self, observer_lat=None, observer_lon=None):
        geotiff_path = '/Users/joehavens/Documents/Documents/pythonProject/Astro/astronomers-toolbox/wa_2015_original/World_Atlas_2015.tif'
        if (observer_lat and observer_lon) is None:
            observer_lat, observer_lon = self.observer_location
        # Open the GEOTiff file
        with rasterio.open(geotiff_path) as dataset:
            # Optionally: Transform the observer's coordinates if necessary
            # ...
            # Get the pixel value at the given coordinates (which corresponds to light pollution)
            value = dataset.read(1)  # Assuming light pollution data is in band 1
            row, col = dataset.index(observer_lon, observer_lat)
            return value[row, col]

    def get_limiting_mag(self, observer_lat=None, observer_lon=None):
        """
        Calculates the limiting magnitude of the sky.
        :param observer_lon: observer's longitude in degrees
        :param observer_lat: observer's latitude in degrees
        :return: limiting magnitude
        """

        def interpolate_limiting_mag(sqm_value, mag_range_low, mag_range_high, sqm_low, sqm_high):
            """
            Interpolates the limiting magnitude based on SQM value.
            :param sqm_value: SQM value to interpolate
            :param mag_range_low: Lower limit of the magnitude range
            :param mag_range_high: Upper limit of the magnitude range
            :param sqm_low: Lower limit of the SQM range
            :param sqm_high: Upper limit of the SQM range
            :return: Interpolated limiting magnitude
            """
            # Proportion of SQM value within its range
            proportion = (sqm_value - sqm_low) / (sqm_high - sqm_low)
            # Interpolated limiting magnitude
            return mag_range_low + proportion * (mag_range_high - mag_range_low)

        if (observer_lat and observer_lon) is None:
            observer_lat, observer_lon = self.observer_location

        # Define the natural sky brightness
        natural_brightness = 0.171168465
        total_brightness = self.get_light_pollution_value(observer_lat=observer_lat,
                                                          observer_lon=observer_lon) + natural_brightness
        SQM = np.log10(total_brightness / 108000000) / -0.4

        # Determine the limiting magnitude based on the SQM value
        if SQM >= 22.00:
            lim_mag = 7.5  # or higher
        elif 21.99 >= SQM > 21.89:
            lim_mag = interpolate_limiting_mag(SQM, 7.0, 7.49, 21.89, 21.99)
        elif 21.89 >= SQM > 21.69:
            lim_mag = interpolate_limiting_mag(SQM, 6.5, 6.99, 21.69, 21.89)
        elif 21.69 >= SQM > 20.49:
            lim_mag = interpolate_limiting_mag(SQM, 6.0, 6.49, 20.49, 21.69)
        elif 20.49 >= SQM > 19.50:
            lim_mag = interpolate_limiting_mag(SQM, 5.5, 5.99, 19.50, 20.49)
        elif 19.50 >= SQM > 18.94:
            lim_mag = interpolate_limiting_mag(SQM, 5.0, 5.49, 18.94, 19.50)
        elif 18.94 >= SQM > 18.38:
            lim_mag = interpolate_limiting_mag(SQM, 4.5, 4.99, 18.38, 18.94)
        elif 18.38 >= SQM:
            lim_mag = 4.49  # or lower
        else:
            lim_mag = 4.0  # or lower, this is an edge case if SQM is below the table's range

        return lim_mag

    def find_objects_by_coordinates(self, ra=None, dec=None, alt=None, az=None, radius=1,
                                    lim_mag=None, observer_lat=None, observer_lon=None):
        """
        Find celestial objects near the given coordinates within a specified limiting magnitude.

        :param coord_type: Type of coordinates provided ('radec' or 'altaz')
        :param ra: Right Ascension in degrees (if coord_type is 'radec')
        :param dec: Declination in degrees (if coord_type is 'radec')
        :param alt: Altitude in degrees (if coord_type is 'altaz')
        :param az: Azimuth in degrees (if coord_type is 'altaz')
        :param radius: Search radius in arcminutes
        :param lim_mag: Limiting magnitude for the search. use an absurdly large number for no limit
        :param observer_lon: observers longitude in degrees
        :param observer_lat: observers latitude in degrees
        :return: List of tuples with object information (ID, RA, Dec, Object Type)
        """
        # Set up SIMBAD with the fields you want to return
        Simbad.add_votable_fields('otype', 'flux(V)')
        current_time = Time.now()

        if (observer_lat and observer_lon) is None:
            observer_lat, observer_lon = self.observer_location
        if (ra and dec) is None:
            ra = self.ra
            dec = self.dec
        if alt and az:
            ra, dec = self.altaz2radec(alt, az)

        # Use the RA and DEC to make a SkyCoord object, which we use to query and image
        sky_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')

        # Query SIMBAD database
        result_table = Simbad.query_region(sky_coord)

        if result_table is None:
            return []

        if lim_mag is None:
            lim_mag = self.get_limiting_mag()

        # Filter results by limiting magnitude if specified
        if lim_mag is not None:
            result_table = result_table[result_table['FLUX_V'] <= lim_mag]

        # Now filter by radius
        result_table = result_table[result_table['RA'].separation(ra * u.deg) <= radius * u.arcmin]

        # Extract the list of objects
        objects_list = []
        for obj in result_table:
            object_id = obj['MAIN_ID']
            object_ra = obj['RA']
            object_dec = obj['DEC']
            object_type = obj['OTYPE']
            object_mag = obj['FLUX_V']  # This is the visual magnitude from SIMBAD
            objects_list.append((object_id, object_ra, object_dec, object_type, object_mag))

        # If limiting magnitude is set, filter the objects
        if lim_mag is not None:
            objects_list = [obj for obj in objects_list if obj[4] is not None and obj[4] <= lim_mag]

            # Now check for planets and add them to the list if they are within the search radius
        with solar_system_ephemeris.set('builtin'):
            observer_location = EarthLocation(lat=observer_lat * u.deg, lon=observer_lon * u.deg)
            for planet_name in ['mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']:
                planet = get_body(planet_name, Time.now(), observer_location)
                planet_coord = SkyCoord(ra=planet.ra, dec=planet.dec, frame='icrs')

                # Check if the planet is within the search radius
                if planet_coord.separation(sky_coord) <= radius * u.arcmin:
                    # Add the planet to the objects list
                    objects_list.append(
                        (planet_name.capitalize(), planet.ra.to_string(), planet.dec.to_string(), 'Planet', None))

        return objects_list

    # Observational Planning
    def observable_time(self, observer_location=None, object_name=None):
        """
        calculates if an object is observable in the next 24 hours
        :param (list) observer_location:  observer's location in lat lon format
        :param object_name: name of the object
        :return: boolean if object is observable
        """
        if observer_location is None:
            observer_location = self.observer_location
        if object_name is None:
            object_name = self.object_name
        observer_location = EarthLocation(lat=observer_location[0] * u.deg, lon=observer_location[1] * u.deg)
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

    def plot_sky(self, observer_location=None, targets=None, observation_time=None):
        """
        Plots the airmass and altitude of multiple objects.
        :param observer_location: observer's location in lat lon format
        :param targets: list of target names or coordinates
        :param observation_time: Time of observation in 'YYYY-MM-DD HH:MM:SS' format
        :return: plot of the airmass and altitude of the objects
        """
        if observer_location is None:
            observer_location = self.observer_location  # Replace with your default observer location

        observer_location = EarthLocation(lat=observer_location[0] * u.deg, lon=observer_location[1] * u.deg,
                                          height=1 * u.m)
        # Create an Observer object
        observer = Observer(location=observer_location)

        if observation_time is None:
            now = datetime.utcnow()
            observation_time = Time(now)

        fig, ax1 = plt.subplots()

        ax2 = ax1.twinx()

        for target in targets:
            if isinstance(target, str):
                target = FixedTarget.from_name(target)
            elif isinstance(target, dict):
                coord = SkyCoord(ra=target['ra'] * u.deg, dec=target['dec'] * u.deg)
                target = FixedTarget(coord=coord, name=target.get('name', 'Unknown'))

            plot_airmass(target, observer, observation_time, ax=ax1)
            plot_altitude(target, observer, observation_time, ax=ax2)

        # Plot a line and shaded region at and above an elevation of 60 degrees
        xmin, xmax = ax1.get_xlim()
        x_line = np.linspace(xmin, xmax, 100)
        y = [60] * 100
        ax2.plot(x_line, y, ls='--')
        ax2.fill_between(x_line, 60, 90, alpha=0.1, color='orange')

        ax1.set_ylabel('Airmass')
        ax2.set_ylabel('Altitude [deg]')
        plt.legend()
        plt.show()

    def get_finder_scope(self, fov, ra=None, dec=None, object_name='object', img_width=40 / 60, img_height=40 / 60):
        if ra is None:
            ra = self.ra
        if dec is None:
            dec = self.dec

        # Use the RA and DEC to make a SkyCoord object, which we use to query and image
        source_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

        ####FIX IMPLEMENTED####
        # Query an image using astroquery with the given coordinates and image width/height
        # xout = SkyView.get_images(source_coord,survey=['DSS'],height=img_height*u.deg,width=img_width*u.deg)
        xout = SkyView.get_images(source_coord, survey=['DSS'], pixels=1200)

        # Make the figure object and handle the fits image appropriately.
        fig, ax = plt.subplots(figsize=(8, 8))
        b = xout[0][0]
        ax.imshow(xout[0][0].data, aspect='equal', cmap='gray_r',
                  extent=[b.header['CRVAL1'] - (b.header['NAXIS1'] - b.header['CRPIX1']) * b.header['CDELT1'],
                          b.header['CRVAL1'] + (b.header['NAXIS1'] - b.header['CRPIX1']) * b.header['CDELT1'],
                          b.header['CRVAL2'] + (b.header['NAXIS2'] - b.header['CRPIX2']) * b.header['CDELT2'],
                          b.header['CRVAL2'] - (b.header['NAXIS2'] - b.header['CRPIX2']) * b.header['CDELT2']])
        # Overlay a rectangle indicating the fov
        rect = plt.Rectangle((ra - 0.5 * (fov), dec - 0.5 * (fov)), fov, fov, linewidth=1, fill=False,
                             color='k')

        # plot the finder scope
        ax.set_xticklabels(
            [f"{round(h, 2)}h {round(m, 2)}m {round(s, 2)}s" for h, m, s in [self.deg2hms(x) for x in ax.get_xticks()]],
            fontsize=10, rotation=30)
        ax.set_yticklabels(
            [f"{round(d, 2)}˚ {round(m, 2)}' {round(s, 2)}s" for d, m, s in [self.deg2dms(y) for y in ax.get_yticks()]],
            fontsize=10, rotation=30)
        ax.add_patch(rect)
        ax.set_xlabel(r'RA [h:m:s]')
        ax.set_ylabel(r'DEC [d:m:s]')
        plt.title(object_name)
        plt.show()

    # Get spectrum data from SDSS
    def get_spectrum(self, object_name=None):
        """
        Gets the spectrum of an object from SDSS.

        :param object_name: Name of the celestial object.
        :return: Spectrum of the object.
        """
        if object_name is None:
            object_name = self.object_name
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

    def load_catalog(self, file_path):
        return pd.read_csv(file_path)

    # Additional Utility Functions
    def angular_distance(self, ra1, dec1, ra2, dec2):
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

    # Add more functions as needed
