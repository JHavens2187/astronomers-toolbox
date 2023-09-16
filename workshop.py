import astropy.units as u
from astropy.coordinates import EarthLocation

import astro_toolbox as tool

# define your location
latitude = 38.9717  # Your latitude in decimal degrees (positive for North, negative for South)
longitude = -95.2353  # Your longitude in decimal degrees (positive for East, negative for West)

# Time and Motion
ra_cen, dec_cen = tool.get_object_coordinates("Pleiades")
print(f'ra: {ra_cen[0]}h:{ra_cen[1]}m:{ra_cen[2]:0.4}s, dec: {dec_cen[0]}˚:{dec_cen[1]}\':{dec_cen[2]:0.4}\"')
prec_ra_hr, prec_ra_min, prec_ra_sec = tool.deg2hms(
    tool.precession(tool.hms2deg(ra_cen[0], ra_cen[1], ra_cen[2]), tool.dms2deg(dec_cen[0], dec_cen[1], dec_cen[2]))[0])
prec_dec_hr, prec_dec_min, prec_dec_sec = tool.deg2dms(
    tool.precession(tool.hms2deg(ra_cen[0], ra_cen[1], ra_cen[2]), tool.dms2deg(dec_cen[0], dec_cen[1], dec_cen[2]))[1])
print(
    f'Precession: {prec_ra_hr}h:{prec_ra_min}m:{prec_ra_sec:0.4}s, dec: {prec_dec_hr}˚:{prec_dec_min}\':{prec_dec_sec:0.4}\"')
# Observational Planning
observer_location = EarthLocation(lat=33.356389 * u.deg, lon=tool.dms2deg(-116, 51, 54) * u.deg, height=1712 * u.m)
print("Observable Time:", tool.observable_time(observer_location, "Pleiades"))
positions = tool.planet_positions([2023, 9, 8], [38.9717, -95.2353])
mars_altitude = positions[2][0]  # Mars is the 3rd element in the list, altitude is the 1st element in the inner list
mars_azimuth = positions[2][1]  # Mars is the 3rd element in the list, azimuth is the 2nd element in the inner list
mars = positions[2]
print(f'Mars: {mars[0]:0.4}˚, {mars[1]:0.4}˚')
mars_ra_dec = tool.altaz2radec(mars[0], mars[1], 38.9717, -95.2353, lst=tool.LST(38.9717, -95.2353)[0])
print(f'Mars: {mars_ra_dec[0]:0.4}˚, {mars_ra_dec[1]:0.4}˚')

