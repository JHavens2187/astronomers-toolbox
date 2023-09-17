import astropy.units as u
from astropy.coordinates import EarthLocation
import astro_toolbox as tool

# Define your location in decimal degrees
latitude = 38.9717  # Positive for North, negative for South
longitude = -95.2353  # Positive for East, negative for West

# --- Time and Motion ---

# Get the coordinates of a celestial object (e.g., Pleiades)
ra_cen, dec_cen = tool.get_object_coordinates("Pleiades")

# Print the coordinates in a readable format
print(f'Coordinates of Pleiades: RA = {ra_cen[0]}h:{ra_cen[1]}m:{ra_cen[2]:0.4}s, Dec = {dec_cen[0]}˚:{dec_cen[1]}\'{dec_cen[2]:0.4}\"')

# Calculate the precession of the object's coordinates
prec_ra, prec_dec = tool.precession(
    tool.hms2deg(*ra_cen),
    tool.dms2deg(*dec_cen)
)

# Convert precession angles from degrees to HMS and DMS
prec_ra_hms = tool.deg2hms(prec_ra)
prec_dec_dms = tool.deg2dms(prec_dec)

# Print the precession values
print(f'Precession: RA = {prec_ra_hms[0]}h:{prec_ra_hms[1]}m:{prec_ra_hms[2]:0.4}s, Dec = {prec_dec_dms[0]}˚:{prec_dec_dms[1]}\'{prec_dec_dms[2]:0.4}\"')

# --- Observational Planning ---

# Define the observer's location
observer_location = EarthLocation(lat=33.356389 * u.deg, lon=tool.dms2deg(-116, 51, 54) * u.deg, height=1712 * u.m)

# Calculate the observable time for the object
print("Observable Time for Pleiades:", tool.observable_time(observer_location, "Pleiades"))

# Get the positions of planets for a specific date and location
positions = tool.planet_positions([2023, 9, 8], [latitude, longitude])

# Extract Mars' altitude and azimuth
mars_altitude, mars_azimuth = positions[2][:2]

# Print Mars' altitude and azimuth
print(f'Mars: Altitude = {mars_altitude:0.4}˚, Azimuth = {mars_azimuth:0.4}˚')

# Convert Mars' Alt-Az to RA-Dec
mars_ra_dec = tool.altaz2radec(mars_altitude, mars_azimuth, latitude, longitude, lst=tool.LST(latitude, longitude)[0])

# Print Mars' RA and Dec
print(f'Mars in RA-Dec: RA = {mars_ra_dec[0]:0.4}˚, Dec = {mars_ra_dec[1]:0.4}˚')
