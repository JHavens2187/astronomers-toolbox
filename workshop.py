
import math
import datetime
from astro_toolbox import AstroToolbox
from astropy.time import Time

# Welcome to the AstroToolbox Workshop!
# This script is designed to help you learn how to use the AstroToolbox for various astronomical calculations.

# Step 1: Initializing the Toolbox
# First, we create an instance of the AstroToolbox class.
# Here, we're using the Andromeda Galaxy as an example, but you can replace these values with your target object's data.
toolbox = AstroToolbox(ra=10.684, dec=41.269, object_name="Andromeda",
                       observer_location=[38.97106, -95.25498],
                       date=[2023, 10, 24], epoch='2000-01-01')
print("AstroToolbox initialized for:", toolbox.object_name)

# Best Practice: Always double-check your inputs, especially units (degrees, radians, etc.), to avoid common errors.

# Step 2: Conversion Functions
# AstroToolbox provides various conversion functions. Here are some examples:
print("Degrees to Radians:", toolbox.deg2rad(180))  # Converting 180 degrees to radians
print("Radians to Degrees:", toolbox.rad2deg(math.pi))  # Converting Pi radians to degrees
# Tip: Use these conversion functions to ensure your data is in the correct format for other functions.

# Step 3: Coordinate Transformations
# These functions are useful for converting between different celestial coordinate systems.
print("Equatorial to Galactic:", toolbox.equatorial_to_galactic())
print("Galactic to Equatorial:", toolbox.galactic_to_equatorial(0, 0))
# Best Practice: When transforming coordinates, verify the input coordinate system to ensure accurate conversions.

# Step 4: Flux Conversions
# Demonstrating how to convert flux values between different units.
print("CGS Angstroms to Janskys:", toolbox.cgs_angs2jy(148.6e-11, 5450))
print("Janskys to CGS Angstroms:", toolbox.jy2cgs_angs(10, 5450))
# Tip: Flux conversions are crucial in astrophysics for comparing the brightness of objects in different units.

# Czernik 44 Example
source_RA = toolbox.hms2deg(22, 52, 6)
source_DEC = toolbox.dms2deg(58, 17, 00)
fov = 27.9139  # Field of View in arcminutes
fov_degrees = fov / 60
# Uncomment the line below to display the finder scope for Czernik 44
# toolbox.get_finder_scope(fov_degrees, ra=source_RA, dec=source_DEC, object_name="Czernik 44")

targets = [
    {'ra': toolbox.hms2deg(22, 52, 6), 'dec': toolbox.dms2deg(58, 17, 00), 'name': 'Czernik 44'},
    'Polaris'
]

# Uncomment the line below to display the sky plot with airmass and altitude
toolbox.plot_sky(targets=targets, observation_time=Time('2024-11-27 23:00:00'))

# Coordinate Transformations
print("Equatorial to Galactic:", toolbox.equatorial_to_galactic())
print("Galactic to Equatorial:", toolbox.galactic_to_equatorial(121.174, -21.5728))

# Local Sidereal Time
print("Local Sidereal Time:", toolbox.LST())

# Altitude and Azimuth
alt, az = toolbox.radec2altaz()
print("Altitude:", alt, "Azimuth:", az)

# Proper Motion
new_ra, new_dec = toolbox.proper_motion(pm_ra=0.0001, pm_dec=0.0002)
print("New RA:", new_ra, "New Dec:", new_dec)

# Precession
prec_ra, prec_dec = toolbox.precession()
print("Precessed RA:", prec_ra, "Precessed Dec:", prec_dec)

# Planet Positions
print("Planet Positions:", toolbox.planet_positions())

# Planetary Phase
planet = 301  # representing the moon
print(f'Planetary Phase of {planet}:, {toolbox.planetary_phase(planet, date=[2023, 10, 28]):0.4}')

limiting_magnitude = toolbox.get_limiting_mag(38.97113, -95.25416)
print("Limiting magnitude Value:", limiting_magnitude)

# What object is that?
perseid_coord = toolbox.get_object_coordinates('Polaris')
print(perseid_coord)
results_radec = toolbox.find_objects_by_coordinates(perseid_coord, radius=60)

# Print the results
print("RA/Dec results:", results_radec)

# End of the Workshop
# Congratulations! You have now explored the key features of the AstroToolbox.
# Remember, practice is key to mastering these tools and concepts. Happy stargazing!
