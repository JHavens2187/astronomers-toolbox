# AstroToolbox Workshop
# J. Havens, 2023

from astro_toolbox import AstroToolbox
import datetime

# Initialize AstroToolbox object
toolbox = AstroToolbox(ra=10.684, dec=41.269, object_name="Pleiades", observer_location=(38.9717, -95.2353), date='2023-09-08')

# Utility Functions
print("Utility Functions:")
print("Degrees to Radians:", toolbox.deg2rad(180))
print("Radians to Degrees:", toolbox.rad2deg(3.14159))
print("Degrees to HMS:", toolbox.deg2hms(180))
print("HMS to Degrees:", toolbox.hms2deg(12, 0, 0))

# Coordinate Transformations
print("\nCoordinate Transformations:")
print("Equatorial to Galactic:", toolbox.equatorial_to_galactic())
print("Galactic to Equatorial:", toolbox.galactic_to_equatorial(121.174, -21.573))

# Time and Motion
print("\nTime and Motion:")
print("Proper Motion:", toolbox.proper_motion(10.684, 41.269))
print("Precession:", toolbox.precession())

# Celestial Object Information from SIMBAD
print("\nCelestial Object Information:")
print("Object Coordinates:", toolbox.get_object_coordinates())

# Observational Planning
print("\nObservational Planning:")
print("Observable Time:", toolbox.observable_time())

# Additional Functions
print("\nAdditional Functions:")
print("LST:", toolbox.LST())
print("Planet Positions:", toolbox.planet_positions())
print("Planetary Phase (Mars):", toolbox.planetary_phase('mars barycenter'))

