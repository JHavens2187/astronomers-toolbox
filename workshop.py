from astro_toolbox import AstroToolbox
import math

# Initialize AstroToolbox object
toolbox = AstroToolbox(ra=10.684, dec=41.269, object_name="Andromeda", observer_location=[38.9717, -95.2353], date=[2023, 9, 21], epoch='2000-01-01')

# Utility Functions
print("Degrees to Radians:", toolbox.deg2rad(180))
print("Radians to Degrees:", toolbox.rad2deg(math.pi))
print("Degrees to HMS:", toolbox.deg2hms(180))
print("HMS to Degrees:", toolbox.hms2deg(12, 0, 0))

# Coordinate Transformations
print("Equatorial to Galactic:", toolbox.equatorial_to_galactic())
print("Galactic to Equatorial:", toolbox.galactic_to_equatorial(121.174, -21.5728))

# Local Sidereal Time
print("Local Sidereal Time:", toolbox.LST())

# Altitude and Azimuth
alt, az = toolbox.radec2altaz(lst=toolbox.LST())
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
print(f'Planetary Phase of Mars:, {toolbox.planetary_phase("mars barycenter"):0.4}')
