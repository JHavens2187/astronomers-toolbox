import math

from astropy.time import Time

from astro_toolbox import AstroToolbox

# Initialize AstroToolbox object
toolbox = AstroToolbox(ra=10.684, dec=41.269, object_name="Andromeda", observer_location=[38.9717, -95.2353],
                       date=[2023, 9, 21], epoch='2000-01-01')

# Utility Functions
print("Degrees to Radians:", toolbox.deg2rad(180))
print("Radians to Degrees:", toolbox.rad2deg(math.pi))
print("Degrees to HMS:", toolbox.deg2hms(180))
print(toolbox.cgs_angs2jy(148.6e-11, 5450), toolbox.cgs_angs2jy(4.2e-11, 5450))
print(toolbox.jy2cgs_angs(1472.2, 5450))

# Berkeley 45
source_RA1 = toolbox.hms2deg(19, 19, 12)
source_DEC1 = toolbox.dms2deg(15, 43, 00)
fov = 27.9139  # in arcmin
fov_degrees = fov/60
toolbox.get_finder_scope(fov_degrees, ra=source_RA1, dec=source_DEC1)

# # Czernik 40
# source_RA4 = toolbox.hms2deg(19, 42, 36)
# source_DEC4 = toolbox.dms2deg(21, 9, 14)
# fov = 27.9139  # in arcmin
# fov_degrees = fov/60
# toolbox.get_finder_scope(fov_degrees, ra=source_RA4, dec=source_DEC4)

# # Berkeley 51
# source_RA5 = toolbox.hms2deg(20, 11, 54)
# source_DEC5 = toolbox.dms2deg(34, 24, 6)
# fov = 27.9139  # in arcmin
# fov_degrees = fov / 60
# toolbox.get_finder_scope(fov_degrees, ra=source_RA5, dec=source_DEC5)

# # Berkeley 49
# source_RA6 = toolbox.hms2deg(19, 54, 31)
# source_DEC6 = toolbox.dms2deg(34, 38, 48)
# fov = 27.9139  # in arcmin
# fov_degrees = fov/60
# toolbox.get_finder_scope(fov_degrees, ra=source_RA6, dec=source_DEC6)

targets = [
    {'ra': 274.987, 'dec': 72.5, 'name': 'Target1'},
    {'ra': 275.987, 'dec': 12.5, 'name': 'Target2'},
    'Polaris'
]
toolbox.plot_sky(targets=targets, observation_time=Time('2023-10-10 18:00:00'))

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
