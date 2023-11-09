import math

from astro_toolbox import AstroToolbox

# Initialize AstroToolbox object
toolbox = AstroToolbox(ra=10.684, dec=41.269, object_name="Andromeda", observer_location=[38.97106, -95.25498],
                       date=[2023, 10, 24], epoch='2000-01-01')

# Utility Functions
print("Degrees to Radians:", toolbox.deg2rad(180))
print("Radians to Degrees:", toolbox.rad2deg(math.pi))
print("Degrees to HMS:", toolbox.deg2hms(180))
print(toolbox.cgs_angs2jy(148.6e-11, 5450), toolbox.cgs_angs2jy(4.2e-11, 5450))
print(toolbox.jy2cgs_angs(1472.2, 5450))

# Czernik 44
source_RA = toolbox.hms2deg(22, 52, 6)
source_DEC = toolbox.dms2deg(58, 17, 00)
fov = 27.9139  # in arcmin
fov_degrees = fov/60
# uncomment this to show the finder scope for the star above
# toolbox.get_finder_scope(fov_degrees, ra=source_RA, dec=source_DEC, object_name="Czernik 44")

targets = [
    {'ra': toolbox.hms2deg(22, 52, 6), 'dec': toolbox.dms2deg(58, 17, 00), 'name': 'Czernik 44'},
    'Polaris'
]

# uncomment this to show the sky plot with airmass and altitude
# toolbox.plot_sky(targets=targets, observation_time=Time('2023-10-10 18:00:00'))

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
planet = "mars barycenter"
planet = 301
print(f'Planetary Phase of {planet}:, {toolbox.planetary_phase(planet, date=[2023, 10, 28]):0.4}')

limiting_magnitude = toolbox.get_limiting_mag(38.97113, -95.25416)
print("Limiting magnitude Value:", limiting_magnitude)

# what object is that?
perseid_coord = toolbox.get_object_coordinates('Perseid')
print(perseid_coord)
results_radec = toolbox.find_objects_by_coordinates(perseid_coord, radius=60)

# Print the results
print("RA/Dec results:", results_radec)
