import datetime
import math


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
    # Define the current date and time
    now = datetime.datetime.utcnow()

    # Calculate the Julian date
    J0 = 367 * now.year - math.floor(7 * (now.year + math.floor((now.month + 9) / 12)) / 4) + math.floor(
        275 * now.month / 9) + now.day + 1721013.5
    J1 = (now.hour + now.minute / 60 + now.second / 3600) / 24
    JD = J0 + J1

    # Calculate the number of centuries since J2000.0
    T = (JD - 2451545) / 36525

    # Calculate the Greenwich mean sidereal time
    GMST = 280.46061837 + 360.98564736629 * (JD - 2451545) + 0.000387933 * T ** 2 - (T ** 3) / 38710000

    # Convert the Greenwich mean sidereal time to degrees
    GMST_deg = GMST % 360

    # Calculate the local sidereal time
    LST = GMST_deg + longitude

    # Convert the local sidereal time to hours, minutes, and seconds
    LST_hours = math.floor(LST / 15)
    LST_minutes = math.floor((LST / 15 - LST_hours) * 60)
    LST_seconds = ((LST / 15 - LST_hours) * 60 - LST_minutes) * 60

    # Print the local sidereal time
    return LST_hours, LST_minutes, LST_seconds


'''
LST_hours, LST_minutes, LST_seconds = LST()
print("Local Sidereal Time: {}h {}m {:.2f}s".format(LST_hours, LST_minutes, LST_seconds))
'''
