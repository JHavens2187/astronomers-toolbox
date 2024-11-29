import unittest
import rasterio
from astronomers_toolbox import AstroToolbox  # Replace with the actual module name


class TestGetLightPollutionValue(unittest.TestCase):
    def setUp(self):
        self.toolbox = AstroToolbox()  # Replace with the actual class name

    def test_get_light_pollution_value_outside_bounds(self):
        observer_lat = 91  # Outside the valid latitude range
        observer_lon = 0

        with self.assertRaises(ValueError):
            self.toolbox.get_light_pollution_value(observer_lat=observer_lat, observer_lon=observer_lon)


if __name__ == '__main__':
    unittest.main()
