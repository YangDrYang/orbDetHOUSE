import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import orbit_propagator_wrapper
import unittest


class TestOrbitPropagatorWrapper(unittest.TestCase):

    def setUp(self):
        # Create an instance of the wrapper
        self.wrapper = orbit_propagator_wrapper.OrbitPropagatorWapper(
            "../yamls/config_orb.yml"
        )

    def test_propagate(self):
        # Call the propagate method and check the result
        result = self.wrapper.propagate()
        # Add assertions to verify the result
        self.assertIsNotNone(result)
        # Add more specific assertions based on expected result structure
        # Example: self.assertEqual(result['some_key'], expected_value)


if __name__ == "__main__":
    unittest.main()
