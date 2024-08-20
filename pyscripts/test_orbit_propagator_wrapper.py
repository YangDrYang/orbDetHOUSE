import sys
import os

# Add the directory containing the .so file to the system path
lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib"))
print(f"Library path: {lib_path}")
sys.path.insert(0, lib_path)

# Print the sys.path for debugging
print("sys.path:", sys.path)

# Check if the .so file exists in the lib_path
so_file = os.path.join(lib_path, "orbit_propagator_wrapper.so")
if os.path.exists(so_file):
    print(f"Found orbit_propagator_wrapper.so at {so_file}")
else:
    print(f"orbit_propagator_wrapper.so not found at {so_file}")

try:
    import orbit_propagator_wrapper

    print("Module imported successfully")
except ModuleNotFoundError as e:
    print(f"Error: {e}")

# Assuming you have a class named OrbitPropagator in your C++ code
propagator = orbit_propagator_wrapper.OrbitPropagatorWapper("yamls/config_orb.yml")
result = propagator.propagate()
print("propagated orbit:", result)
