import sys
import os


def main(lib_path):
    # Add the directory containing the .so file to the system path
    lib_path = os.path.abspath(lib_path)
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
        return

    # Assuming you have a class named OrbitPropagator in your C++ code
    propagator = orbit_propagator_wrapper.OrbitPropagatorWrapper("yamls/config_orb.yml")
    results = propagator.propagateOrbit()

    # Define the headers and results file name
    headerTraj = ["tSec", "x", "y", "z", "vx", "vy", "vz"]
    resultsFileName = "prop_results_py.csv"
    # Save the results
    propagator.saveResults(results, headerTraj, resultsFileName)

    # print("propagated orbit:", results)


# Run the main function: python3.10 pyscripts/test_orbit_propagator_wrapper.py --lib_path /path/to/your/lib
if __name__ == "__main__":
    # Default lib_path if not provided
    default_lib_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../lib")
    )
    main(default_lib_path)
