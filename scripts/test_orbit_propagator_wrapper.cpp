#include "orbit_propagator_wrapper.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: " << argv[0] << " <config_file>" << endl;
        return 1;
    }

    string configFilename = argv[1];

    try
    {
        OrbitPropagatorWrapper propagator(configFilename);
        MatrixXd results = propagator.propagateOrbit();
        // Save the results to a CSV file
        vector<string> headerTraj({"tSec", "x", "y", "z", "vx", "vy", "vz"});
        string resultsFileName = "/prop_results.csv";
        propagator.saveResults(results, headerTraj, resultsFileName);
    }
    catch (const exception &e)
    {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }

    return 0;
}