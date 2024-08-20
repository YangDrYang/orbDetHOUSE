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
        OrbitPropagatorWapper propagator(configFilename);
        MatrixXd result = propagator.propagate();

        vector<string> headerTraj({"tSec", "x", "y", "z", "vx", "vy", "vz"});
        string propFile = "out/out_prop/prop_results.csv";
        EigenCSV::write(result, headerTraj, propFile);
    }
    catch (const std::exception &e)
    {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }

    return 0;
}