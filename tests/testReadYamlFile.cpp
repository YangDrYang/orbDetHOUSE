#include "testReadYamlFile.hpp"
using namespace Eigen;
using namespace std;

// propagator variables
// erp_t erpt;
VectorXd initialState(6), groundStation(6);
double leapSec;
// IERS iersInstance;
ForceModels forceModelsOpt = {};
// EGMCoef egm;
// void *pJPLEph;
// Propagator orbitProp;
struct EpochInfo epoch;

void readConfigFile(string fileName, ForceModels &options, struct EpochInfo &epoch, VectorXd &initialState,
                    VectorXd &groundStation, struct Filters &filters, int &numTrials, struct Errors &errorStatistics,
                    string &iniatialStateType)
{
    // load file
    YAML::Node config = YAML::LoadFile(fileName);
    YAML::Node parameter;

    // read filter options (required)
    YAML::Node filterOpts = config["filter_options"];
    filters.house = filterOpts["HOUSE"].as<bool>();
    filters.ukf = filterOpts["UKF"].as<bool>();
    filters.cut4 = filterOpts["CUT4"].as<bool>();
    filters.cut6 = filterOpts["CUT6"].as<bool>();
    numTrials = filterOpts["num_trials"].as<int>();

    // read orbital parameters (required)
    YAML::Node orbitParams = config["initial_orbtial_parameters"];
    std::vector<double> tempVec;
    double test = orbitParams["MJD_start"].as<double>();
    epoch.startMJD = orbitParams["MJD_start"].as<double>();
    epoch.endMJD = orbitParams["MJD_end"].as<double>();
    epoch.timeStep = orbitParams["time_step"].as<double>();
    errorStatistics.azimuthErr = orbitParams["elevation_error"].as<double>();
    errorStatistics.elevationErr = orbitParams["azimuth_error"].as<double>();
    errorStatistics.rangeErr = orbitParams["range_error"].as<double>();
    errorStatistics.rangeRateErr = orbitParams["range_rate_error"].as<double>();
    iniatialStateType = orbitParams["initial_state_type"].as<string>();

    // read params as standard vector, convert to eigen vector
    tempVec = orbitParams["initial_state"].as<std::vector<double>>();
    initialState = stdVec2EigenVec(tempVec);
    tempVec = orbitParams["ground_station"].as<std::vector<double>>();
    groundStation = stdVec2EigenVec(tempVec);

    // read propagator settings (optional)
    YAML::Node propSettings = config["propagator_truth_settings"];
    if (parameter = propSettings["earth_gravaity"])
        options.earth_gravity = parameter.as<bool>();
    if (parameter = propSettings["solid_earth_tide"])
        options.solid_earth_tide = parameter.as<bool>();
    if (parameter = propSettings["ocean_tide_loading"])
        options.ocean_tide_loading = parameter.as<bool>();
    if (parameter = propSettings["third_body_attraction"])
        options.third_body_attraction = parameter.as<bool>();
    if (parameter = propSettings["third_body_sun"])
        options.third_body_sun = parameter.as<bool>();
    if (parameter = propSettings["third_body_moon"])
        options.third_body_moon = parameter.as<bool>();
    if (parameter = propSettings["third_body_planet"])
        options.third_body_planet = parameter.as<bool>();
    if (parameter = propSettings["relativity_effect"])
        options.relativity_effect = parameter.as<bool>();
    if (parameter = propSettings["solar_radiation_pressure"])
        options.solar_radiation_pressure = parameter.as<bool>();
    if (parameter = propSettings["thermal_emission"])
        options.thermal_emission = parameter.as<bool>();
    if (parameter = propSettings["earth_albedo"])
        options.earth_albedo = parameter.as<bool>();
    if (parameter = propSettings["infrared_radiation"])
        options.infrared_radiation = parameter.as<bool>();
    if (parameter = propSettings["antenna_thrust"])
        options.antenna_thrust = parameter.as<bool>();
    if (parameter = propSettings["empirical_acceleration"])
        options.empirical_acceleration = parameter.as<bool>();
    if (parameter = propSettings["satellite_manoeuvre"])
        options.satellite_manoeuvre = parameter.as<bool>();
    if (parameter = propSettings["satMass"])
        options.satMass = parameter.as<double>();
    if (parameter = propSettings["srpArea"])
        options.srpArea = parameter.as<double>();
    if (parameter = propSettings["srpCoef"])
        options.srpCoef = parameter.as<double>();
}

VectorXd stdVec2EigenVec(const std::vector<double> &stdVec)
{
    VectorXd eigenVec(stdVec.size());
    for (int i = 0; i < 6; i++)
    {
        eigenVec(i) = stdVec[i];
    }
    return eigenVec;
}

int main(int argc, char *argv[])
{
    string configFilename;
    // input config .yaml file
    switch (argc)
    {
    case 1: // Read default file
        configFilename = DEFAULT_CONFIG_FILENAME;
        break;
    case 2: // Read supplied file
        configFilename = argv[1];
        break;
    default: // wrong number of args
        cerr << "Accepts up to 1 argument (.yaml input file).\nIf no argument given, the default config file will be read (config.yaml)" << endl;
        exit(1);
        break;
    }
    cout << "Reading configuration from file: " << configFilename << endl;

    struct Errors errorStatistics;
    struct Filters filters;
    int numTrials;
    string initialStateType;
    // read parameter/settings from config file
    readConfigFile(configFilename, forceModelsOpt, epoch, initialState, groundStation, filters, numTrials, errorStatistics, initialStateType);

    cout << "earth_gravity: " << forceModelsOpt.earth_gravity << endl;
    cout << "solid_earth_tide: " << forceModelsOpt.solid_earth_tide << endl;
    cout << "satMass: " << forceModelsOpt.satMass << endl;
    cout << "srpArea: " << forceModelsOpt.srpArea << endl;
    cout << "srpCoef: " << forceModelsOpt.srpCoef << endl;

    return 0;
}