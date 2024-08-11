#include "orbit_propagator_wrapper.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <yaml-cpp/yaml.h>

OrbitPropagatorWapper::OrbitPropagatorWapper(const std::string &configFilename, double param1, int param2)
    : param1_(param1), param2_(param2)
{
    readConfigFile(configFilename);
    initGlobalVariables();
}

std::vector<double> OrbitPropagatorWapper::propagate(double initial_position, double initial_velocity, double time)
{
    Propagator propagator;
    // Initialize the propagator with necessary parameters
    propagator.setPropOption(forceModelsPropOpt_);
    propagator.printPropOption();
    propagator.initPropagator(initialState_.initialStateVec, epoch_.startMJD, leapSec_, &erpt_, egm_, pJPLEph_);

    double absErr = 1E-6;
    double relErr = 1E-6;
    DynamicModel::stf accMdl = accelerationModel;
    DynamicModel orbFun(accMdl, initialState_.dimState, absErr, relErr);

    double dt = epoch_.timeStep;
    int nTotalSteps = (epoch_.endMJD - epoch_.startMJD) * 86400 / dt + 1;
    std::cout << "total steps:\t" << nTotalSteps << std::endl;

    // Linear spaced times
    Eigen::VectorXd tSec;
    tSec.setLinSpaced(nTotalSteps, 0, (nTotalSteps - 1) * dt);
    Eigen::MatrixXd tableTrajTruth(nTotalSteps, initialState_.dimState + 1);
    tableTrajTruth.col(0) = tSec;
    std::cout << "initial state type:\t" << initialState_.initialStateType << std::endl;

    if (initialState_.initialStateType == "MEE")
    {
        // VectorXd meeSat = initialState_.initialStateVec;
        // tableTrajTruth.row(0).tail(initialState_.dimState) = coe2eci(mee2coe(meeSat), GM_Earth);
    }
    else
    {
        tableTrajTruth.row(0).tail(initialState_.dimState) = initialState_.initialStateVec;
    }

    Eigen::VectorXd propStateVec = initialState_.initialStateVec;
    Timer timer;
    timer.tick();
    for (int k = 1; k < nTotalSteps; k++)
    {
        propStateVec = orbFun(time, time + dt, propStateVec, Eigen::VectorXd::Zero(6));
        time += dt;
        if (initialState_.initialStateType == "MEE")
        {
            // VectorXd meeSat = propStateVec;
            // tableTrajTruth.row(k).tail(initialState_.dimState) = coe2eci(mee2coe(meeSat), GM_Earth);
            // // cout << "meeSat:\t" << coe2eci(mee2coe(meeSat), GM_Earth) << endl;
        }
        else
        {
            tableTrajTruth.row(k).tail(initialState_.dimState) = propStateVec;
        }

        std::cout << "The " << k + 1 << "th time step" << std::endl;
    }
    std::cout << "The total time consumption is:\t" << timer.tock() << std::endl;

    // Convert the result to a std::vector<double>
    std::vector<double> result;
    for (int i = 0; i < tableTrajTruth.rows(); ++i)
    {
        for (int j = 0; j < tableTrajTruth.cols(); ++j)
        {
            result.push_back(tableTrajTruth(i, j));
        }
    }

    // Save the results to a CSV file
    std::vector<std::string> headerTraj({"tSec", "x", "y", "z", "vx", "vy", "vz"});
    std::string propFile = snrInfo_.outDir + "/prop_results.csv";
    EigenCSV::write(tableTrajTruth, headerTraj, propFile);

    return result;
}

void OrbitPropagatorWapper::readConfigFile(const std::string &fileName)
{
    // load file
    YAML::Node config = YAML::LoadFile(fileName);
    YAML::Node parameter;

    // read scenario parameters (required)
    YAML::Node snrParams = config["scenario_parameters"];
    snrInfo_.epoch.startMJD = snrParams["MJD_start"].as<double>();
    snrInfo_.epoch.endMJD = snrParams["MJD_end"].as<double>();
    snrInfo_.epoch.timeStep = snrParams["time_step"].as<double>();
    snrInfo_.outDir = snrParams["output_directory"].as<std::string>();

    // read orbital parameters (required)
    YAML::Node orbitParams = config["initial_orbtial_parameters"];
    int dimState = orbitParams["dim_state"].as<int>();
    initialState_.dimState = dimState;
    std::vector<double> tempVec;
    initialState_.initialStateType = orbitParams["initial_state_type"].as<std::string>();

    // read params as standard vector, convert to eigen vector
    tempVec = orbitParams["initial_state"].as<std::vector<double>>();
    initialState_.initialStateVec = stdVec2EigenVec(tempVec);

    // read propagator settings for filters (optional)
    YAML::Node propFilterSettings = config["propagator_truth_settings"];
    if (parameter = propFilterSettings["earth_gravaity"])
        forceModelsPropOpt_.earth_gravity = parameter.as<bool>();
    if (parameter = propFilterSettings["earth_gravity_model_order"])
        forceModelsPropOpt_.egmAccOrd = parameter.as<int>();
    if (parameter = propFilterSettings["earth_gravity_model_degree"])
        forceModelsPropOpt_.egmAccDeg = parameter.as<int>();
    if (parameter = propFilterSettings["solid_earth_tide"])
        forceModelsPropOpt_.solid_earth_tide = parameter.as<bool>();
    if (parameter = propFilterSettings["ocean_tide_loading"])
        forceModelsPropOpt_.ocean_tide_loading = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_attraction"])
        forceModelsPropOpt_.third_body_attraction = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_sun"])
        forceModelsPropOpt_.third_body_sun = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_moon"])
        forceModelsPropOpt_.third_body_moon = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_planet"])
        forceModelsPropOpt_.third_body_planet = parameter.as<bool>();
    if (parameter = propFilterSettings["relativity_effect"])
        forceModelsPropOpt_.relativity_effect = parameter.as<bool>();
    if (parameter = propFilterSettings["atmospheric_drag"])
        forceModelsPropOpt_.atmospheric_drag = parameter.as<bool>();
    if (parameter = propFilterSettings["solar_radiation_pressure"])
        forceModelsPropOpt_.solar_radiation_pressure = parameter.as<bool>();
    if (parameter = propFilterSettings["thermal_emission"])
        forceModelsPropOpt_.thermal_emission = parameter.as<bool>();
    if (parameter = propFilterSettings["earth_albedo"])
        forceModelsPropOpt_.earth_albedo = parameter.as<bool>();
    if (parameter = propFilterSettings["infrared_radiation"])
        forceModelsPropOpt_.infrared_radiation = parameter.as<bool>();
    if (parameter = propFilterSettings["antenna_thrust"])
        forceModelsPropOpt_.antenna_thrust = parameter.as<bool>();
    if (parameter = propFilterSettings["empirical_acceleration"])
        forceModelsPropOpt_.empirical_acceleration = parameter.as<bool>();
    if (parameter = propFilterSettings["satellite_manoeuvre"])
        forceModelsPropOpt_.satellite_manoeuvre = parameter.as<bool>();
    if (parameter = propFilterSettings["satMass"])
        forceModelsPropOpt_.satMass = parameter.as<double>();
    if (parameter = propFilterSettings["srpArea"])
        forceModelsPropOpt_.srpArea = parameter.as<double>();
    if (parameter = propFilterSettings["srpCoef"])
        forceModelsPropOpt_.srpCoef = parameter.as<double>();
    if (parameter = propFilterSettings["dragArea"])
        forceModelsPropOpt_.dragArea = parameter.as<double>();
    if (parameter = propFilterSettings["dragCoef"])
        forceModelsPropOpt_.dragCoef = parameter.as<double>();

    // read file options for info/data that are relied on
    YAML::Node fileOpt = config["supporting_files"];
    suppFiles_.erpFile = fileOpt["ERP_file"].as<std::string>();
    suppFiles_.grvFile = fileOpt["gravity_file"].as<std::string>();
    suppFiles_.ephFile = fileOpt["ephemeris_file"].as<std::string>();
}

void OrbitPropagatorWapper::initGlobalVariables()
{
    if (!initEGMCoef(suppFiles_.grvFile))
    {
        std::cerr << "Failed to initialize EGM coefficients from " << suppFiles_.grvFile << std::endl;
        return;
    }

    erpt_ = {.n = 0};
    if (!readerp(suppFiles_.erpFile, &erpt_))
    {
        std::cerr << "Failed to read ERP file: " << suppFiles_.erpFile << std::endl;
        return;
    }

    leapSec_ = -getLeapSecond(convertMJD2Time_T(epoch_.startMJD));
    getIERS(epoch_.startMJD, erpt_, iersInstance_);

    const char *ephFile = suppFiles_.ephFile.c_str();
    pJPLEph_ = jpl_init_ephemeris(ephFile, nullptr, nullptr);
    if (!pJPLEph_)
    {
        std::cerr << "Failed to initialize JPL ephemeris from " << suppFiles_.ephFile << std::endl;
        return;
    }
}

Eigen::VectorXd OrbitPropagatorWapper::stdVec2EigenVec(const std::vector<double> &stdVec)
{
    Eigen::VectorXd eigenVec(stdVec.size());
    for (int i = 0; i < stdVec.size(); i++)
    {
        eigenVec(i) = stdVec[i];
    }
    return eigenVec;
}

bool OrbitPropagatorWapper::initEGMCoef(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    int n, m;
    double cmn, smn;

    while (std::getline(file, line))
    {
        std::istringstream buffer(line);
        buffer >> m >> n >> cmn >> smn;

        if (m < GRAVITY_DEG_M && n < GRAVITY_DEG_M)
        {
            egm_.cmn(m, n) = cmn;
            egm_.smn(m, n) = smn;
        }
    }

    return true;
}