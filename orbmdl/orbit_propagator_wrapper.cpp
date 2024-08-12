#include "orbit_propagator_wrapper.h"

// OrbitPropagatorWapper::OrbitPropagatorWapper(const string &configFilename, double param1, int param2)
//     : param1_(param1), param2_(param2)
// {
//     readConfigFile(configFilename);
//     initGlobalVariables();
// }
OrbitPropagatorWapper::OrbitPropagatorWapper(const string &configFilename)
{
    readConfigFile(configFilename);
    initGlobalVariables();
}

vector<double> OrbitPropagatorWapper::propagate()
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

    double time = 0;
    double dt = epoch_.timeStep;
    int nTotalSteps = (epoch_.endMJD - epoch_.startMJD) * 86400 / dt + 1;
    cout << "total steps:\t" << nTotalSteps << endl;

    // Linear spaced times
    VectorXd tSec;
    tSec.setLinSpaced(nTotalSteps, 0, (nTotalSteps - 1) * dt);
    MatrixXd tableTrajTruth(nTotalSteps, initialState_.dimState + 1);
    tableTrajTruth.col(0) = tSec;
    cout << "initial state type:\t" << initialState_.initialStateType << endl;

    if (initialState_.initialStateType == "MEE")
    {
        // VectorXd meeSat = initialState_.initialStateVec;
        // tableTrajTruth.row(0).tail(initialState_.dimState) = coe2eci(mee2coe(meeSat), GM_Earth);
    }
    else
    {
        tableTrajTruth.row(0).tail(initialState_.dimState) = initialState_.initialStateVec;
    }

    VectorXd propStateVec = initialState_.initialStateVec;
    Timer timer;
    timer.tick();
    for (int k = 1; k < nTotalSteps; k++)
    {
        propStateVec = orbFun(time, time + dt, propStateVec, VectorXd::Zero(6));
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

        cout << "The " << k + 1 << "th time step" << endl;
    }
    cout << "The total time consumption is:\t" << timer.tock() << endl;

    // Convert the result to a vector<double>
    vector<double> result;
    for (int i = 0; i < tableTrajTruth.rows(); ++i)
    {
        for (int j = 0; j < tableTrajTruth.cols(); ++j)
        {
            result.push_back(tableTrajTruth(i, j));
        }
    }

    // // Save the results to a CSV file
    // vector<string> headerTraj({"tSec", "x", "y", "z", "vx", "vy", "vz"});
    // string propFile = snrInfo_.outDir + "/prop_results.csv";
    // EigenCSV::write(tableTrajTruth, headerTraj, propFile);

    return result;
}

void OrbitPropagatorWapper::readConfigFile(const string &fileName)
{
    // load file
    YAML::Node config = YAML::LoadFile(fileName);
    YAML::Node parameter;

    // read scenario parameters (required)
    YAML::Node snrParams = config["scenario_parameters"];
    snrInfo_.epoch.startMJD = snrParams["MJD_start"].as<double>();
    snrInfo_.epoch.endMJD = snrParams["MJD_end"].as<double>();
    snrInfo_.epoch.timeStep = snrParams["time_step"].as<double>();
    snrInfo_.outDir = snrParams["output_directory"].as<string>();

    // read orbital parameters (required)
    YAML::Node orbitParams = config["initial_orbtial_parameters"];
    int dimState = orbitParams["dim_state"].as<int>();
    initialState_.dimState = dimState;
    vector<double> tempVec;
    initialState_.initialStateType = orbitParams["initial_state_type"].as<string>();

    // read params as standard vector, convert to eigen vector
    tempVec = orbitParams["initial_state"].as<vector<double>>();
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
    suppFiles_.erpFile = fileOpt["ERP_file"].as<string>();
    suppFiles_.grvFile = fileOpt["gravity_file"].as<string>();
    suppFiles_.ephFile = fileOpt["ephemeris_file"].as<string>();
}

void OrbitPropagatorWapper::initGlobalVariables()
{
    if (!initEGMCoef(suppFiles_.grvFile))
    {
        cerr << "Failed to initialize EGM coefficients from " << suppFiles_.grvFile << endl;
        return;
    }

    erpt_ = {.n = 0};
    if (!readerp(suppFiles_.erpFile, &erpt_))
    {
        cerr << "Failed to read ERP file: " << suppFiles_.erpFile << endl;
        return;
    }

    leapSec_ = -getLeapSecond(convertMJD2Time_T(epoch_.startMJD));
    getIERS(epoch_.startMJD, erpt_, iersInstance_);

    const char *ephFile = suppFiles_.ephFile.c_str();
    pJPLEph_ = jpl_init_ephemeris(ephFile, nullptr, nullptr);
    if (!pJPLEph_)
    {
        cerr << "Failed to initialize JPL ephemeris from " << suppFiles_.ephFile << endl;
        return;
    }
}

VectorXd OrbitPropagatorWapper::stdVec2EigenVec(const vector<double> &stdVec)
{
    VectorXd eigenVec(stdVec.size());
    for (int i = 0; i < stdVec.size(); i++)
    {
        eigenVec(i) = stdVec[i];
    }
    return eigenVec;
}

bool OrbitPropagatorWapper::initEGMCoef(const string &filename)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed to open file: " << filename << endl;
        return false;
    }

    string line;
    int n, m;
    double cmn, smn;

    while (getline(file, line))
    {
        istringstream buffer(line);
        buffer >> m >> n >> cmn >> smn;

        if (m < GRAVITY_DEG_M && n < GRAVITY_DEG_M)
        {
            egm_.cmn(m, n) = cmn;
            egm_.smn(m, n) = smn;
        }
    }

    return true;
}