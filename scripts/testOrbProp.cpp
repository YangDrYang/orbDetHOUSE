#include "testOrbDet.hpp"

// propagator variables
erp_t erpt;
double leapSec;
IERS iersInstance;
EGMCoef egm;
void *pJPLEph;
Propagator orbitProp;
struct EpochInfo epoch;

VectorXd simpleAccerationModel(double t, const VectorXd &X, const VectorXd &fd);
VectorXd simpleAccerationModel(double t, const VectorXd &X, const VectorXd &fd)
{
    VectorXd Xf(6);
    Vector3d r, v, a;
    r = X.head(3);
    v = X.tail(3);

    // Earth's central body gravity acceleration
    a = -MU / pow(r.norm(), 3) * r;

    Xf.head(3) = v;
    Xf.tail(3) = a;
    return Xf;
}

VectorXd accelerationModel(double tSec, const VectorXd &X, const VectorXd &fd)
{
    VectorXd Xf(6);

    Matrix3d mECI2ECEF = Matrix3d::Identity();
    Matrix3d mdECI2ECEF = Matrix3d::Identity();

    // get leap seconds from the table
    leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD + tSec / 86400));
    // set up the IERS instance
    getIERS(epoch.startMJD + tSec / 86400, erpt, iersInstance);

    eci2ecef_sofa(epoch.startMJD + tSec / 86400, iersInstance, mECI2ECEF, mdECI2ECEF);

    Vector3d acceleration;
    orbitProp.updPropagator(epoch.startMJD + tSec / 86400, leapSec, &erpt);

    if (abs(X(1)) < 1 && abs(X(2)) < 1) // modified equnoctial elements
    {
        // // calculate acceleration
        // Xf = orbitProp.calculateTimeDerivativeMEE(X, mECI2ECEF);
    }
    else // Cartesian elements
    {
        // calculate acceleration
        acceleration = orbitProp.calculateAcceleration(X.head(3), X.tail(3), mECI2ECEF);
        // set state vector
        Xf.head(3) = X.tail(3);
        Xf.tail(3) = acceleration;
        // cout << "acceleration:\t" << acceleration << endl;
    }

    return Xf;
}

void readConfigFile(string fileName, ForceModels &optProp, struct ScenarioInfo &snrInfo, struct InitialState &initialState,
                    struct MeasModel &measMdl, struct Filters &filters, struct FileInfo &suppFiles)
{
    // load file
    YAML::Node config = YAML::LoadFile(fileName);
    YAML::Node parameter;

    // read scenario parameters (required)
    YAML::Node snrParams = config["scenario_parameters"];
    snrInfo.epoch.startMJD = snrParams["MJD_start"].as<double>();
    snrInfo.epoch.endMJD = snrParams["MJD_end"].as<double>();
    snrInfo.epoch.timeStep = snrParams["time_step"].as<double>();
    snrInfo.outDir = snrParams["output_directory"].as<string>();

    // read orbital parameters (required)
    YAML::Node orbitParams = config["initial_orbtial_parameters"];
    int dimState = orbitParams["dim_state"].as<int>();
    initialState.dimState = dimState;
    vector<double> tempVec;
    initialState.initialStateType = orbitParams["initial_state_type"].as<string>();

    // read params as standard vector, convert to eigen vector
    tempVec = orbitParams["initial_state"].as<vector<double>>();
    initialState.initialStateVec = stdVec2EigenVec(tempVec);

    // read propagator settings for filters (optional)
    YAML::Node propFilterSettings = config["propagator_truth_settings"];
    if (parameter = propFilterSettings["earth_gravaity"])
        optProp.earth_gravity = parameter.as<bool>();
    if (parameter = propFilterSettings["earth_gravity_model_order"])
        optProp.egmAccOrd = parameter.as<int>();
    if (parameter = propFilterSettings["earth_gravity_model_degree"])
        optProp.egmAccDeg = parameter.as<int>();
    if (parameter = propFilterSettings["solid_earth_tide"])
        optProp.solid_earth_tide = parameter.as<bool>();
    if (parameter = propFilterSettings["ocean_tide_loading"])
        optProp.ocean_tide_loading = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_attraction"])
        optProp.third_body_attraction = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_sun"])
        optProp.third_body_sun = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_moon"])
        optProp.third_body_moon = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_planet"])
        optProp.third_body_planet = parameter.as<bool>();
    if (parameter = propFilterSettings["relativity_effect"])
        optProp.relativity_effect = parameter.as<bool>();
    if (parameter = propFilterSettings["atmospheric_drag"])
        optProp.atmospheric_drag = parameter.as<bool>();
    if (parameter = propFilterSettings["solar_radiation_pressure"])
        optProp.solar_radiation_pressure = parameter.as<bool>();
    if (parameter = propFilterSettings["thermal_emission"])
        optProp.thermal_emission = parameter.as<bool>();
    if (parameter = propFilterSettings["earth_albedo"])
        optProp.earth_albedo = parameter.as<bool>();
    if (parameter = propFilterSettings["infrared_radiation"])
        optProp.infrared_radiation = parameter.as<bool>();
    if (parameter = propFilterSettings["antenna_thrust"])
        optProp.antenna_thrust = parameter.as<bool>();
    if (parameter = propFilterSettings["empirical_acceleration"])
        optProp.empirical_acceleration = parameter.as<bool>();
    if (parameter = propFilterSettings["satellite_manoeuvre"])
        optProp.satellite_manoeuvre = parameter.as<bool>();
    if (parameter = propFilterSettings["satMass"])
        optProp.satMass = parameter.as<double>();
    if (parameter = propFilterSettings["srpArea"])
        optProp.srpArea = parameter.as<double>();
    if (parameter = propFilterSettings["srpCoef"])
        optProp.srpCoef = parameter.as<double>();
    if (parameter = propFilterSettings["dragArea"])
        optProp.dragArea = parameter.as<double>();
    if (parameter = propFilterSettings["dragCoef"])
        optProp.dragCoef = parameter.as<double>();

    // read file options for info/data that are relied on
    YAML::Node fileOpt = config["supporting_files"];
    suppFiles.erpFile = fileOpt["ERP_file"].as<string>();
    suppFiles.grvFile = fileOpt["gravity_file"].as<string>();
    suppFiles.ephFile = fileOpt["ephemeris_file"].as<string>();
}

VectorXd stdVec2EigenVec(const vector<double> &stdVec)
{
    VectorXd eigenVec(stdVec.size());
    for (int i = 0; i < 6; i++)
    {
        eigenVec(i) = stdVec[i];
    }
    return eigenVec;
}

bool initEGMCoef(const string &filename)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed to open file: " << filename << endl;
        return false; // Indicate failure to open the file
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
            egm.cmn(m, n) = cmn;
            egm.smn(m, n) = smn;
        }
    }

    return true; // Indicate success
}

void initGlobalVariables(VectorXd &initialStateVec, string stateType, struct FileInfo &suppFiles)
{
    if (!initEGMCoef(suppFiles.grvFile))
    {
        cerr << "Failed to initialize EGM coefficients from " << suppFiles.grvFile << endl;
        return; // Assuming initEGMCoef returns a bool indicating success/failure
    }

    // cout << suppFiles.erpFile << endl;
    erpt.n = 0;
    if (!readerp(suppFiles.erpFile, &erpt))
    {
        cerr << "Failed to read ERP file: " << suppFiles.erpFile << endl;
        return; // Assuming readerp returns a bool indicating success/failure
    }
    // cout << "erpt mjd\t" << erpt.data->mjd << endl;

    // get leap seconds from the table
    leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD));
    // set up the IERS instance
    getIERS(epoch.startMJD, erpt, iersInstance);

    const char *ephFile = suppFiles.ephFile.c_str();
    pJPLEph = jpl_init_ephemeris(ephFile, nullptr, nullptr);
    if (!pJPLEph)
    {
        cerr << "Failed to initialize JPL ephemeris from " << suppFiles.ephFile << endl;
        return; // Check if jpl_init_ephemeris returns nullptr on failure
    }

    // todo: fix the error zsh: abort      bin/scripts/testOrbProp yamls/config_orb.yml
    //  cout << "initial state type:\t" << stateType << endl;
    //  if (stateType == "ECEF")
    //  {
    //      cout << "Convert state from ECEF to ECI\n";
    //      VectorXd rvECI = VectorXd::Zero(6);
    //      ecef2eciVec_sofa(epoch.startMJD, iersInstance, initialStateVec, rvECI);
    //      initialStateVec = rvECI;
    //      cout << "ECI initial state:\t" << initialStateVec << endl;
    //  }
    //  else if (stateType == "MEE")
    //  {
    //      cout << "Convert state from ECI to MEE\n";
    //      VectorXd rvECI = VectorXd::Zero(6);
    //      rvECI = initialStateVec;
    //      initialStateVec = coe2mee(eci2coe(rvECI, GM_Earth));
    //      cout << "MEE initial state:\t" << initialStateVec << endl;
    //  }
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

    struct ForceModels forceModelsPropOpt;
    struct ScenarioInfo snrInfo;
    struct MeasModel measMdl;
    struct Filters filters;
    struct InitialState initialState;
    struct FileInfo suppFiles;
    // read parameter/settings from config file
    readConfigFile(configFilename, forceModelsPropOpt, snrInfo, initialState, measMdl, filters, suppFiles);
    epoch = snrInfo.epoch;
    const int dimState = initialState.dimState;
    string initialStateType = initialState.initialStateType;
    VectorXd initialStateVec = initialState.initialStateVec;

    // initialise
    initGlobalVariables(initialStateVec, initialStateType, suppFiles);

    // setup orbit propagator
    orbitProp.setPropOption(forceModelsPropOpt);
    orbitProp.printPropOption();
    orbitProp.initPropagator(initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);
    double absErr = 1E-6;
    double relErr = 1E-6;
    DynamicModel::stf accMdl = accelerationModel;
    DynamicModel orbFun(accMdl, dimState, absErr, relErr);

    double time = 0, dt = epoch.timeStep;
    int nTotalSteps = (epoch.endMJD - epoch.startMJD) * 86400 / dt + 1;
    cout << "total steps:\t" << nTotalSteps << endl;
    // linear spaced times
    VectorXd tSec;
    tSec.setLinSpaced(nTotalSteps, 0, (nTotalSteps - 1) * dt);
    MatrixXd tableTrajTruth(nTotalSteps, dimState + 1);
    tableTrajTruth.col(0) = tSec;
    cout << "initial state type:\t" << initialStateType << endl;
    if (initialStateType == "MEE")
    {
        // VectorXd meeSat = initialStateVec;
        // tableTrajTruth.row(0).tail(dimState) = coe2eci(mee2coe(meeSat), GM_Earth);
    }
    else
    {
        tableTrajTruth.row(0).tail(dimState) = initialStateVec;
    }

    VectorXd propStateVec = initialStateVec;
    Timer timer;
    timer.tick();
    for (int k = 1; k < nTotalSteps; k++)
    {
        propStateVec = orbFun(time, time + dt, propStateVec, VectorXd::Zero(6));
        time += dt;
        if (initialStateType == "MEE")
        {
            // VectorXd meeSat = propStateVec;
            // tableTrajTruth.row(k).tail(dimState) = coe2eci(mee2coe(meeSat), GM_Earth);
            // // cout << "meeSat:\t" << coe2eci(mee2coe(meeSat), GM_Earth) << endl;
        }
        else
        {
            tableTrajTruth.row(k).tail(dimState) = propStateVec;
        }

        cout << "The " << k + 1 << "th time step" << endl;
    }
    cout << "The total time consumption is:\t" << timer.tock() << endl;
    // header for the saved file
    vector<string> headerTraj({"tSec", "x", "y", "z", "vx", "vy", "vz"});
    string propFile = snrInfo.outDir + "/prop_results.csv";
    EigenCSV::write(tableTrajTruth, headerTraj, propFile);
}