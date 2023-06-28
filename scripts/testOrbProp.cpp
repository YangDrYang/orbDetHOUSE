#include "testOrbDet.hpp"

#define R_EARTH 6371E3
#define DEG M_PI / 180
#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)
#define MAXLEAPS 18
const double leaps[MAXLEAPS + 1][7] =
    {
        /* leap seconds (y,m,d,h,m,s,utc-gpst) */
        {2017, 1, 1, 0, 0, 0, -18},
        {2015, 7, 1, 0, 0, 0, -17},
        {2012, 7, 1, 0, 0, 0, -16},
        {2009, 1, 1, 0, 0, 0, -15},
        {2006, 1, 1, 0, 0, 0, -14},
        {1999, 1, 1, 0, 0, 0, -13},
        {1997, 7, 1, 0, 0, 0, -12},
        {1996, 1, 1, 0, 0, 0, -11},
        {1994, 7, 1, 0, 0, 0, -10},
        {1993, 7, 1, 0, 0, 0, -9},
        {1992, 7, 1, 0, 0, 0, -8},
        {1991, 1, 1, 0, 0, 0, -7},
        {1990, 1, 1, 0, 0, 0, -6},
        {1988, 1, 1, 0, 0, 0, -5},
        {1985, 7, 1, 0, 0, 0, -4},
        {1983, 7, 1, 0, 0, 0, -3},
        {1982, 7, 1, 0, 0, 0, -2},
        {1981, 7, 1, 0, 0, 0, -1},
        {0}};

// propagator variables
erp_t erpt;
double leapSec;
IERS iersInstance;
EGMCoef egm;
void *pJPLEph;
Propagator orbitProp;
struct EpochInfo epoch;

time_t convertMJD2Time_T(double mjd)
{
    // Convert MJD to Unix timestamp
    double unix_time = (mjd - 40587) * 86400.0;

    // Convert Unix timestamp to time_t
    time_t t = static_cast<time_t>(unix_time);

    return t;
}

double getLeapSecond(time_t t)
{
    // convert to tm
    tm tm = *gmtime(&t);
    int year = tm.tm_year + 1900;

    // find the latest leap second that is earlier than t
    int i;
    for (i = 0; i < MAXLEAPS; i++)
    {
        if (leaps[i][0] < year || (leaps[i][0] == year && leaps[i][1] < tm.tm_mon + 1) ||
            (leaps[i][0] == year && leaps[i][1] == tm.tm_mon + 1 && leaps[i][2] <= tm.tm_mday))
        {
            break;
        }
    }

    // return the leap second value
    return leaps[i][6];
}

void getIERS(double mjd)
{
    // get leap seconds from the table
    // double leapSec = -getLeapSecond(convertMJD2Time_T(mjd));

    double erpv[4] = {};
    geterp_from_utc(&erpt, leapSec, mjd, erpv);

    double dUT1_UTC = erpv[2];
    double dUTC_TAI = -(19 + leapSec);
    double xp = erpv[0];
    double yp = erpv[1];
    double lod = erpv[3];

    iersInstance.Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
}

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
    Vector3d r;
    Vector3d v;
    VectorXd rvECI(6);
    rvECI = X;

    Matrix3d mECI2ECEF = Matrix3d::Identity();
    Matrix3d mdECI2ECEF = Matrix3d::Identity();

    // get leap seconds from the table
    leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD + tSec / 86400));
    // set up the IERS instance
    getIERS(epoch.startMJD + tSec / 86400);

    eci2ecef_sofa(epoch.startMJD + tSec / 86400, iersInstance, mECI2ECEF, mdECI2ECEF);

    Vector3d acceleration;
    orbitProp.updPropagator(epoch.startMJD + tSec / 86400, leapSec, &erpt);

    // calculate acceleration
    acceleration = orbitProp.calculateAcceleration(X.head(3), X.tail(3), mECI2ECEF);

    // set state vector
    Xf.head(3) = X.tail(3);
    Xf.tail(3) = acceleration;

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

void initEGMCoef(string filename)
{
    ifstream file(filename);
    string line;
    int n, m;
    double cmn, smn;

    while (getline(file, line))
    {
        // line structure:
        // m        n       Cnm     Snm     0      0
        istringstream buffer(line);
        buffer >> m >> n >> cmn >> smn;

        // gravity model degree
        if (m < GRAVITY_DEG_M && n < GRAVITY_DEG_M)
        {
            egm.cmn(m, n) = cmn;
            egm.smn(m, n) = smn;
        }
    }
}

void initGlobalVariables(VectorXd &initialStateVec, string stateType, struct FileInfo &suppFiles)
{
    initEGMCoef(suppFiles.grvFile);

    erpt = {.n = 0};
    // cout << suppFiles.erpFile << endl;
    readerp(suppFiles.erpFile, &erpt);
    // cout << "erpt mjd\t" << erpt.data->mjd << endl;

    // get leap seconds from the table
    leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD));
    // set up the IERS instance
    getIERS(epoch.startMJD);

    const char *ephFile = suppFiles.ephFile.c_str();
    pJPLEph = jpl_init_ephemeris(ephFile, nullptr, nullptr);

    VectorXd rvECI = VectorXd::Zero(6);
    string ecefTag = "ECEF";
    if (stateType == ecefTag)
    {
        cout << "Converted state from ECEF to ECI\n";
        VectorXd rvECI = VectorXd::Zero(6);
        ecef2eciVec_sofa(epoch.startMJD, iersInstance, initialStateVec, rvECI);
        initialStateVec = rvECI;
    }
    else
    {
        rvECI = initialStateVec;
    }
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
    MatrixXd initialCov = initialState.initialCovarianceMat;

    // initialise
    initGlobalVariables(initialStateVec, initialStateType, suppFiles);
    // // get leap seconds from the table
    // double leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD));

    // setup orbit propagator
    orbitProp.setPropOption(forceModelsPropOpt);
    orbitProp.printPropOption();
    orbitProp.initPropagator(initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);
    double absErr = 1E-6;
    double relErr = 1E-6;
    DynamicModel::stf g = accelerationModel;
    DynamicModel f(g, dimState, absErr, relErr);

    VectorXd propStateVec = initialStateVec;
    double time = 0, dt = epoch.timeStep;
    int nTotalSteps = (epoch.endMJD - epoch.startMJD) * 86400 / dt + 1;
    cout << "total steps:\t" << nTotalSteps << endl;
    // linear spaced times
    VectorXd tSec;
    tSec.setLinSpaced(nTotalSteps, 0, (nTotalSteps - 1) * dt);
    MatrixXd tableTrajTruth(nTotalSteps, dimState + 1);
    tableTrajTruth.col(0) = tSec;
    tableTrajTruth.row(0).tail(dimState) = initialStateVec;

    Timer timer;
    timer.tick();
    for (int k = 0; k < nTotalSteps - 1; k++)
    {
        propStateVec = f(time, time + dt, propStateVec, Vector3d::Zero());
        time += dt;
        tableTrajTruth.row(k + 1).tail(dimState) = propStateVec;
        // cout << "The " << k + 1 << "th time step" << endl;
    }
    cout << "The total time consumption is:\t" << timer.tock() << endl;
    // header for the saved file
    vector<string> headerTraj({"tSec", "x", "y", "z", "vx", "vy", "vz"});
    string propFile = snrInfo.outDir + "/prop_results.csv";
    EigenCSV::write(tableTrajTruth, headerTraj, propFile);
}
