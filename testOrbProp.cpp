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
ForceModels forceModelsTruthOpt = {}, forceModelsFilterOpt = {};
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
    // // get leap seconds from the table
    // double leapSec = getLeapSecond(convertMJD2Time_T(mjd));

    double erpv[4] = {};
    geterp_from_utc(&erpt, leapSec, mjd, erpv);

    double dUT1_UTC = erpv[2];
    double dUTC_TAI = -19 + leapSec;
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

Vector2d simpleMeasurementModel(double tSec, const VectorXd &satECI, const VectorXd &stnECEF)
{
    Vector2d z;
    Vector3d p;
    VectorXd rsECEF(6); // = VectorXd::Zero(6);
    VectorXd stnECI = VectorXd::Zero(6);
    VectorXd stnECEF_ = stnECEF; // define a new variable that is not a const vector

    ecef2eciVec_sofa(epoch.startMJD + tSec / 86400, iersInstance, stnECEF_, stnECI);
    // end transformation code
    p = satECI.head(3) - stnECI.head(3);

    if (p.dot(stnECI.head(3)) >= 0)
    {
        // range
        // azimuth angle
        z(0) = atan2(p(1), p(0));
        // elevation angle
        z(1) = asin(p(2) / p.norm());
    }
    // not visible
    else
    {
        z(0) = NO_MEASUREMENT;
        z(1) = NO_MEASUREMENT;
    }
    return z;
}

Vector4d measurementModel(double tSec, const VectorXd &satECI, const VectorXd &stnECEF)
{
    Vector4d z;
    Vector3d p, r, v;
    VectorXd stnECI = VectorXd::Zero(6);
    VectorXd stnECEF_ = stnECEF; // define a new variable that is not a const vector

    // get leap seconds from the table
    leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD + tSec / 86400));
    // set up the IERS instance
    getIERS(epoch.startMJD + tSec / 86400);

    ecef2eciVec_sofa(epoch.startMJD + tSec / 86400, iersInstance, stnECEF_, stnECI);
    // cout << "epoch" << tSec << endl;
    // cout << "ground station in ECI" << stnECI.transpose() << endl;

    // end transformation code
    p = satECI.head(3) - stnECI.head(3);
    v = satECI.tail(3) - stnECI.tail(3);

    if (p.dot(stnECI.head(3)) >= 0)
    {
        // azimuth angle
        z(0) = atan2(p(1), p(0));
        // elevation angle
        z(1) = asin(p(2) / p.norm());
        // range
        z(2) = p.norm();
        // range rate
        z(3) = p.dot(v) / p.norm();
    }
    // not visible
    else
    {
        z(0) = NO_MEASUREMENT;
        z(1) = NO_MEASUREMENT;
        z(2) = NO_MEASUREMENT;
        z(3) = NO_MEASUREMENT;
    }

    // // azimuth angle
    // z(0) = atan2(p(1), p(0));
    // // elevation angle
    // z(1) = asin(p(2) / p.norm());
    // // range
    // z(2) = p.norm();
    // // range rate
    // z(3) = p.dot(v) / p.norm();

    return z;
}

// // Function to generate a noise matrix for a window of epochs
MatrixXd generateNoiseMatrix(int seed, int nEpoch, const MatrixXd &covariance)
{
    // Create a random number generator with the given seed
    mt19937_64 gen(seed);

    // Create a normal distribution with mean 0 and covariance matrix
    MatrixXd L = covariance.llt().matrixL();
    normal_distribution<double> dist(0.0, 1.0);

    // Generate a noise matrix with the specified number of epochs
    MatrixXd noise(covariance.rows(), nEpoch);
    for (int i = 0; i < nEpoch; i++)
    {
        for (int j = 0; j < covariance.rows(); j++)
        {
            noise(j, i) = dist(gen);
        }
        noise.col(i) = L * noise.col(i);
    }

    return noise;
}

MatrixXd generateTrueResults(DynamicModel &f, struct EpochInfo epoch, VectorXd initialStateVec)
{
    // create results matrix
    int length = floor((epoch.endMJD - epoch.startMJD) / (epoch.timeStep / 86400.0)) + 1;
    MatrixXd results(length, 7);

    double time = 0;
    double step = epoch.timeStep;

    Vector3d W = Vector3d::Zero();
    VectorXd X = initialStateVec;
    VectorXd result(7);

    // Save initial state
    result << time, X;
    results.row(0) = result;
    // TODO: remove variation is forces model (the w vector)
    for (int i = 1; i < length; i++)
    {
        X = f(time, time + step, X, W);
        time += step;
        result << time, X;
        results.row(i) = result;
    }

    return results;
}

void readConfigFile(string fileName, ForceModels &optTruth, ForceModels &optFilter, struct SimInfo &simInfo, struct InitialState &initialState,
                    struct MeasModel &measMdl, struct Filters &filters)
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
    filters.numTrials = filterOpts["num_trials"].as<int>();

    // read simulation parameters (required)
    YAML::Node simParams = config["simulation_parameters"];
    simInfo.epoch.startMJD = simParams["MJD_start"].as<double>();
    simInfo.epoch.endMJD = simParams["MJD_end"].as<double>();
    simInfo.epoch.timeStep = simParams["time_step"].as<double>();
    simInfo.epoch.timePass = simParams["time_pass"].as<double>();
    simInfo.file.outDir = simParams["output_directory"].as<string>();

    // read orbital parameters (required)
    YAML::Node orbitParams = config["initial_orbtial_parameters"];
    int dimState = orbitParams["dim_state"].as<int>();
    initialState.dimState = dimState;
    vector<double> tempVec;
    initialState.initialStateType = orbitParams["initial_state_type"].as<string>();
    // read params as standard vector, convert to eigen vector
    tempVec = orbitParams["initial_state"].as<vector<double>>();
    initialState.initialStateVec = stdVec2EigenVec(tempVec);
    // Read matrix from YAML file
    MatrixXd tempMat = MatrixXd::Zero(dimState, dimState);
    const YAML::Node &covInitState = orbitParams["initial_covariance"];
    for (int i = 0; i < dimState; ++i)
    {
        const YAML::Node &row = covInitState[i];
        for (int j = 0; j < dimState; ++j)
        {
            tempMat(i, j) = row[j].as<double>();
        }
    }
    initialState.initialCovarianceMat = tempMat;

    tempMat = MatrixXd::Zero(dimState, dimState);
    const YAML::Node &covProNoise = orbitParams["process_noise_covariance"];
    for (int i = 0; i < dimState; ++i)
    {
        const YAML::Node &row = covProNoise[i];
        for (int j = 0; j < dimState; ++j)
        {
            tempMat(i, j) = row[j].as<double>();
        }
    }
    initialState.processNoiseCovarianceMat = tempMat;

    // read measurement characristics (noise standard deviations)
    YAML::Node measParams = config["measurement_parameters"];
    tempVec = measParams["ground_station"].as<vector<double>>();
    measMdl.groundStation = stdVec2EigenVec(tempVec);
    measMdl.dimMeas = measParams["dim_meas"].as<int>();
    measMdl.errorStd.azimuthErr = measParams["elevation_error"].as<double>();
    measMdl.errorStd.elevationErr = measParams["azimuth_error"].as<double>();
    measMdl.errorStd.rangeErr = measParams["range_error"].as<double>();
    measMdl.errorStd.rangeRateErr = measParams["range_rate_error"].as<double>();

    // read propagator settings (optional)
    YAML::Node propTruthSettings = config["propagator_truth_settings"];
    if (parameter = propTruthSettings["earth_gravaity"])
        optTruth.earth_gravity = parameter.as<bool>();
    if (parameter = propTruthSettings["solid_earth_tide"])
        optTruth.solid_earth_tide = parameter.as<bool>();
    if (parameter = propTruthSettings["ocean_tide_loading"])
        optTruth.ocean_tide_loading = parameter.as<bool>();
    if (parameter = propTruthSettings["third_body_attraction"])
        optTruth.third_body_attraction = parameter.as<bool>();
    if (parameter = propTruthSettings["third_body_sun"])
        optTruth.third_body_sun = parameter.as<bool>();
    if (parameter = propTruthSettings["third_body_moon"])
        optTruth.third_body_moon = parameter.as<bool>();
    if (parameter = propTruthSettings["third_body_planet"])
        optTruth.third_body_planet = parameter.as<bool>();
    if (parameter = propTruthSettings["relativity_effect"])
        optTruth.relativity_effect = parameter.as<bool>();
    if (parameter = propTruthSettings["atmospheric_drag"])
        optTruth.atmospheric_drag = parameter.as<bool>();
    if (parameter = propTruthSettings["solar_radiation_pressure"])
        optTruth.solar_radiation_pressure = parameter.as<bool>();
    if (parameter = propTruthSettings["thermal_emission"])
        optTruth.thermal_emission = parameter.as<bool>();
    if (parameter = propTruthSettings["earth_albedo"])
        optTruth.earth_albedo = parameter.as<bool>();
    if (parameter = propTruthSettings["infrared_radiation"])
        optTruth.infrared_radiation = parameter.as<bool>();
    if (parameter = propTruthSettings["antenna_thrust"])
        optTruth.antenna_thrust = parameter.as<bool>();
    if (parameter = propTruthSettings["empirical_acceleration"])
        optTruth.empirical_acceleration = parameter.as<bool>();
    if (parameter = propTruthSettings["satellite_manoeuvre"])
        optTruth.satellite_manoeuvre = parameter.as<bool>();
    if (parameter = propTruthSettings["satMass"])
        optTruth.satMass = parameter.as<double>();
    if (parameter = propTruthSettings["srpArea"])
        optTruth.srpArea = parameter.as<double>();
    if (parameter = propTruthSettings["srpCoef"])
        optTruth.srpCoef = parameter.as<double>();

    // read propagator settings for filters (optional)
    YAML::Node propFilterSettings = config["propagator_filter_settings"];
    if (parameter = propFilterSettings["earth_gravaity"])
        optFilter.earth_gravity = parameter.as<bool>();
    if (parameter = propFilterSettings["solid_earth_tide"])
        optFilter.solid_earth_tide = parameter.as<bool>();
    if (parameter = propFilterSettings["ocean_tide_loading"])
        optFilter.ocean_tide_loading = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_attraction"])
        optFilter.third_body_attraction = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_sun"])
        optFilter.third_body_sun = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_moon"])
        optFilter.third_body_moon = parameter.as<bool>();
    if (parameter = propFilterSettings["third_body_planet"])
        optFilter.third_body_planet = parameter.as<bool>();
    if (parameter = propFilterSettings["relativity_effect"])
        optFilter.relativity_effect = parameter.as<bool>();
    if (parameter = propFilterSettings["atmospheric_drag"])
        optFilter.atmospheric_drag = parameter.as<bool>();
    if (parameter = propFilterSettings["solar_radiation_pressure"])
        optFilter.solar_radiation_pressure = parameter.as<bool>();
    if (parameter = propFilterSettings["thermal_emission"])
        optFilter.thermal_emission = parameter.as<bool>();
    if (parameter = propFilterSettings["earth_albedo"])
        optFilter.earth_albedo = parameter.as<bool>();
    if (parameter = propFilterSettings["infrared_radiation"])
        optFilter.infrared_radiation = parameter.as<bool>();
    if (parameter = propFilterSettings["antenna_thrust"])
        optFilter.antenna_thrust = parameter.as<bool>();
    if (parameter = propFilterSettings["empirical_acceleration"])
        optFilter.empirical_acceleration = parameter.as<bool>();
    if (parameter = propFilterSettings["satellite_manoeuvre"])
        optFilter.satellite_manoeuvre = parameter.as<bool>();
    if (parameter = propFilterSettings["satMass"])
        optFilter.satMass = parameter.as<double>();
    if (parameter = propFilterSettings["srpArea"])
        optFilter.srpArea = parameter.as<double>();
    if (parameter = propFilterSettings["srpCoef"])
        optFilter.srpCoef = parameter.as<double>();
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

void initGlobalVariables(VectorXd &initialStateVec, string stateType)
{
    // START MOD
    initEGMCoef("./auxdata/GGM03S.txt");
    erpt = {.n = 14};

    readerp("./auxdata/cod15657.erp", &erpt);
    // cout << "erpt mjd\t" << erpt.data->mjd << endl;

    // // transform ground station from
    // leapSec = 15;

    // double erpv[4] = {};
    // geterp_from_utc(&erpt, leapSec, epoch.startMJD, erpv);
    // // cout << "erpv[3]\t" << erpv[3] << endl;

    // double dUT1_UTC = erpv[2];
    // double dUTC_TAI = -(19 + leapSec);
    // double xp = erpv[0];
    // double yp = erpv[1];
    // double lod = erpv[3];

    // iersInstance.Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
    // // double mjdTT = mjdUTC + iersInstance.TT_UTC(mjdUTC) / 86400;

    // get leap seconds from the table
    leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD));
    // set up the IERS instance
    getIERS(epoch.startMJD);

    pJPLEph = jpl_init_ephemeris(JPL_EPHEMERIS_FILENAME, nullptr, nullptr);

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

    struct SimInfo simInfo;
    struct MeasModel measMdl;
    struct Filters filters;
    struct InitialState initialState;
    // read parameter/settings from config file
    readConfigFile(configFilename, forceModelsTruthOpt, forceModelsFilterOpt, simInfo, initialState, measMdl, filters);
    epoch = simInfo.epoch;
    int numTrials = filters.numTrials;
    const int dimState = initialState.dimState;
    string initialStateType = initialState.initialStateType;
    VectorXd initialStateVec = initialState.initialStateVec;
    MatrixXd initialCov = initialState.initialCovarianceMat;
    const VectorXd groundStation = measMdl.groundStation;

    // initialise
    initGlobalVariables(initialStateVec, initialStateType);

    // setup orbit propagator
    orbitProp.setPropOption(forceModelsTruthOpt);
    orbitProp.initPropagator(initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);

    DynamicModel::stf g = accelerationModel;
    DynamicModel f(g, dimState, 1E-6, 1E-6);

    Timer timer;

    // simulate ground-truth trajectory and generate non-corrupted measurement vectors
    orbitProp.setPropOption(forceModelsTruthOpt);
    orbitProp.initPropagator(initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);
    VectorXd propStateVec = initialStateVec;
    double time = 0, dt = epoch.timeStep;
    int nTotalSteps = (epoch.endMJD - epoch.startMJD) * 86400 / dt + 1;
    // linear spaced times
    VectorXd tSec;
    tSec.setLinSpaced(nTotalSteps, 0, (nTotalSteps - 1) * dt);
    MatrixXd tableTrajTruth(nTotalSteps, dimState + 1);
    tableTrajTruth.col(0) = tSec;
    tableTrajTruth.row(0).tail(dimState) = initialStateVec;
    timer.tick();
    for (int k = 0; k < nTotalSteps - 1; k++)
    {
        propStateVec = f(time, time + dt, propStateVec, Vector3d::Zero());
        time += dt;
        tableTrajTruth.row(k + 1).tail(dimState) = propStateVec;
    }
    // header for the saved file
    vector<string> headerTraj({"tSec", "x", "y", "z", "vx", "vy", "vz"});
    string trajTruthFile = simInfo.file.outDir + "/trajectory_truth.csv";
    EigenCSV::write(tableTrajTruth, headerTraj, trajTruthFile);
}
