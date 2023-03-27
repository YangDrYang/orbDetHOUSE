#include "testOrbDet.hpp"

#define R_EARTH 6371E3
#define DEG M_PI / 180
#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)

// propagator variables
erp_t erpt;
VectorXd rvPhiS(6), groundStation(6);
double leapSec;
IERS iersInstance;
ForceModels forceModelsTruthOpt = {}, forceModelsFilterOpt = {};
EGMCoef egm;
void *pJPLEph;
Propagator orbitProp;
struct EpochInfo epoch;

using namespace Eigen;
using namespace std;
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

VectorXd accelerationModel(double t, const VectorXd &X, const VectorXd &fd)
{
    VectorXd Xf(6);
    Vector3d r;
    Vector3d v;
    VectorXd rvECI(6);
    rvECI = X;

    Matrix3d mECI2ECEF = Matrix3d::Identity();
    Matrix3d mdECI2ECEF = Matrix3d::Identity();

    eci2ecef_sofa(epoch.startMJD + t / 86400, iersInstance, mECI2ECEF, mdECI2ECEF);

    Vector3d acceleration;
    orbitProp.updPropagator(epoch.startMJD + t / 86400, leapSec, &erpt);

    // calculate acceleration
    acceleration = orbitProp.calculateAcceleration(X.head(3), X.tail(3), mECI2ECEF);

    // set state vector
    Xf.head(3) = X.tail(3);
    Xf.tail(3) = acceleration;

    return Xf;
}

Vector2d simpleMeasurementModel(double t, const VectorXd &satECI)
{
    Vector2d z;
    Vector3d p;
    VectorXd rsECEF(6); // = VectorXd::Zero(6);
    VectorXd stnECI = VectorXd::Zero(6);

    // ground station
    rsECEF = groundStation;

    ecef2eciVec_sofa(epoch.startMJD + t, iersInstance, rsECEF, stnECI);
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

Vector4d measurementModel(double t, const VectorXd &satECI, const VectorXd &stnECEF)
{
    Vector4d z;
    Vector3d p, r, v;
    VectorXd stnECI = VectorXd::Zero(6);
    VectorXd stnECEF_ = stnECEF; // define a new variable that is not a const vector

    ecef2eciVec_sofa(epoch.startMJD + t, iersInstance, stnECEF_, stnECI);
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
    return z;
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
    struct Errors errorStd;
    struct Filters filters;
    struct InitialState initialState;
    // read parameter/settings from config file
    readConfigFile(configFilename, forceModelsTruthOpt, forceModelsFilterOpt, simInfo, initialState, groundStation, filters, errorStd);
    epoch = simInfo.epoch;
    int numTrials = filters.numTrials;
    // cout << "Initial covariance matrix:\n"
    //      << initialCov << std::endl;
    const int nDim = initialState.dimState;
    string initialStateType = initialState.initialStateType;
    VectorXd initialStateVec = initialState.initialStateVec;
    MatrixXd initialCov = initialState.initialCovarianceMat;

    // initialise
    initGlobalVariables(initialStateVec, initialStateType);

    // setup orbit propagator
    orbitProp.setPropOption(forceModelsTruthOpt);
    orbitProp.initPropagator(initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);

    // UKF state & measurement models
    DynamicModel::stf g = accelerationModel;
    DynamicModel f(g, 6, 1E-6, 1E-6);
    UKF::meas_model h = [](double t, const VectorXd &x)
        -> VectorXd
    {
        return measurementModel(t, x, groundStation);
    };

    // HOUSE measurement model
    HOUSE::meas_model hh = [](double t, const VectorXd &x, const VectorXd &n)
        -> VectorXd
    {
        return measurementModel(t, x, groundStation) + n;
    };

    Matrix4d R;
    R << pow(errorStd.azimuthErr * ARC_SEC, 2), 0, 0, 0,
        0, pow(errorStd.azimuthErr * ARC_SEC, 2), 0, 0,
        0, 0, pow(errorStd.rangeErr, 2), 0,
        0, 0, 0, pow(errorStd.rangeRateErr, 2);

    MatrixXd Q = MatrixXd::Zero(nDim, nDim);
    // rx, ry, rz, vx, vy, vz
    for (int i = 0; i < nDim; i++)
    {
        Q(i, i) = 1e-9;
    }

    // Normal noise generator
    mt19937_64 mt;
    normal_distribution<double> nd;
    MatrixXd runTimesMC(numTrials, 4);

    const int numThreads = 8;
    const int trialsperThread = numTrials / numThreads;

    // Initialize UKF & CUT filters
    UKF ukf(f, h, true, 0, initialStateVec, initialCov, Q, R, UKF::sig_type::JU, 1);
    UKF cut4(f, h, true, 0, initialStateVec, initialCov, Q, R, UKF::sig_type::CUT4, 1);
    UKF cut6(f, h, true, 0, initialStateVec, initialCov, Q, R, UKF::sig_type::CUT6, 1);
    UKF cut8(f, h, true, 0, initialStateVec, initialCov, Q, R, UKF::sig_type::CUT8, 1);

    Timer timer;

    // simulate ground-truth trajectory
    vector<VectorXd> trueData;
    trueData.clear();
    orbitProp.setPropOption(forceModelsTruthOpt);
    orbitProp.initPropagator(initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);
    VectorXd X = initialStateVec;
    double time = 0, dt = epoch.timeStep;
    trueData.push_back(X);
    timer.tick();
    while (time < (epoch.endMJD - epoch.startMJD) * 86400)
    {
        Vector3d fd;
        fd << 0, 0, 0;
        X = f(time, time + epoch.timeStep, X, fd);
        time += epoch.timeStep;
        trueData.push_back(X);
    }
    MatrixXd runTimeTruth(1, 1);
    runTimeTruth(0, 0) = timer.tock();

    int steps = trueData.size();
    // linear spaced times
    VectorXd t;
    t.setLinSpaced(steps, 0, (steps - 1) * dt);

    // Save true vectors
    MatrixXd tableTrue(steps, 7);
    tableTrue.col(0) = t;
    for (int k = 0; k < steps; k++)
    {
        tableTrue.row(k).tail(6) = trueData[k];
    };

    // header for the saved file
    vector<string> header({"t", "x", "y", "z", "vx", "vy", "vz"});
    string xtrufile = simInfo.file.outDir + "/trajectory_truth";
    xtrufile += ".csv";
    EigenCSV::write(tableTrue, header, xtrufile);
    // Save truth generation run times
    vector<string> truthString;
    truthString.push_back("true trajectory");

    string time_file = simInfo.file.outDir + "/run_times_truth";
    time_file += ".csv";
    EigenCSV::write(runTimeTruth, truthString, time_file);

    // Define a random number generator
    random_device rd;
    mt19937_64 gen(rd());
    normal_distribution<double> dist(0.0, 1.0);

    // HOUSE distributions for state
    HOUSE::Dist distXi(initialCov);
    // HOUSE distributions for state noise
    HOUSE::Dist distw(Q);
    // HOUSE distributions for measurement noise
    HOUSE::Dist distn(R);

    // Initialize HOUSE
    HOUSE house(f, hh, 4, 0, distXi, distw, distn, 0);
    // perform trials
    for (int j = 1; j <= numTrials; j++)
    {
        cout << "Trial " << j << endl;

        // Generate a random vector based on the mean and covariance matrix
        MatrixXd L = initialCov.llt().matrixL();
        VectorXd z(nDim);
        for (int k = 0; k < nDim; k++)
        {
            z(k) = dist(gen);
        }
        VectorXd initialState_ = initialStateVec + L * z;
        distXi.mean = initialState_;
        // cout << (L * z).transpose() << endl;
        // cout << distXi.mean.transpose() << endl;

        house.reset(0, distXi);
        ukf.reset(0, initialState_, initialCov);
        cut4.reset(0, initialState_, initialCov);
        cut6.reset(0, initialState_, initialCov);

        // Generate true vectors
        // TODO: add montecarlo simulation
        orbitProp.setPropOption(forceModelsFilterOpt);
        orbitProp.initPropagator(initialState_, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);

        // Generate Measurement vectors
        int nMeasure = R.rows();
        vector<string> headerMeasure({"t", "range", "azimuth", "elevation"});
        MatrixXd tableMeasure(steps, nMeasure + 1);
        int numMeasure = 0;
        MatrixXd Ztru(nMeasure, steps);

        // for each time step
        for (int k = 0; k < steps; k++)
        {
            // planar condition for visibility
            // generate true measurement
            Ztru.col(k) = h(t(k), trueData[k]);

            // if !NO_MEASUREMENT
            if (Ztru(0, k) != NO_MEASUREMENT)
            {
                // corrupt measurement with noise
                // Ztru(0,k) += errorStd.rangeErr * nd(mt);
                Ztru(0, k) += errorStd.azimuthErr * ARC_SEC * nd(mt);
                Ztru(1, k) += errorStd.elevationErr * ARC_SEC * nd(mt);
                Ztru(2, k) += errorStd.rangeErr * nd(mt);
                Ztru(3, k) += errorStd.rangeRateErr * nd(mt);
                // store corrupted measurement & time
                tableMeasure(numMeasure, 0) = t(k);
                tableMeasure(numMeasure, 1) = Ztru(0, k);
                tableMeasure(numMeasure, 2) = Ztru(1, k);
                tableMeasure(numMeasure, 3) = Ztru(2, k);
                tableMeasure(numMeasure, 4) = Ztru(3, k);
                // tableMeasure(numMeasure, 3) = Ztru(2,k);

                numMeasure++;
            }
        }
        // save corrupted measurment data
        string xmeasurefile = simInfo.file.outDir + "/meas_";
        xmeasurefile += to_string(j);
        xmeasurefile += ".csv";
        EigenCSV::write(tableMeasure.topRows(numMeasure), headerMeasure, xmeasurefile);

        string outputFile;
        if (filters.house)
        {
            cout << "\tHOUSE" << '\n';
            orbitProp.initPropagator(initialState_, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            house.run(t, Ztru);
            runTimesMC(j - 1, 0) = timer.tock();

            outputFile = simInfo.file.outDir + "/house_";
            outputFile += to_string(j);
            outputFile += ".csv";
            house.save(outputFile);
        }

        // UKF Filter
        if (filters.ukf)
        {
            cout << "\tUKF" << '\n';
            orbitProp.initPropagator(initialState_, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            ukf.run(t, Ztru);
            runTimesMC(j - 1, 1) = timer.tock();

            outputFile = simInfo.file.outDir + "/ukf_";
            outputFile += to_string(j);
            outputFile += ".csv";
            ukf.save(outputFile);
        }

        // CUT-4 Filter
        if (filters.cut4)
        {
            cout << "\tCUT-4" << '\n';
            orbitProp.initPropagator(initialState_, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            cut4.run(t, Ztru);
            runTimesMC(j - 1, 2) = timer.tock();

            outputFile = simInfo.file.outDir + "/cut4_";
            outputFile += to_string(j);
            outputFile += ".csv";
            cut4.save(outputFile);
        }

        // Cut-6 Filter
        if (filters.cut6)
        {
            cout << "\tCUT-6" << '\n';
            orbitProp.initPropagator(initialState_, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            cut6.run(t, Ztru);
            runTimesMC(j - 1, 3) = timer.tock();

            outputFile = simInfo.file.outDir + "/cut6_";
            outputFile += to_string(j);
            outputFile += ".csv";
            cut6.save(outputFile);
        }
    }

    // Save Filter run times
    vector<string> filterStrings;
    filterStrings.push_back("true");
    filterStrings.push_back("house");
    filterStrings.push_back("ukf");
    filterStrings.push_back("cut4");
    filterStrings.push_back("cut6");

    time_file = simInfo.file.outDir + "/run_times_.csv";
    EigenCSV::write(runTimesMC, filterStrings, time_file);
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
                    VectorXd &groundStation, struct Filters &filters, struct Errors &errorStd)
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
    MatrixXd tempMat(dimState, dimState);
    const YAML::Node &covariance = orbitParams["initial_covariance"];
    for (int i = 0; i < dimState; ++i)
    {
        const YAML::Node &row = covariance[i];
        for (int j = 0; j < dimState; ++j)
        {
            tempMat(i, j) = row[j].as<double>();
        }
    }
    initialState.initialCovarianceMat = tempMat;

    // read measurement characristics (noise standard deviations)
    YAML::Node measParams = config["measurement_parameters"];
    tempVec = measParams["ground_station"].as<vector<double>>();
    groundStation = stdVec2EigenVec(tempVec);
    errorStd.azimuthErr = measParams["elevation_error"].as<double>();
    errorStd.elevationErr = measParams["azimuth_error"].as<double>();
    errorStd.rangeErr = measParams["range_error"].as<double>();
    errorStd.rangeRateErr = measParams["range_rate_error"].as<double>();

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

    // transform ground station from
    leapSec = 0;

    double erpv[4] = {};
    geterp_from_utc(&erpt, leapSec, epoch.startMJD, erpv);

    double dUT1_UTC = erpv[2];
    double dUTC_TAI = -(19 + leapSec);
    double xp = erpv[0];
    double yp = erpv[1];
    double lod = erpv[3];

    iersInstance.Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
    // double mjdTT = mjdUTC + iersInstance.TT_UTC(mjdUTC) / 86400;

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

    // END MOD
    int nState = 6;
    int nVar = 7; // 6 state variables and one orbital parameter Cr
    rvPhiS = VectorXd::Zero(nState + nState * nVar);
    // Create combined vector from epoch state, epoch transition matrix (=1) and epoch sensitivity matrix (=0)
    for (int i = 0; i < nState; i++)
    {
        rvPhiS(i) = rvECI(i);
        for (int j = 0; j < nVar; j++)
        {
            rvPhiS(nState * (j + 1) + i) = (i == j ? 1 : 0);
        }
    }
}