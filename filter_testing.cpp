#include "filter_testing.hpp"

#define R_EARTH 6.3781E6 // in metres
#define DEG M_PI / 180
#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)

// propagator variables
erp_t erpt;
VectorXd rvPhiS(6), initialState(6), groundStation(6);
double leapSec;
IERS iersInstance;
ForceModels forceModelsOpt = {};
EGMCoef egm;
void *pJPLEph;
Propagator orbitProp;
struct EpochInfo epoch;

using namespace Eigen;
using namespace std;

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
    Xf.tail(3) = acceleration; // + dragAcceleration;

    return Xf;
}

Vector2d simpleMeasurementModel(double t, const VectorXd &X)
{
    Vector2d z;
    Vector3d p;
    VectorXd rsECEF(6); // = VectorXd::Zero(6);
    VectorXd rsECI = VectorXd::Zero(6);

    // ground station
    rsECEF = groundStation;

    ecef2eciVec_sofa(epoch.startMJD + t, iersInstance, rsECEF, rsECI);
    // end transformation code
    p = X.head(3) - rsECI.head(3);

    if (p.dot(rsECI.head(3)) >= 0)
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

Vector4d measurementModel(double t, const VectorXd &X)
{
    Vector4d z;
    Vector3d p, r, v;
    VectorXd rsECEF(6); // = VectorXd::Zero(6);
    VectorXd rsECI = VectorXd::Zero(6);

    // ground station
    rsECEF = groundStation;

    ecef2eciVec_sofa(epoch.startMJD + t, iersInstance, rsECEF, rsECI);
    // end transformation code
    p = X.head(3) - rsECI.head(3);
    v = X.tail(3) - rsECI.tail(3);
    if (p.dot(rsECI.head(3)) >= 0)
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

    struct Errors errorStd;
    struct Filters filters;
    int numTrials;
    string initialStateType;
    // read parameter/settings from config file
    readConfigFile(configFilename, forceModelsOpt, epoch, initialState, groundStation, filters, numTrials, errorStd, initialStateType);

    // initialise
    initGlobalVariables(initialState, initialStateType);

    // setup orbit propagator
    orbitProp.setPropOption(forceModelsOpt);
    orbitProp.initPropagator(initialState, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);

    // UKF state & measurement models
    DynamicModel::stf g = accelerationModel;
    DynamicModel f(g, 6, 1E-6, 1E-6);
    // UKF::dyn_model f = accelerationModel;
    UKF::meas_model h = measurementModel;

    // HOUSE measurement model
    HOUSE::meas_model hh = [](double t, const VectorXd &x, const VectorXd &n)
        -> VectorXd
    {
        return measurementModel(t, x) + n;
    };

    // Constants
    // double sr, sth, l1, l2;
    // sr  = 100;
    // sth = 1 * deg;
    // l1  = 0.16;
    // l2  = 0.01;

    // Measurement noise covariance
    // Matrix3d R;
    // R << pow(errorStd.rangeErr * ARC_SEC, 2), 0, 0,
    //      0, pow(errorStd.azimuthErr  * ARC_SEC, 2), 0,
    //      0, 0, pow(errorStd.elevationErr, 2);

    Matrix4d R;
    R << pow(errorStd.azimuthErr * ARC_SEC, 2), 0, 0, 0,
        0, pow(errorStd.azimuthErr * ARC_SEC, 2), 0, 0,
        0, 0, pow(errorStd.rangeErr, 2), 0,
        0, 0, 0, pow(errorStd.rangeRateErr, 2);

    // Prior mean & covariance
    VectorXd mxi(6), cxx(6);
    cxx << 1, 1, 1, 1e-5, 1e-5, 1e-5; // prior state standard deviation
    mxi = initialState;               // prior state mean
    MatrixXd Pxxi(6, 6);
    Pxxi << 1.481e2, 0, 0, 0, -9.237e-2, -5.333e-2,
        0, 2.885e1, 9.994, -3.121e-2, 0, 0,
        0, 9.994, 5.770, -1.242, 0, 0,
        0, -3.121e-2, -1.242e-2, 3.687e-5, 0, 0,
        -9.237e-3, 0, 0, 0, 6.798e-5, 3.145e-5,
        -5333e-3, 0, 0, 0, 3.145e-5, 3.166e-5;
    Pxxi *= 1e6;

    // Pxxi = MatrixXd::Identity(6,6) * cxx;
    // HOUSE distributions
    HOUSE::Dist distXi(Pxxi), distn(R);
    distXi.mean = mxi;

    // Normal noise generator
    mt19937_64 mt;
    normal_distribution<double> nd;
    vector<string> header({"t", "x", "y", "z", "vx", "vy", "vz"});
    MatrixXd run_times(numTrials, 5);

    MatrixXd Q = MatrixXd::Zero(6, 6);
    // rx, ry, rz, vx, vy, vz
    for (int i = 0; i < 6; i++)
    {
        Q(i, i) = 1e-9;
    }

    HOUSE::Dist distw(Q);

    // Initialize UKF & CUT filters
    UKF ukf(f, h, true, 0, mxi, Pxxi, Q, R, UKF::sig_type::JU, 1);
    UKF cut4(f, h, true, 0, mxi, Pxxi, Q, R, UKF::sig_type::CUT4, 1);
    UKF cut6(f, h, true, 0, mxi, Pxxi, Q, R, UKF::sig_type::CUT6, 1);
    UKF cut8(f, h, true, 0, mxi, Pxxi, Q, R, UKF::sig_type::CUT8, 1);

    // Initialize HOUSE
    HOUSE house(f, hh, 4, 0, distXi, distw, distn, 0);

    Timer timer;
    // perform trials
    for (int j = 1; j <= numTrials; j++)
    {
        cout << "Trial " << j << endl;
        house.reset(0, distXi);
        ukf.reset(0, initialState, Pxxi);
        cut4.reset(0, initialState, Pxxi);
        cut6.reset(0, initialState, Pxxi);

        // Generate true vectors
        // TODO: add montecarlo simulation
        vector<VectorXd> trueData;
        trueData.clear();
        orbitProp.setPropOption(forceModelsOpt);
        orbitProp.initPropagator(initialState, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);

        VectorXd X = initialState;
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
        run_times(j - 1, 0) = timer.tock();

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
        string xtrufile = "out/tru_";
        xtrufile += to_string(j);
        xtrufile += ".csv";
        EigenCSV::write(tableTrue, header, xtrufile);

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
        string xmeasurefile = "out/meas_";
        xmeasurefile += to_string(j);
        xmeasurefile += ".csv";
        EigenCSV::write(tableMeasure.topRows(numMeasure), headerMeasure, xmeasurefile);

        string outputFile;
        if (filters.house)
        {
            cout << "\tHOUSE" << '\n';
            orbitProp.initPropagator(initialState, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            house.run(t, Ztru);
            run_times(j - 1, 1) = timer.tock();

            outputFile = "out";
            outputFile += "/house_";
            outputFile += to_string(j);
            outputFile += ".csv";
            house.save(outputFile);
        }

        // UKF Filter
        if (filters.ukf)
        {
            cout << "\tUKF" << '\n';
            orbitProp.initPropagator(initialState, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            ukf.run(t, Ztru);
            run_times(j - 1, 2) = timer.tock();

            outputFile = "out";
            outputFile += "/ukf_";
            outputFile += to_string(j);
            outputFile += ".csv";
            ukf.save(outputFile);
        }

        // CUT-4 Filter
        if (filters.cut4)
        {
            cout << "\tCUT-4" << '\n';
            orbitProp.initPropagator(initialState, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            cut4.run(t, Ztru);
            run_times(j - 1, 3) = timer.tock();

            outputFile = "out";
            outputFile += "/cut4_";
            outputFile += to_string(j);
            outputFile += ".csv";
            cut4.save(outputFile);
        }

        // Cut-6 Filter
        if (filters.cut6)
        {
            cout << "\tCUT-6" << '\n';
            orbitProp.initPropagator(initialState, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            cut6.run(t, Ztru);
            run_times(j - 1, 4) = timer.tock();

            outputFile = "out";
            outputFile += "/cut6_";
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

    string time_file = "out/run_times_";
    time_file += ".csv";

    EigenCSV::write(run_times, filterStrings, time_file);
}

Eigen::MatrixXd generateTrueResults(DynamicModel &f, struct EpochInfo epoch, VectorXd initialState)
{
    // create results matrix
    int length = floor((epoch.endMJD - epoch.startMJD) / (epoch.timeStep / 86400.0)) + 1;
    MatrixXd results(length, 7);

    double time = 0;
    double step = epoch.timeStep;

    Vector3d W = Vector3d::Zero();
    VectorXd X = initialState;
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

void readConfigFile(string fileName, ForceModels &options, struct EpochInfo &epoch, Eigen::VectorXd &initialState,
                    Eigen::VectorXd &groundStation, struct Filters &filters, int &numTrials, struct Errors &errorStd,
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
    YAML::Node orbitParams = config["orbtial_parameters"];
    std::vector<double> tempVec;
    double test = orbitParams["MJD_start"].as<double>();
    epoch.startMJD = orbitParams["MJD_start"].as<double>();
    epoch.endMJD = orbitParams["MJD_end"].as<double>();
    epoch.timeStep = orbitParams["time_step"].as<double>();
    errorStd.azimuthErr = orbitParams["elevation_error"].as<double>();
    errorStd.elevationErr = orbitParams["azimuth_error"].as<double>();
    errorStd.rangeErr = orbitParams["range_error"].as<double>();
    errorStd.rangeRateErr = orbitParams["range_rate_error"].as<double>();
    iniatialStateType = orbitParams["initial_state_type"].as<string>();

    // read params as standard vector, convert to eigen vector
    tempVec = orbitParams["initial_state"].as<std::vector<double>>();
    initialState = stdVec2EigenVec(tempVec);
    tempVec = orbitParams["ground_station"].as<std::vector<double>>();
    groundStation = stdVec2EigenVec(tempVec);

    // read propagator settings (optional)
    YAML::Node propSettings = config["propagator_settings"];
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

Eigen::VectorXd stdVec2EigenVec(const std::vector<double> &stdVec)
{
    Eigen::VectorXd eigenVec(stdVec.size());
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

void initGlobalVariables(VectorXd &initialState, string stateType)
{
    // START MOD
    initEGMCoef("GGM03S.txt");
    erpt = {.n = 14};

    readerp("cod15657.erp", &erpt);

    // transform ground station from
    leapSec = 0;

    double erpv[0] = {};
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
        ecef2eciVec_sofa(epoch.startMJD, iersInstance, initialState, rvECI);
        initialState = rvECI;
    }
    else
    {
        rvECI = initialState;
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
