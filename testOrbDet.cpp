#include "testOrbDet.hpp"

#define R_EARTH 6371E3
#define DEG M_PI / 180
#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)

// propagator variables
erp_t erpt;
double leapSec;
IERS iersInstance;
ForceModels forceModelsTruthOpt = {}, forceModelsFilterOpt = {};
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

// Vector2d simpleMeasurementModel(double t, const VectorXd &satECI)
// {
//     Vector2d z;
//     Vector3d p;
//     VectorXd rsECEF(6); // = VectorXd::Zero(6);
//     VectorXd stnECI = VectorXd::Zero(6);

//     // ground station
//     rsECEF = groundStation;

//     ecef2eciVec_sofa(epoch.startMJD + t, iersInstance, rsECEF, stnECI);
//     // end transformation code
//     p = satECI.head(3) - stnECI.head(3);

//     if (p.dot(stnECI.head(3)) >= 0)
//     {
//         // range
//         // azimuth angle
//         z(0) = atan2(p(1), p(0));
//         // elevation angle
//         z(1) = asin(p(2) / p.norm());
//     }
//     // not visible
//     else
//     {
//         z(0) = NO_MEASUREMENT;
//         z(1) = NO_MEASUREMENT;
//     }
//     return z;
// }

Vector4d measurementModel(double t, const VectorXd &satECI, const VectorXd &stnECEF)
{
    Vector4d z;
    Vector3d p, r, v;
    VectorXd stnECI = VectorXd::Zero(6);
    VectorXd stnECEF_ = stnECEF; // define a new variable that is not a const vector

    ecef2eciVec_sofa(epoch.startMJD + t / 86400, iersInstance, stnECEF_, stnECI);
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
    struct MeasModel measMdl;
    struct Filters filters;
    struct InitialState initialState;
    // read parameter/settings from config file
    readConfigFile(configFilename, forceModelsTruthOpt, forceModelsFilterOpt, simInfo, initialState, measMdl, filters);
    epoch = simInfo.epoch;
    int numTrials = filters.numTrials;
    const VectorXd groundStation = measMdl.groundStation;

    // initialise
    initGlobalVariables(initialState.initialStateVec, initialState.initialStateType);

    // setup orbit propagator
    orbitProp.setPropOption(forceModelsTruthOpt);
    orbitProp.initPropagator(initialState.initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);

    // UKF state & measurement models
    DynamicModel::stf g = accelerationModel;
    DynamicModel f(g, 6, 1E-6, 1E-6);
    UKF::meas_model h = [&groundStation](double t, const VectorXd &x)
        -> VectorXd
    {
        return measurementModel(t, x, groundStation);
    };

    // HOUSE measurement model
    HOUSE::meas_model hh = [&groundStation](double t, const VectorXd &x, const VectorXd &n)
        -> VectorXd
    {
        return measurementModel(t, x, groundStation) + n;
    };

    Matrix4d measNoiseCov;
    measNoiseCov << pow(measMdl.errorStd.azimuthErr * ARC_SEC, 2), 0, 0, 0,
        0, pow(measMdl.errorStd.azimuthErr * ARC_SEC, 2), 0, 0,
        0, 0, pow(measMdl.errorStd.rangeErr, 2), 0,
        0, 0, 0, pow(measMdl.errorStd.rangeRateErr, 2);

    // Prior mean & covariance
    VectorXd mxi = initialState.initialStateVec; // prior state mean
    MatrixXd Pxxi(6, 6);
    Pxxi << 1.481e2, 0, 0, 0, -9.237e-2, -5.333e-2,
        0, 2.885e1, 9.994, -3.121e-2, 0, 0,
        0, 9.994, 5.770, -1.242e-2, 0, 0,
        0, -3.121e-2, -1.242e-2, 3.687e-5, 0, 0,
        -9.237e-2, 0, 0, 0, 6.798e-5, 3.145e-5,
        -5.333e-2, 0, 0, 0, 3.145e-5, 3.166e-5;
    Pxxi *= 1e6;

    // HOUSE distributions for state
    HOUSE::Dist distXi(Pxxi);
    distXi.mean = mxi;

    MatrixXd procNoiseCov = MatrixXd::Zero(6, 6);
    // rx, ry, rz, vx, vy, vz
    for (int i = 0; i < 6; i++)
    {
        procNoiseCov(i, i) = 1e-9;
    }
    // HOUSE distributions for state noise
    HOUSE::Dist distw(procNoiseCov);

    // HOUSE distributions for measurement noise
    HOUSE::Dist distn(measNoiseCov);

    // run time
    MatrixXd run_times(numTrials, 4);

    Timer timer;
    // simulate ground-truth trajectory
    vector<VectorXd> trueData;
    trueData.clear();
    orbitProp.setPropOption(forceModelsTruthOpt);
    orbitProp.initPropagator(initialState.initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);

    VectorXd X = initialState.initialStateVec;
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

    int nTotalSteps = trueData.size();
    // linear spaced times
    VectorXd t;
    t.setLinSpaced(nTotalSteps, 0, (nTotalSteps - 1) * dt);

    MatrixXd tableTrue(nTotalSteps, 7);
    tableTrue.col(0) = t;
    // generate non-corrupted measurement vectors
    int dimMeas = measMdl.dimMeas;
    MatrixXd measTruth(dimMeas, nTotalSteps);
    // for each time step
    for (int k = 0; k < nTotalSteps; k++)
    {
        // Save true vectors
        tableTrue.row(k).tail(6) = trueData[k];
        // generate true measurement
        measTruth.col(k) = h(t(k), trueData[k]);
    };

    // header for the saved file
    vector<string> headerTraj({"t", "x", "y", "z", "vx", "vy", "vz"});
    string trajTruthFile = simInfo.file.outDir + "/trajectory_truth";
    trajTruthFile += ".csv";
    EigenCSV::write(tableTrue, headerTraj, trajTruthFile);

    // Initialize UKF & CUT filters
    UKF ukf(f, h, true, 0, mxi, Pxxi, procNoiseCov, measNoiseCov, UKF::sig_type::JU, 1);
    UKF cut4(f, h, true, 0, mxi, Pxxi, procNoiseCov, measNoiseCov, UKF::sig_type::CUT4, 1);
    UKF cut6(f, h, true, 0, mxi, Pxxi, procNoiseCov, measNoiseCov, UKF::sig_type::CUT6, 1);
    UKF cut8(f, h, true, 0, mxi, Pxxi, procNoiseCov, measNoiseCov, UKF::sig_type::CUT8, 1);
    // Initialize HOUSE
    HOUSE house(f, hh, 4, 0, distXi, distw, distn, 0);

    vector<string> headerMeas({"t", "azimuth", "elevation", "range", "range rate"});

    // Normal noise generator
    mt19937_64 gen;
    normal_distribution<double> dist;

    // perform trials
    for (int j = 1; j <= numTrials; j++)
    {
        cout << "Trial " << j << endl;
        house.reset(0, distXi);
        // cout << "initialStateVec\n"
        //      << initialState.initialStateVec << endl;
        // cout << "initialCov\n"
        //      << Pxxi << endl;
        ukf.reset(0, initialState.initialStateVec, Pxxi);
        cut4.reset(0, initialState.initialStateVec, Pxxi);
        cut6.reset(0, initialState.initialStateVec, Pxxi);

        MatrixXd tableMeasure(nTotalSteps, dimMeas + 1);
        int numMeasure = 0;

        // for each time step
        for (int k = 0; k < nTotalSteps; k++)
        {
            // if !NO_MEASUREMENT
            if (measTruth(0, k) != NO_MEASUREMENT)
            {
                // corrupt measurement with noise
                // measTruth(0,k) += errorStd.rangeErr * nd(mt);
                measTruth(0, k) += measMdl.errorStd.azimuthErr * ARC_SEC * dist(gen);
                measTruth(1, k) += measMdl.errorStd.elevationErr * ARC_SEC * dist(gen);
                measTruth(2, k) += measMdl.errorStd.rangeErr * dist(gen);
                measTruth(3, k) += measMdl.errorStd.rangeRateErr * dist(gen);
                // store corrupted measurement & time
                tableMeasure(numMeasure, 0) = t(k);
                tableMeasure(numMeasure, 1) = measTruth(0, k);
                tableMeasure(numMeasure, 2) = measTruth(1, k);
                tableMeasure(numMeasure, 3) = measTruth(2, k);
                tableMeasure(numMeasure, 4) = measTruth(3, k);

                numMeasure++;
            }
        }
        // save corrupted measurment data
        string measFile = "out/meas_";
        measFile += to_string(j);
        measFile += ".csv";
        EigenCSV::write(tableMeasure.topRows(numMeasure), headerMeas, measFile);

        string outputFile;
        if (filters.house)
        {
            cout << "\tHOUSE" << '\n';
            orbitProp.initPropagator(initialState.initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            house.run(t, measTruth);
            run_times(j - 1, 0) = timer.tock();

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
            // cout << "startMJD\t" << epoch.startMJD << endl;
            // cout << "leapSec\t" << leapSec << endl;
            orbitProp.initPropagator(initialState.initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            ukf.run(t, measTruth);
            run_times(j - 1, 1) = timer.tock();

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
            orbitProp.initPropagator(initialState.initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            cut4.run(t, measTruth);
            run_times(j - 1, 2) = timer.tock();

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
            orbitProp.initPropagator(initialState.initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            timer.tick();
            cut6.run(t, measTruth);
            run_times(j - 1, 3) = timer.tock();

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

// void readConfigFile(string fileName, ForceModels &optTruth, ForceModels &optFilter, struct EpochInfo &epoch, Eigen::VectorXd &initialState,
//                     Eigen::VectorXd &groundStation, struct Filters &filters, int &numTrials, struct Errors &errorStd,
//                     string &iniatialStateType)
// {
//     // load file
//     YAML::Node config = YAML::LoadFile(fileName);
//     YAML::Node parameter;

//     // read filter options (required)
//     YAML::Node filterOpts = config["filter_options"];
//     filters.house = filterOpts["HOUSE"].as<bool>();
//     filters.ukf = filterOpts["UKF"].as<bool>();
//     filters.cut4 = filterOpts["CUT4"].as<bool>();
//     filters.cut6 = filterOpts["CUT6"].as<bool>();
//     numTrials = filterOpts["num_trials"].as<int>();

//     // read simulation parameters (required)
//     YAML::Node simParams = config["simulation_parameters"];
//     epoch.startMJD = simParams["MJD_start"].as<double>();
//     epoch.endMJD = simParams["MJD_end"].as<double>();
//     epoch.timeStep = simParams["time_step"].as<double>();

//     // read orbital parameters (required)
//     YAML::Node orbitParams = config["initial_orbtial_parameters"];
//     std::vector<double> tempVec;
//     iniatialStateType = orbitParams["initial_state_type"].as<string>();
//     // read params as standard vector, convert to eigen vector
//     tempVec = orbitParams["initial_state"].as<std::vector<double>>();
//     initialState = stdVec2EigenVec(tempVec);

//     // read measurement characristics (noise standard deviations)
//     YAML::Node measParams = config["measurement_parameters"];
//     tempVec = measParams["ground_station"].as<std::vector<double>>();
//     groundStation = stdVec2EigenVec(tempVec);
//     errorStd.azimuthErr = measParams["elevation_error"].as<double>();
//     errorStd.elevationErr = measParams["azimuth_error"].as<double>();
//     errorStd.rangeErr = measParams["range_error"].as<double>();
//     errorStd.rangeRateErr = measParams["range_rate_error"].as<double>();

//     // read propagator settings (optional)
//     YAML::Node propTruthSettings = config["propagator_truth_settings"];
//     if (parameter = propTruthSettings["earth_gravaity"])
//         optTruth.earth_gravity = parameter.as<bool>();
//     if (parameter = propTruthSettings["solid_earth_tide"])
//         optTruth.solid_earth_tide = parameter.as<bool>();
//     if (parameter = propTruthSettings["ocean_tide_loading"])
//         optTruth.ocean_tide_loading = parameter.as<bool>();
//     if (parameter = propTruthSettings["third_body_attraction"])
//         optTruth.third_body_attraction = parameter.as<bool>();
//     if (parameter = propTruthSettings["third_body_sun"])
//         optTruth.third_body_sun = parameter.as<bool>();
//     if (parameter = propTruthSettings["third_body_moon"])
//         optTruth.third_body_moon = parameter.as<bool>();
//     if (parameter = propTruthSettings["third_body_planet"])
//         optTruth.third_body_planet = parameter.as<bool>();
//     if (parameter = propTruthSettings["relativity_effect"])
//         optTruth.relativity_effect = parameter.as<bool>();
//     if (parameter = propTruthSettings["atmospheric_drag"])
//         optTruth.atmospheric_drag = parameter.as<bool>();
//     if (parameter = propTruthSettings["solar_radiation_pressure"])
//         optTruth.solar_radiation_pressure = parameter.as<bool>();
//     if (parameter = propTruthSettings["thermal_emission"])
//         optTruth.thermal_emission = parameter.as<bool>();
//     if (parameter = propTruthSettings["earth_albedo"])
//         optTruth.earth_albedo = parameter.as<bool>();
//     if (parameter = propTruthSettings["infrared_radiation"])
//         optTruth.infrared_radiation = parameter.as<bool>();
//     if (parameter = propTruthSettings["antenna_thrust"])
//         optTruth.antenna_thrust = parameter.as<bool>();
//     if (parameter = propTruthSettings["empirical_acceleration"])
//         optTruth.empirical_acceleration = parameter.as<bool>();
//     if (parameter = propTruthSettings["satellite_manoeuvre"])
//         optTruth.satellite_manoeuvre = parameter.as<bool>();
//     if (parameter = propTruthSettings["satMass"])
//         optTruth.satMass = parameter.as<double>();
//     if (parameter = propTruthSettings["srpArea"])
//         optTruth.srpArea = parameter.as<double>();
//     if (parameter = propTruthSettings["srpCoef"])
//         optTruth.srpCoef = parameter.as<double>();

//     // read propagator settings for filters (optional)
//     YAML::Node propFilterSettings = config["propagator_filter_settings"];
//     if (parameter = propFilterSettings["earth_gravaity"])
//         optFilter.earth_gravity = parameter.as<bool>();
//     if (parameter = propFilterSettings["solid_earth_tide"])
//         optFilter.solid_earth_tide = parameter.as<bool>();
//     if (parameter = propFilterSettings["ocean_tide_loading"])
//         optFilter.ocean_tide_loading = parameter.as<bool>();
//     if (parameter = propFilterSettings["third_body_attraction"])
//         optFilter.third_body_attraction = parameter.as<bool>();
//     if (parameter = propFilterSettings["third_body_sun"])
//         optFilter.third_body_sun = parameter.as<bool>();
//     if (parameter = propFilterSettings["third_body_moon"])
//         optFilter.third_body_moon = parameter.as<bool>();
//     if (parameter = propFilterSettings["third_body_planet"])
//         optFilter.third_body_planet = parameter.as<bool>();
//     if (parameter = propFilterSettings["relativity_effect"])
//         optFilter.relativity_effect = parameter.as<bool>();
//     if (parameter = propFilterSettings["atmospheric_drag"])
//         optFilter.atmospheric_drag = parameter.as<bool>();
//     if (parameter = propFilterSettings["solar_radiation_pressure"])
//         optFilter.solar_radiation_pressure = parameter.as<bool>();
//     if (parameter = propFilterSettings["thermal_emission"])
//         optFilter.thermal_emission = parameter.as<bool>();
//     if (parameter = propFilterSettings["earth_albedo"])
//         optFilter.earth_albedo = parameter.as<bool>();
//     if (parameter = propFilterSettings["infrared_radiation"])
//         optFilter.infrared_radiation = parameter.as<bool>();
//     if (parameter = propFilterSettings["antenna_thrust"])
//         optFilter.antenna_thrust = parameter.as<bool>();
//     if (parameter = propFilterSettings["empirical_acceleration"])
//         optFilter.empirical_acceleration = parameter.as<bool>();
//     if (parameter = propFilterSettings["satellite_manoeuvre"])
//         optFilter.satellite_manoeuvre = parameter.as<bool>();
//     if (parameter = propFilterSettings["satMass"])
//         optFilter.satMass = parameter.as<double>();
//     if (parameter = propFilterSettings["srpArea"])
//         optFilter.srpArea = parameter.as<double>();
//     if (parameter = propFilterSettings["srpCoef"])
//         optFilter.srpCoef = parameter.as<double>();
// }
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
    initEGMCoef("./auxdata/GGM03S.txt");
    erpt = {.n = 14};

    readerp("./auxdata/cod15657.erp", &erpt);
    // cout << "erpt mjd\t" << erpt.data->mjd << endl;

    // transform ground station from
    leapSec = 15;

    double erpv[4] = {};
    geterp_from_utc(&erpt, leapSec, epoch.startMJD, erpv);
    // cout << "erpv[3]\t" << erpv[3] << endl;

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
}