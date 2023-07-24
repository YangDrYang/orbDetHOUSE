#include "testOrbDet.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>

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
// double leapSec;
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
    // get leap seconds from the table
    double leapSec = -getLeapSecond(convertMJD2Time_T(mjd));

    double erpv[4] = {};
    // cout << "leap second:   " << leapSec << "   "
    //      << "mjd:   " << mjd << endl;
    geterp_from_utc(&erpt, leapSec, mjd, erpv);

    double dUT1_UTC = erpv[2];
    double dUTC_TAI = -(19 + leapSec);
    double xp = erpv[0];
    double yp = erpv[1];
    double lod = erpv[3];

    iersInstance.Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);

    // cout << "dUT1_UTC:  " << dUT1_UTC << "dUTC_TAI: " << dUTC_TAI << "xp:   " << xp << endl;
}

VectorXd accelerationModel(double tSec, const VectorXd &X, const VectorXd &fd)
{
    VectorXd Xf(6);

    Matrix3d mECI2ECEF = Matrix3d::Identity();
    Matrix3d mdECI2ECEF = Matrix3d::Identity();

    // get leap seconds from the table
    double leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD + tSec / 86400));
    // set up the IERS instance
    getIERS(epoch.startMJD + tSec / 86400);

    eci2ecef_sofa(epoch.startMJD + tSec / 86400, iersInstance, mECI2ECEF, mdECI2ECEF);

    Vector3d acceleration;
    orbitProp.updPropagator(epoch.startMJD + tSec / 86400, leapSec, &erpt);

    // calculate acceleration
    if (abs(X(1)) < 1 && abs(X(2)) < 1) // modified equnoctial elements
    {
        Xf = orbitProp.calculateTimeDerivativeMEE(X, mECI2ECEF);

        // cout << "xf in orbit prediction: \t" << Xf << endl;
    }
    else // Cartesian elements
    {

        acceleration = orbitProp.calculateAcceleration(X.head(3), X.tail(3), mECI2ECEF);
        // set state vector
        Xf.head(3) = X.tail(3);
        Xf.tail(3) = acceleration;
    }

    return Xf;
}

Vector2d measurementModel(double tSec, const VectorXd &satECI, const VectorXd &stnECEF)
{
    Vector2d z;
    Vector3d p, r, v;
    VectorXd stnECI = VectorXd::Zero(6);
    VectorXd stnECEF_ = stnECEF; // define a new variable that is not a const vector

    // // get leap seconds from the table
    // double leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD + tSec / 86400));

    // cout << "tSec:  " << tSec << endl;
    // cout << "leapSec:   " << leapSec << endl;

    // set up the IERS instance
    getIERS(epoch.startMJD + tSec / 86400);

    // cout << "starting mjd:  " << epoch.startMJD << endl;
    ecef2eciVec_sofa(epoch.startMJD + tSec / 86400, iersInstance, stnECEF_, stnECI);
    // cout << "epoch" << tSec << endl;
    // cout << "ground station in ECEF " << stnECEF_.transpose() << endl;
    // cout << "ground station in ECI  " << stnECI.transpose() << endl;

    // end transformation code
    p = satECI.head(3) - stnECI.head(3);
    v = satECI.tail(3) - stnECI.tail(3);

    // right ascension angle
    z(0) = atan2(p(1), p(0));
    if (z(0) < 0)
    {
        z(0) += 2 * M_PI;
    };

    // declination angle
    z(1) = asin(p(2) / p.norm());

    // cout << "calculated measurements: " << z(0) << "\t" << z(1) << endl;

    return z;
}

Vector2d measurementMEEModel(double tSec, const VectorXd &satMEE, const VectorXd &stnECEF)
{
    Vector2d z;
    Vector3d p, r, v;
    VectorXd satECI = VectorXd::Zero(6), stnECI = VectorXd::Zero(6);
    VectorXd stnECEF_ = stnECEF; // define a new variable that is not a const vector

    // set up the IERS instance
    getIERS(epoch.startMJD + tSec / 86400);

    // transfer the station coordinate from ECEF to ECI
    ecef2eciVec_sofa(epoch.startMJD + tSec / 86400, iersInstance, stnECEF_, stnECI);

    // convert the modified equinoctial elements to ECI coordiantes
    satECI = coe2eci(mee2coe(satMEE), GM_Earth);
    // cout << "satMEE in measurement\t" << satMEE << endl;
    // cout << "satECI in measurement\t" << satECI << endl;

    // end transformation code
    p = satECI.head(3) - stnECI.head(3);
    v = satECI.tail(3) - stnECI.tail(3);

    // right ascension angle
    z(0) = atan2(p(1), p(0));
    if (z(0) < 0)
    {
        z(0) += 2 * M_PI;
    };

    // declination angle
    z(1) = asin(p(2) / p.norm());

    // cout << "calculated measurements: " << z(0) << "\t" << z(1) << endl;

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

MatrixXd readCSV(const string &filename, int headerLinesToSkip)
{
    ifstream file(filename);
    vector<vector<double>> data;

    string line;
    int linesSkipped = 0;

    while (getline(file, line))
    {
        if (linesSkipped < headerLinesToSkip)
        {
            linesSkipped++;
            continue; // Skip header lines
        }

        istringstream iss(line);
        string value;
        vector<double> row;

        while (getline(iss, value, ','))
        {
            try
            {
                row.push_back(stod(value));
            }
            catch (const invalid_argument &e)
            {
                // Handle conversion errors, e.g., non-numeric values
                cerr << "Invalid argument: " << e.what() << endl;
                // You can choose to skip this value or handle it differently
                row.push_back(0.0); // Default value if conversion fails
            }
        }

        data.push_back(row);
    }

    int rows = static_cast<int>(data.size());
    int cols = static_cast<int>(data[0].size());

    MatrixXd matrix(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            matrix(i, j) = data[i][j];
        }
    }

    return matrix;
}

void readConfigFile(string fileName, ForceModels &optFilter, struct ScenarioInfo &snrInfo, struct InitialState &initialState,
                    struct MeasModel &measMdl, struct Filters &filters, struct FileInfo &suppFiles)
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

    // read scenario parameters (required)
    YAML::Node snrParams = config["scenario_parameters"];
    snrInfo.epoch.startMJD = snrParams["MJD_start"].as<double>();
    snrInfo.epoch.endMJD = snrParams["MJD_end"].as<double>();
    snrInfo.epoch.maxTimeStep = snrParams["max_time_step"].as<double>();
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
    initialState.initialSkewness = orbitParams["initial_skewness"].as<double>();
    initialState.initialKurtosis = orbitParams["initial_kurtosis"].as<double>();

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
    measMdl.measFile = measParams["meas_file"].as<string>();
    tempVec = measParams["ground_station"].as<vector<double>>();
    measMdl.groundStation = stdVec2EigenVec(tempVec);
    measMdl.dimMeas = measParams["dim_meas"].as<int>();
    measMdl.errorStd.rightAscensionErr = measParams["right_ascension_error"].as<double>();
    measMdl.errorStd.declinationErr = measParams["declination_error"].as<double>();

    // read propagator settings for filters (optional)
    YAML::Node propFilterSettings = config["propagator_filter_settings"];
    if (parameter = propFilterSettings["earth_gravaity"])
        optFilter.earth_gravity = parameter.as<bool>();
    if (parameter = propFilterSettings["earth_gravity_model_order"])
        optFilter.egmAccOrd = parameter.as<int>();
    if (parameter = propFilterSettings["earth_gravity_model_degree"])
        optFilter.egmAccDeg = parameter.as<int>();
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
    if (parameter = propFilterSettings["dragArea"])
        optFilter.dragArea = parameter.as<double>();
    if (parameter = propFilterSettings["dragCoef"])
        optFilter.dragCoef = parameter.as<double>();

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

// void initGlobalVariables(struct InitialState &initialState, struct FileInfo &suppFiles)
void initGlobalVariables(VectorXd &initialStateVec, MatrixXd &initialCov, MatrixXd &procNoiseCov, string stateType, struct FileInfo &suppFiles)
{
    // string stateType = initialState.initialStateType;
    // VectorXd initialStateVec = initialState.initialStateVec;
    // MatrixXd initialCov = initialState.initialCovarianceMat;
    // MatrixXd procNoiseCov = initialState.processNoiseCovarianceMat;

    initEGMCoef(suppFiles.grvFile);

    erpt = {.n = 0};
    // cout << suppFiles.erpFile << endl;
    readerp(suppFiles.erpFile, &erpt);

    // set up the IERS instance
    getIERS(epoch.startMJD);

    const char *ephFile = suppFiles.ephFile.c_str();
    pJPLEph = jpl_init_ephemeris(ephFile, nullptr, nullptr);

    double dimState = initialStateVec.size();
    if (stateType == "ecef")
    {
        cout << "Converted state from ECEF to ECI\n";
        VectorXd satECI = VectorXd::Zero(dimState);
        ecef2eciVec_sofa(epoch.startMJD, iersInstance, initialStateVec, satECI);
        initialStateVec = satECI;
    }
    else if (stateType == "mee")
    {
        cout << "Convert state from ECI to MEE\n";

        MatrixXd rvCov = initialCov;

        const double mu = GM_Earth;
        const VectorXd rvECI = initialStateVec;

        UT::trans_model coorTransECI2MEE = [&mu](const VectorXd &satECI) -> VectorXd
        {
            // Perform coordinate transformation and return the result
            return eci2mee(satECI, mu);
        };
        // convert procNoiseCov from RIC to ECI first

        // calcualte the process noise cov first as the initialStateVec will be overwritten
        UT utECI2MEENoise(coorTransECI2MEE, false, 0, initialStateVec, rvCov, procNoiseCov, UT::sig_type::JU, 1);
        utECI2MEENoise(procNoiseCov);
        UT utECI2MEEStateCov(coorTransECI2MEE, false, 0, initialStateVec, rvCov, MatrixXd::Zero(dimState, dimState), UT::sig_type::JU, 1);
        utECI2MEEStateCov(initialStateVec, initialCov);
    }
}

int findClosestIndex(const VectorXd &tSec, double target)
{
    double minDiff = std::numeric_limits<double>::max();
    int closestIndex = -1;

    for (int i = 0; i < tSec.size(); i++)
    {
        double diff = std::abs(tSec[i] - target);
        if (diff < minDiff)
        {
            minDiff = diff;
            closestIndex = i;
        }
    }

    return closestIndex;
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
    // Convert string to lowercase
    for (char &c : initialStateType)
        c = tolower(c);
    VectorXd initialStateVec = initialState.initialStateVec;
    MatrixXd initialCov = initialState.initialCovarianceMat;
    double initialSkewness = initialState.initialSkewness;
    double initialKurtosis = initialState.initialKurtosis;
    // process noise covariance
    MatrixXd procNoiseCov = initialState.processNoiseCovarianceMat;
    const VectorXd groundStation = measMdl.groundStation;
    // Find the position of the dot (file type extension)
    size_t dotPos = measMdl.measFile.find_last_of('.');
    // Extract the last five characters without the file type extension
    string noradID = measMdl.measFile.substr(dotPos - 5, 5);

    // initGlobalVariables(initialState, suppFiles);
    initGlobalVariables(initialStateVec, initialCov, procNoiseCov, initialStateType, suppFiles);
    cout << "initialStateVec:\t\n"
         << initialStateVec << endl;
    cout << "initialCov:\t\n"
         << initialCov << endl;
    cout << "procNoiseCov:\t\n"
         << procNoiseCov << endl;

    // setup orbit propagator
    double leapSec = -getLeapSecond(convertMJD2Time_T(epoch.startMJD));
    orbitProp.setPropOption(forceModelsPropOpt);
    orbitProp.printPropOption();
    orbitProp.initPropagator(initialStateVec, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);
    // UKF state & measurement models
    double absErr = 1E-6;
    double relErr = 1E-6;
    DynamicModel::stf accMdl = accelerationModel;
    DynamicModel orbFun(accMdl, dimState, absErr, relErr);

    // measurement model
    UKF::meas_model h;
    HOUSE::meas_model hh;
    if (initialStateType == "eci")
    {
        h = [&groundStation](double t, const VectorXd &x) -> VectorXd
        {
            return measurementModel(t, x, groundStation);
        };

        hh = [&groundStation](double t, const VectorXd &x, const VectorXd &n)
            -> VectorXd
        {
            return measurementModel(t, x, groundStation) + n;
        };
    }
    else if (initialStateType == "mee")
    {
        h = [&groundStation](double t, const VectorXd &x) -> VectorXd
        {
            return measurementMEEModel(t, x, groundStation);
        };

        hh = [&groundStation](double t, const VectorXd &x, const VectorXd &n)
            -> VectorXd
        {
            return measurementMEEModel(t, x, groundStation) + n;
        };
    }

    // measurement noise covariance
    Matrix2d measNoiseCov;
    measNoiseCov << pow(measMdl.errorStd.rightAscensionErr * ARC_SEC, 2), 0,
        0, pow(measMdl.errorStd.declinationErr * ARC_SEC, 2);

    // read angular measurements from the exiting file
    string filename = measMdl.measFile;
    int headerLinesToSkip = 1; // Number of lines to skip as header
    MatrixXd matMeas = readCSV(filename, headerLinesToSkip);
    // cout << "measurement : \n"
    //      << matMeas << endl;
    int dimMeas = measMdl.dimMeas;
    int closestStartIndex = findClosestIndex(matMeas.col(6).array(), epoch.startMJD);
    int closestEndIndex = findClosestIndex(matMeas.col(6).array(), epoch.endMJD);
    // cout << closestStartIndex << "\t" << closestEndIndex << endl;
    int nRows = closestEndIndex - closestStartIndex + 1;
    VectorXd tSec = matMeas.col(6).array().segment(closestStartIndex, nRows);
    // time intervals relative to the first epoch in MJD
    tSec = (tSec.array() - tSec(0)) * 86400;
    // Extract two angular measurements in the last two columns
    MatrixXd angMeas = matMeas.block(closestStartIndex, matMeas.cols() - dimMeas, nRows, dimMeas);
    // Convert from degress into raidans
    angMeas = angMeas * (M_PI / 180.0);
    // Transpose of angMeas to ensure the measurement in a vector for each time epoch
    angMeas = angMeas.transpose().eval();
    // vector<string> measStrings({"angular measurements"});
    // EigenCSV::write(angMeas, measStrings, "angles.csv");
    // // test the prediction step only
    // MatrixXd angMeas = MatrixXd::Constant(nRows, dimMeas, 4 * M_PI);
    // angMeas = angMeas.transpose().eval();

    double dtMax = epoch.maxTimeStep;
    cout << "max time step:\t" << dtMax << endl;
    // Initialize UKF & CUT filters
    UKF ukf(orbFun, h, true, 0, dtMax, initialStateVec, initialCov, procNoiseCov, measNoiseCov, UKF::sig_type::JU, 1);
    UKF cut4(orbFun, h, true, 0, dtMax, initialStateVec, initialCov, procNoiseCov, measNoiseCov, UKF::sig_type::CUT4, 1);
    UKF cut6(orbFun, h, true, 0, dtMax, initialStateVec, initialCov, procNoiseCov, measNoiseCov, UKF::sig_type::CUT6, 1);
    UKF cut8(orbFun, h, true, 0, dtMax, initialStateVec, initialCov, procNoiseCov, measNoiseCov, UKF::sig_type::CUT8, 1);

    // HOUSE distributions for state
    HOUSE::Dist distXi(initialCov);
    distXi.mean = initialStateVec;
    distXi.skew.setConstant(dimState, initialSkewness);
    distXi.kurt.setConstant(dimState, initialKurtosis);

    // HOUSE distributions for state noise
    HOUSE::Dist distw(procNoiseCov);
    // HOUSE distributions for measurement noise
    HOUSE::Dist distn(measNoiseCov);

    string outputFile;
    MatrixXd runTimesMC(1, 4);
    Timer timer;

    if (filters.house)
    {
        cout << "\tHOUSE" << '\n';

        for (int j = 1; j <= filters.numTrials; j++)
        {
            cout << "HOUSE Trial " << j << endl;
            timer.tick();
            // Initialize HOUSE with different delta
            double delta = 0.1 / filters.numTrials * (j - 1);
            HOUSE house(orbFun, hh, dimMeas, 0, dtMax, distXi, distw, distn, delta);
            // HOUSE house(orbFun, hh, dimMeas, 0, dtMax, distXi, distw, distn, 0.05);
            house.run(tSec, angMeas);
            runTimesMC(0) = timer.tock();

            if (filters.numTrials == 1)
                outputFile = snrInfo.outDir + "/house_id_" + noradID + "_" + initialStateType + ".csv";
            else
                outputFile = snrInfo.outDir + "/house_id_" + noradID + "_" + initialStateType + "_" + to_string(j) + ".csv";
            house.save(outputFile, initialStateType);

            // Save Filter run times
            if (j == 1)
            {
                vector<string> filterStrings({"house"});
                string timeFile = snrInfo.outDir + "/run_times_house_id_" + noradID + "_" + initialStateType + ".csv";
                EigenCSV::write(runTimesMC.col(0), filterStrings, timeFile);
            }
        }
    }

    // UKF Filter
    if (filters.ukf)
    {
        cout << "\tUKF" << '\n';
        // VectorXd propStateVec = orbFun(0, tSec(1), initialStateVec, Vector3d::Zero());
        // cout << "initial state: \t" << setw(24) << setprecision(20) << initialStateVec << endl;
        // cout << "predicted state: \t" << setw(24) << setprecision(20) << propStateVec << endl;
        timer.tick();
        ukf.run(tSec, angMeas);
        runTimesMC(1) = timer.tock();

        outputFile = snrInfo.outDir + "/ukf_id_" + noradID + "_" + initialStateType + ".csv";
        ukf.save(outputFile, initialStateType);

        // Save Filter run times
        vector<string> filterStrings({"ukf"});
        string timeFile = snrInfo.outDir + "/run_times_ukf_id_" + noradID + "_" + initialStateType + ".csv";
        EigenCSV::write(runTimesMC.col(1), filterStrings, timeFile);
    }

    // CUT-4 Filter
    if (filters.cut4)
    {
        cout << "\tCUT-4" << '\n';
        timer.tick();
        cut4.run(tSec, angMeas);
        runTimesMC(2) = timer.tock();

        outputFile = snrInfo.outDir + "/cut4_id_" + noradID + "_" + initialStateType + ".csv";
        cut4.save(outputFile, initialStateType);

        // Save Filter run times
        vector<string> filterStrings({"cut4"});
        string timeFile = snrInfo.outDir + "/run_times_cut4_id_" + noradID + "_" + initialStateType + ".csv";
        EigenCSV::write(runTimesMC.col(2), filterStrings, timeFile);
    }

    // Cut-6 Filter
    if (filters.cut6)
    {
        cout << "\tCUT-6" << '\n';
        timer.tick();
        cut6.run(tSec, angMeas);
        runTimesMC(3) = timer.tock();

        outputFile = snrInfo.outDir + "/cut6_id_" + noradID + "_" + initialStateType + ".csv";
        cut6.save(outputFile, initialStateType);

        // Save Filter run times
        vector<string> filterStrings({"cut6"});
        string timeFile = snrInfo.outDir + "/run_times_cut6_id_" + noradID + "_" + initialStateType + ".csv";
        EigenCSV::write(runTimesMC.col(3), filterStrings, timeFile);
    }
}
