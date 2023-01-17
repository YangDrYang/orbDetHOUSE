#include "filter_testing.hpp"

#define R_EARTH 6371E3
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
void run(bool gauss) {


}

VectorXd accelerationModel(double t, const VectorXd& X, const VectorXd& fd) {
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

VectorXd measurementModel(double t, const VectorXd& X){
    Vector3d z, p; 
    VectorXd rsECEF(6); // = VectorXd::Zero(6);
    VectorXd rsECI = VectorXd::Zero(6);

    // ground station
    rsECEF << 4.33781e+6, -2.01181e+6, -4.21011e+6, 0, 0, 0;

    ecef2eciVec_sofa(epoch.startMJD + t, iersInstance, rsECEF, rsECI);
    // end transformation code
    p = X.head(3) - rsECI.head(3);

    if (p.dot(rsECI.head(3)) >= 0){
        // range
        z(0) = p.norm();
        // azimuth angle
        z(1) = atan2(p(1), -p(0));
        // elevation angle
        z(2) = asin(p(2)/z(0));
    }
    // not visible
    else {
        z(0) = NO_MEASUREMENT;
        z(1) = NO_MEASUREMENT;
        z(2) = NO_MEASUREMENT;
    }
    return z;
}

int main(int argc, char *argv[]){
    using namespace Eigen;
    using namespace std;
    
    string configFilename;
    // input config .yaml file
    switch (argc) {
    case 1:     // Read default file
        configFilename = DEFAULT_CONFIG_FILENAME;
        break;
    case 2:     // Read supplied file
        configFilename = argv[1];
        break;
    default:    // wrong number of args
        cerr << "Accepts up to 1 argument (.yaml input file).\nIf no argument given, the default config file will be read (config.yaml)" << endl;
        exit(1);
        break;
    }
    
    struct Errors errorStd;
    struct Filters filters;
    int numTrials;
    // read parameter/settings from config file
    readConfigFile(configFilename, forceModelsOpt, epoch, initialState, groundStation, filters, numTrials, errorStd);

    // initialise
    initGlobalVariables(initialState);

    // setup orbit propagator
    orbitProp.setPropOption(forceModelsOpt);
    orbitProp.initPropagator(initialState, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);
    
    


    // ---------- START ----------
        // initialize global variables
    bool gauss = true;
    DynamicModel::stf g = &accelerationModel;
    // errors were previously 1E-9
    DynamicModel f(g, 6, 1E-3, 1E-3);

    UKF::meas_model h = &measurementModel;
    HOUSE::meas_model hh = [h] (double t, const VectorXd& X, const VectorXd& n)
        -> VectorXd {
            return h(t, X) + n;
        };

    double stdw, stdn;
    stdw = 0; // disturbing force not currently being considered so this won't affect anthing
    // sigma for both measurement values, ideally <0.4, 0.07> arc seconds.
    stdn = ARC_MIN;

    Matrix3d Pww = Matrix3d::Identity() * stdw * stdw;

    // Pnn is measurement noise matrix
    Matrix3d Pnn;
    Pnn <<  errorStd.rangeErr,  0, 0,
            0, errorStd.azimuthErr * ARC_SEC, 0,
            0, 0, errorStd.rangeErr * ARC_SEC;

    double skeww, skewn, kurtw, kurtn, stdx0, stdv0, skew0, kurt0;
    skeww = 0;
    kurtw = 3;

    skewn = 0;
    kurtn = 3;

    skew0 = 0;
    kurt0 = 3;

    // standard deviation for initial state (assume to be the same as measurement)
    stdx0 = errorStd.rangeErr;
    stdv0 = errorStd.rangeErr;

    VectorXd X0m(6), X0std(6), X0skew(6), X0kurt(6);

    // initializes X0m with these constants
    X0m = initialState;

    // sets front half to stdx0
    X0std.head(3).setConstant(stdx0);
    // sets the tail to stdv0
    X0std.tail(3).setConstant(stdv0);
    // initializes vector to given values.
    X0skew.setConstant(skew0);
    X0kurt.setConstant(kurt0);

    MatrixXd Pxx0 = X0std.array().square().matrix().asDiagonal();

    Pearsonator::TypeIV genw_p (0, stdw,  skeww, kurtw),
                        genn_p (0, stdn,  skewn, kurtn),
                        genx0_p(0, stdx0, skew0, kurt0),
                        genv0_p(0, stdv0, skew0, kurt0);

    normal_distribution<double> genw_g (0, stdw),
                                genn_g (0, stdn),
                                genx0_g(0, stdx0),
                                genv0_g(0, stdv0);

    mt19937_64 mt(0);

    typedef function<double(void)> noisemaker;
    noisemaker genw  = [&] () -> double {return genw_g (mt);};
    noisemaker genn  = [&] () -> double {return genn_g (mt);};
    noisemaker genx0 = [&] () -> double {return genx0_g(mt);};
    noisemaker genv0 = [&] () -> double {return genv0_g(mt);};


    UKF ukf (f, h, false, 0, X0m, Pxx0, Pww, Pnn, UKF::sig_type::JU,   1);
    UKF cut4(f, h, false, 0, X0m, Pxx0, Pww, Pnn, UKF::sig_type::CUT4, 1);
    UKF cut6(f, h, false, 0, X0m, Pxx0, Pww, Pnn, UKF::sig_type::CUT6, 1);

    HOUSE::Dist distx0(Pxx0), distw(Pww), distn(Pnn);

    distx0.mean = X0m;
    distx0.skew = X0skew;
    distx0.kurt = X0kurt;
    distw.skew.setConstant(skeww);
    distw.kurt.setConstant(kurtw);
    distn.skew.setConstant(skewn);
    distn.kurt.setConstant(kurtn);

    HOUSE house(f, hh, 2, 0, distx0, distw, distn, 0);
    vector<string> header({"t", "x", "y", "z", "vx", "vy", "vz"});
    MatrixXd run_times(numTrials, 3);
    Timer timer;

    initGlobalVariables(X0m);


    // perform trials
    for (int j = 1; j <= numTrials; j++) {
        cout << "Trial " << j << endl;
        house.reset(0, distx0);
        ukf.reset(0, X0m, Pxx0);
        cut4.reset(0, X0m, Pxx0);
        cut6.reset(0, X0m, Pxx0);

        // Generate true vectors
        // TODO: add montecarlo simulation
        vector<VectorXd> trueData;
        trueData.clear();
        orbitProp.setPropOption(forceModelsOpt);
        orbitProp.initPropagator(X0m, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph);

        
        VectorXd X = initialState;
        double time = 0, dt = epoch.timeStep;        
        trueData.push_back(X);
        while (time < (epoch.endMJD - epoch.startMJD) * 86400) {
            Vector3d fd;
            fd << 0, 0, 0;
            X = f(time, time + epoch.timeStep / 86400, X, fd);
            time += epoch.timeStep;
            trueData.push_back(X);        
        } 

        int steps = trueData.size();
        // linear spaced times
        VectorXd t;
        t.setLinSpaced(steps, 0, (steps - 1)*dt);

        // Save true vectors
        MatrixXd tableTrue(steps, 7);
        tableTrue.col(0) = t;
        for (int k = 0; k < steps; k++){
            tableTrue.row(k).tail(6) = trueData[k];
        };
        string xtrufile = "out/tru_";
        xtrufile += to_string(j);
        xtrufile += ".csv";
        EigenCSV::write(tableTrue, header, xtrufile);

        // Generate Measurement vectors
        vector<string> headerMeasure({"t", "alpha", "delta"});
        MatrixXd tableMeasure(steps, 4);
        int numMeasure = 0;
        MatrixXd Ztru(3, steps);
        
        // for each time step
        for (int k = 0; k < steps; k++) {
            // planar condition for visibility 
            // generate true measurement
            Ztru.col(k) = h(t(k), trueData[k]);

            // if !NO_MEASUREMENT
            if (Ztru(0, k) != NO_MEASUREMENT){
                // corrupt measurement with noise
                Ztru(0,k) += genn();
                Ztru(1,k) += genn();
                // store corrupted measurement & time
                tableMeasure(numMeasure, 0) = t(k);
                tableMeasure(numMeasure, 1) = Ztru(0,k);
                tableMeasure(numMeasure, 2) = Ztru(1,k);
                numMeasure++;
            }
        }
        // save corrupted measurment data
        string xmeasurefile = "out/meas_";
        xmeasurefile += to_string(j);
        xmeasurefile += ".csv";
        EigenCSV::write(tableMeasure.topRows(numMeasure), headerMeasure, xmeasurefile);

        string outputFile;
        if (filters.house){
            cout << "\tHOUSE" << '\n';
            orbitProp.initPropagator(X0m, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            house.run(t, Ztru);
            outputFile = "out";
            outputFile += "/house_";
            outputFile += to_string(j);
            outputFile += ".csv";
            house.save(outputFile);

        }

        // UKF Filter
        if (filters.ukf){
            cout << "UKF" << '\n';
            orbitProp.initPropagator(X0m, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            ukf.run(t, Ztru);
            outputFile = "out";
            outputFile += "/ukf_";
            outputFile += to_string(j);
            outputFile += ".csv";
            house.save(outputFile);
        }

        // CUT-4 Filter
        if (filters.cut4){
            cout << "CUT-4" << '\n';
            orbitProp.initPropagator(X0m, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            cut4.run(t, Ztru);
            outputFile = "out";
            outputFile += "/cut4_";
            outputFile += to_string(j);
            outputFile += ".csv";
            house.save(outputFile);
        }

        // Cut-6 Filter
        if (filters.cut6){
            cout << "CUT-6" << '\n';
            orbitProp.initPropagator(X0m, rvPhiS, epoch.startMJD, leapSec, &erpt, egm, pJPLEph); // reset propagator
            cut6.run(t, Ztru);
            outputFile = "out";
            outputFile += "/cut6_";
            outputFile += to_string(j);
            outputFile += ".csv";
            house.save(outputFile);
        }
    }

    // // Save Filter run times
    // vector<string> filterStrings;
    // filterStrings.push_back("house");
    // filterStrings.push_back("ukf");
    // // filters.push_back("cut4");
    // // filters.push_back("cut6");

    // string time_file = "out/run_times_";
    // time_file += (gauss ? "gauss" : "pearson");
    // time_file += ".csv";

    // EigenCSV::write(run_times, filterStrings, time_file);
    
}



Eigen::MatrixXd generateTrueResults(DynamicModel& f, struct EpochInfo epoch, VectorXd initialState){
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
    for (int i = 1; i < length; i++){
        X = f(time, time + step, X, W);
        time += step;
        result << time, X;
        results.row(i) = result;
        
    }

    return results;
}

void readConfigFile(string fileName, ForceModels& options, struct EpochInfo& epoch, Eigen::VectorXd& initialState,
                    Eigen::VectorXd& groundStation, struct Filters& filters, int& numTrials, struct Errors& errorStd){
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

Eigen::VectorXd stdVec2EigenVec(const std::vector<double>& stdVec){
    Eigen::VectorXd eigenVec(stdVec.size());
    for (int i = 0; i < 6; i++)
    {
        eigenVec(i) = stdVec[i];
    }
    return eigenVec;
}

void initEGMCoef(string filename){
    ifstream file(filename);
    string line;
    int n, m;
    double cmn, smn;
    

    while(getline(file, line)){
        // line structure:
        // m        n       Cnm     Snm     0      0 
        istringstream buffer(line);
        buffer >> m >> n >> cmn >> smn;

        // gravity model degree
        if (m < GRAVITY_DEG_M && n < GRAVITY_DEG_M){
            egm.cmn(m, n) = cmn;
            egm.smn(m, n) = smn;
        }     
    }


}

void initGlobalVariables(VectorXd rvECI) {
// START MOD
    initEGMCoef("GGM03S.txt");
    erpt = {.n = 14};

    readerp("cof19037.erp", &erpt);

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
