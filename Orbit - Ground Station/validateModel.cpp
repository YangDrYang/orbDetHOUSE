#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

// model headers
#include "auxillaryData.hpp"
#include "satRefSys.hpp"
#include "common.hpp"
#include "forceModels.hpp"
#include "jplEph.hpp"
#include "jpl_eph.hpp"
#include "config.hpp"

// filter headers
#include "house.hpp"
#include "ukf.hpp"
#include "dyn.hpp"
#include "filter_aux.hpp"
#include "pearsonator.hpp"
#include "eigen_csv.hpp"
#include "timer.hpp"

// #define MJD_EPOCH_START 5.757099980324074e+04     // days
// #define MJD_EPOCH_END   5.757212410879634e+04   //days
#define TIME_STEP       60                      // seconds

#define R_EARTH 6371E3
#define DEG M_PI / 180
#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)

double MJD_EPOCH_START;
double MJD_EPOCH_END;

erp_t erpt;
double leapSec;
double mjdUTC;
IERS iersInstance;
ForceModels forceModelsOpt = {};
EGMCoef egm;
void *pJPLEph;
VectorXd rvPhiS(6);

Propagator orbitProp;

// prototypes
VectorXd orbitModel(double t, const VectorXd& X);
void initEGMCoef(string filename);
void initGlobalVariables(VectorXd rvECI);


// space seperated input parameters
// MJD_START MJD_END X0 X1 X2 X3 X4 X5 
int main(int argc, char *argv[]){
    if (argc != 9) {
        cerr << "Wrong number of arguments" << endl;
        exit(1);
    }
    // initialise times
    MJD_EPOCH_START = stod(argv[1]);
    MJD_EPOCH_END = stod(argv[2]);

    // initial satillite state
    VectorXd initialState(6);
    for (int i = 3; i < 9; i++){
        initialState(i - 3) = stod(argv[i]);
    }

    ForceModels options = {};
    // Set specific force model options
    options.third_body_attraction = true;
    options.third_body_sun = true;
    options.third_body_moon = true;
    options.solid_earth_tide = false;
    options.ocean_tide_loading = false;
    options.solar_radiation_pressure = true; // should be true
    options.relativity_effect = true;  // should be true
    options.srpCoef = 1.23;
    options.srpArea = 2 * M_PI  * pow(0.974 / 2, 2) * sin(0) + 1.034 * 0.132;
    options.satMass  = 61.14;


    vector<string> header(7);
    header[0] = "t";
    header[1] = "x";
    header[2] = "y";
    header[3] = "z";
    header[4] = "vx";
    header[5] = "vy";
    header[6] = "vz";


    initGlobalVariables(initialState);

    // setup propagator
    orbitProp.setPropOption(forceModelsOpt);
    orbitProp.initPropagator(initialState, rvPhiS, mjdUTC, leapSec, &erpt, egm, pJPLEph);
    
    // create results matrix
    int length = floor((MJD_EPOCH_END - MJD_EPOCH_START) / (TIME_STEP / 86400.0)) + 1;
    
    MatrixXd results(length, 7);

    VectorXd xECEF = VectorXd::Zero(6);
    VectorXd xECI = VectorXd::Zero(6);
    double time = 0;

    DynamicModel::stf g = [] (double t, const VectorXd& X, const VectorXd& fd)
        -> VectorXd {
            VectorXd Xf(6);
            Vector3d r;
            Vector3d v;
            VectorXd rvECI(6);
            rvECI = X;

            Matrix3d mECI2ECEF = Matrix3d::Identity();
            Matrix3d mdECI2ECEF = Matrix3d::Identity();

            eci2ecef_sofa(mjdUTC + t/86400, iersInstance, mECI2ECEF, mdECI2ECEF);

            Vector3d acceleration;
            orbitProp.updPropagator(mjdUTC + t/86400, leapSec, &erpt);

            acceleration = orbitProp.calculateAcceleration(X.head(3), X.tail(3), mECI2ECEF);
            // acceleration = -MU / pow(X.head(3).norm(), 3) * X.head(3);
            // (d/dt) r = v

            // cout << acceleration << endl;
            Xf.head(3) = X.tail(3);
            Xf.tail(3) = acceleration;

            // (d/dt) v = -mu/|r|^3 * r

            return Xf;
        };
    // errors were previously 1E-9
    DynamicModel f(g, 6, 1E-9, 1E-9);

    VectorXd result(7);
    // xECEF = initialState;
    result << MJD_EPOCH_START, initialState;
    results.row(0) = result;
    xECEF = initialState;
    ecef2eciVec_sofa(MJD_EPOCH_START, iersInstance, xECEF, xECI);
    time += TIME_STEP;
    // pass every timestep into model

    for (int i = 1; i < length; i++) {
        VectorXd result(7);
        Vector3d fd;
        fd << 0, 0, 0;

        xECI = f(time, time + TIME_STEP, xECI, fd);

        eci2ecefVec_sofa(MJD_EPOCH_START + time / 86400, iersInstance, xECI, xECEF);
        result << MJD_EPOCH_START + time / 86400, xECEF;
        
        results.row(i) = result; 

        time += TIME_STEP;
        // cout << result << '\n';
    }
    cout << "written to file" << endl;
    EigenCSV::write(results, header, "validate_data.csv");

}

VectorXd orbitModel(double t, const VectorXd& X){
    VectorXd Xf(6);
    Vector3d r;
    Vector3d v;
    VectorXd rvECI(6);
    rvECI = X;

    Matrix3d mECI2ECEF = Matrix3d::Identity();
    Matrix3d mdECI2ECEF = Matrix3d::Identity();

    eci2ecef_sofa(mjdUTC + t, iersInstance, mECI2ECEF, mdECI2ECEF);

    Vector3d acceleration;
    orbitProp.updPropagator(mjdUTC + t/86400, leapSec, &erpt);

    acceleration = orbitProp.calculateAcceleration(X.head(3), X.tail(3), mECI2ECEF);
    // acceleration = -MU / pow(X.head(3).norm(), 3) * X.head(3);
    // (d/dt) r = v
    Xf.head(3) = X.tail(3);
    Xf.tail(3) = acceleration;

    // (d/dt) v = -mu/|r|^3 * r

    return Xf;
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
    mjdUTC = MJD_EPOCH_START;

    double erpv[0] = {};
    geterp_from_utc(&erpt, leapSec, mjdUTC, erpv);

    double dUT1_UTC = erpv[2];
    double dUTC_TAI = -(19 + leapSec);
    double xp = erpv[0];
    double yp = erpv[1];
    double lod = erpv[3];
    
    iersInstance.Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
    // double mjdTT = mjdUTC + iersInstance.TT_UTC(mjdUTC) / 86400;

    pJPLEph = jpl_init_ephemeris("./unxp2000.405", nullptr, nullptr);

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
