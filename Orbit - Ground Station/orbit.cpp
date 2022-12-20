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
// #include "jplEph.hpp"
// #include "config.hpp"

// filter headers
#include "house.hpp"
#include "ukf.hpp"
#include "dyn.hpp"
#include "filter_aux.hpp"
#include "pearsonator.hpp"
#include "eigen_csv.hpp"
#include "timer.hpp"

#define MJD_EPOCH 59945


// #define MU 6.6743E-11 * 5.972E24
#define NUM_ORBITS 10

#define R_EARTH 6371E3
#define DEG M_PI / 180
#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)

// 4.3378   -2.0118   -4.2101
#define GROUND_STATION_X 4.33781E6
#define GROUND_STATION_Y -2.01181E6
#define GROUND_STATION_Z -4.21011E6

#define BB 1E-3

#define THRESHOLD 1E-20

using namespace Eigen;
using namespace std;



int is_visible(VectorXd& X){
    Vector3d x = X.head(3);
    Vector3d rs;
    rs << GROUND_STATION_X, GROUND_STATION_Y, GROUND_STATION_Z;
    Vector3d rel = x - rs;

    
    if (rs.dot(rel) >= 0){
        return 1;
    }
    else{
        return 0;
    }
}

string filter_file(bool gauss, const string& filter, int trial) {
    string filename = "out/";
    filename += filter;
    filename += "_est_";
    filename += (gauss ? "gauss_" : "pearson_");
    filename += to_string(trial);
    filename += ".csv";
    return filename;
}

void run(bool gauss) {
    DynamicModel::stf g = [] (double t, const VectorXd& X, const VectorXd& fd)
        -> VectorXd {
            VectorXd Xf(6);
            VectorXd r(3);
            VectorXd v(3);
            r = X.head(3);
            v = X.tail(3);

            // (d/dt) r = v
            Xf.head(3) = v;

            // (d/dt) v = -mu/|r|^3 * r
            Xf(3) = -MU / pow(r.norm(), 3) * r(0);
            Xf(4) = -MU / pow(r.norm(), 3) * r(1);
            Xf(5) = -MU / pow(r.norm(), 3) * r(2);
            return Xf;
        };
    // errors were previously 1E-9
    DynamicModel f(g, 6, 1E-9, 1E-9);

    UKF::meas_model h = [] (double t, const VectorXd& X)
        -> VectorXd {
            Vector2d z;
            Vector3d p;
            VectorXd rsECEF = VectorXd::Zero(6);
            VectorXd rsECI = VectorXd::Zero(6);

        	erp_t *erpt;

            // ground station
            rsECEF << GROUND_STATION_X, GROUND_STATION_Y, GROUND_STATION_Z;

            // transform ground station from 
            double leapSec = 32;
	        double mjdUTC = 53300;

            double *erpv;
            geterp_from_utc(erpt, leapSec, mjdUTC, erpv);

            double dUT1_UTC = erpv[2];
            double dUTC_TAI = -(19 + leapSec);
            double xp = erpv[0];
            double yp = erpv[1];
            double lod = erpv[3];
            IERS iersInstance;
            iersInstance.Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
            // double mjdTT = mjdUTC + iersInstance.TT_UTC(mjdUTC) / 86400;

            eci2ecefVec_sofa(mjdUTC + t, iersInstance, rsECEF, rsECI);
            // end transformation code


            p = X.head(3) - rsECI.head(3);

            z(0) = atan2(p(1), -p(0));
            z(1) = asin(p(2)/p.norm());


            return z;
        };

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
    Matrix2d Pnn;
    Pnn <<  0.4 * ARC_SEC,  0,
            0,              0.07 * ARC_SEC;

    double skeww, skewn, kurtw, kurtn, stdx0, stdv0, skew0, kurt0;

    if (gauss) {

        skeww = 0;
        kurtw = 3;

        skewn = 0;
        kurtn = 3;

        skew0 = 0;
        kurt0 = 3;

     } else {

        skeww = 1;
        kurtw = 30;

        skewn =  -1;
        kurtn =  30;

        skew0 = 1;
        kurt0 = 30;

    }


    // standard deviation for initial state (assume to be the same as measurement)
    stdx0 = 0.4 * ARC_SEC * R_EARTH;
    stdv0 = 0.07 * ARC_SEC * R_EARTH;

    VectorXd X0m(6), X0std(6), X0skew(6), X0kurt(6);

    // initializes X0m with these constants
    X0m  << 5052.160, -2343.109, -4903.437, 0.477, 6.671, -2.888;
    X0m *= 1000;

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

    noisemaker genw  = [&] () -> double {return gauss ? genw_g (mt) : genw_p (mt);};
    noisemaker genn  = [&] () -> double {return gauss ? genn_g (mt) : genn_p (mt);};
    noisemaker genx0 = [&] () -> double {return gauss ? genx0_g(mt) : genx0_p(mt);};
    noisemaker genv0 = [&] () -> double {return gauss ? genv0_g(mt) : genv0_p(mt);};

    UKF::cut_dir = "../CUT/";

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

    // time steps of 1s for N orbits
    double dt = 30;
    double tmin = 6360 * NUM_ORBITS;

    int trials = 1;

    vector<string> header(7);
    header[0] = "t";
    header[1] = "x";
    header[2] = "y";
    header[3] = "z";
    header[4] = "vx";
    header[5] = "vy";
    header[6] = "vz";

    MatrixXd run_times(trials, 3);
    Timer timer;

    for (int j = 1; j <= trials; j++) {

        cout << "Running trial " << j << endl;

        house.reset(0, distx0);
        ukf.reset(0, X0m, Pxx0);
        cut4.reset(0, X0m, Pxx0);
        cut6.reset(0, X0m, Pxx0);

        vector<VectorXd> Xtru;
        double time;
        // vector<double> times;
        do {
            cout << "       Trying" << endl;
            time = 0;
            VectorXd X(6);
            X = X0m;
            // MONTE CARLO SIMULATION GENERATION
            // for (int i = 0; i < 3; i++) {
            //     X(i)   += genx0();
            //     X(i+3) += genv0();
            // }
            Xtru.clear();
            
            Xtru.push_back(X);
            // times.push_back(time);
            do {
                Vector3d fd;
                fd << genw(), genw(), genw();
                X = f(time, time+dt, X, fd);
                time += dt;
                
                Xtru.push_back(X);
                // times.push_back(time);
            
            } while (time < tmin);
        } while  (0);//(time < tmin);

        int steps = Xtru.size();
        // Map<VectorXd> t(&times[0], times.size());

        // linear spaced times
        VectorXd t;
        t.setLinSpaced(steps, 0, (steps-1)*dt);




        MatrixXd table(steps, 7);
        table.col(0) = t;
        for (int k = 0; k < steps; k++)
            table.row(k).tail(6) = Xtru[k];

        string xtrufile = "out/tru_";
        xtrufile += (gauss ? "gauss_" : "pearson_");
        xtrufile += to_string(j);
        xtrufile += ".csv";

        EigenCSV::write(table, header, xtrufile);

        // matrix of measurements
        MatrixXd Ztru(2, steps);
        for (int k = 0; k < steps; k++) {
            if (is_visible(Xtru[k])){
                Ztru.col(k) = h(t(k), Xtru[k]);
                Ztru(0,k) += genn();
                Ztru(1,k) += genn();
            }
            else {
                Ztru(0,k) = NO_MEASUREMENT; // std::numeric_limits<double>::quiet_NaN();
                Ztru(1,k) = NO_MEASUREMENT; // std::numeric_limits<double>::quiet_NaN();

            }
        }


        cout << "   HOUSE" << endl;
        timer.tick();
        house.run(t, Ztru);
        run_times(j-1, 0) = timer.tock();
        house.save(filter_file(gauss, "house", j));

        cout << "   UKF" << endl;
        timer.tick();
        ukf.run(t, Ztru);
        run_times(j-1, 1) = timer.tock();
        ukf.save(filter_file(gauss, "ukf", j));

        cout << "   CUT-4" << endl;
        timer.tick();
        cut4.run(t, Ztru);
        run_times(j-1, 2) = timer.tock();
        cut4.save(filter_file(gauss, "cut4", j));

        // cout << "   CUT-6" << endl;
        // timer.tick();
        // cut6.run(t, Ztru);
        // run_times(j-1, 3) = timer.tock();
        // cut6.save(filter_file(gauss, "cut6", j));

    }

    vector<string> filters;
    filters.push_back("house");
    filters.push_back("ukf");
    filters.push_back("cut4");
    // filters.push_back("cut6");

    string time_file = "out/run_times_";
    time_file += (gauss ? "gauss" : "pearson");
    time_file += ".csv";

    EigenCSV::write(run_times, filters, time_file);

}

int main() {

    cout << "-- Gaussian Distributions --" << endl;
    run(true);

    cout << "-- Pearson Type IV Distributions --" << endl;
    run(false);

    return 0;

}
