#include "house.hpp"
#include "ukf.hpp"
#include "dyn.hpp"
#include "filter_aux.hpp"
#include "pearsonator.hpp"
#include "timer.hpp"
#include "eigen_csv.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

using namespace Eigen;
using namespace std;

void run(bool gauss) {

    Matrix3d I;

    I << 500+400, 0, 0,
         0, 500+300, 0,
         0, 0, 400+300;

    DynamicModel::stf g = [&I] (double t, const Vector3d& w, const Vector3d& Td)
        -> VectorXd {
            return I.inverse() * (Td - w.cross(I * w));
        };

    DynamicModel f(g, 3, 1E-9, 1E-9);

    UKF::meas_model h = [] (double t, const Vector3d& w)
        -> VectorXd {
            return w.head(1);
        };

    HOUSE::meas_model hh = [] (double t, const VectorXd& w, const VectorXd& v)
        -> VectorXd {
            return w.head(1) + v;
        };

    double stdx0, stdw, stdn, dt;

    stdx0 = 0.01;
    stdw  = 0.001;
    stdn  = 0.001;

    dt = 0.1;

    Matrix3d Pxx0, Pww;
    MatrixXd Pnn(1,1);

    Pxx0 = Matrix3d::Identity() * stdx0 * stdx0;
    Pww  = Matrix3d::Identity() * stdw  * stdw;
    Pnn(0,0) = stdn * stdn;

    Vector3d wm0;
    wm0 << 0, 0, 0;

    Vector3d w0;

    int trials = 100;
    int steps = 6001;

    VectorXd t;
    t.setLinSpaced(steps, 0, (steps-1)*dt);

    vector<vector<VectorXd>> xtru, xest_ukf, xest_cut4, xest_cut6, xest_cut8,
        xest_house;

    double skew0, skeww, skewn, kurt0, kurtw, kurtn;

    if (gauss) {

        skew0 = 0;
        kurt0 = 3;

        skeww = 0;
        kurtw = 3;

        skewn = 0;
        kurtn = 3;

    } else {

        skew0 = -1;
        kurt0 = 30;

        skeww = -1;
        kurtw = 30;

        skewn = -1;
        kurtn = 30;

    }

    Pearsonator::TypeIV gen0_p(0, stdx0, skew0, kurt0), genw_p(0, stdw, skeww, kurtw),
        genn_p(0, stdn, skewn, kurtn);

    normal_distribution<double> gen0_g(0, stdx0), genw_g(0, stdw), genn_g(0, stdn);

    mt19937_64 mt(0);

    typedef function<double(void)> noisemaker;

    noisemaker gen0 = [&] () -> double {return gauss ? gen0_g(mt) : gen0_p(mt);};
    noisemaker genw = [&] () -> double {return gauss ? genw_g(mt) : genw_p(mt);};
    noisemaker genn = [&] () -> double {return gauss ? genn_g(mt) : genn_p(mt);};

    UKF::cut_dir = "../CUT/";

    UKF ukf (f, h, false, 0, wm0, Pxx0, Pww, Pnn, UKF::sig_type::JU,   1);
    UKF cut4(f, h, false, 0, wm0, Pxx0, Pww, Pnn, UKF::sig_type::CUT4, 1);
    UKF cut6(f, h, false, 0, wm0, Pxx0, Pww, Pnn, UKF::sig_type::CUT6, 1);
    UKF cut8(f, h, false, 0, wm0, Pxx0, Pww, Pnn, UKF::sig_type::CUT8, 1);

    HOUSE::Dist distx0(Pxx0), distw(Pww), distn(Pnn);

    distx0.mean = wm0;

    distx0.skew.setConstant(skew0);
    distx0.kurt.setConstant(kurt0);

    distw.skew.setConstant(skeww);
    distw.kurt.setConstant(kurtw);

    distn.skew.setConstant(skewn);
    distn.kurt.setConstant(kurtn);

    HOUSE house(f, hh, 1, 0, distx0, distw, distn, 0);

    MatrixXd run_times(trials,5);
    Timer timer;

    for (int j = 0; j < trials; j++) {

        ukf. reset(0, wm0, Pxx0);
        cut4.reset(0, wm0, Pxx0);
        cut6.reset(0, wm0, Pxx0);
        cut8.reset(0, wm0, Pxx0);
        house.reset(0, distx0);

        cout << "Running Trial " << j+1 << endl;

        MatrixXd T(3, steps), v(1, steps);
        for (int k = 0; k < steps; k++) {
            v(k) = genn();
            T(0,k) = genw();
            T(1,k) = genw();
            T(2,k) = genw();
        }

        w0 = wm0;
        w0(0) += gen0();
        w0(1) += gen0();
        w0(2) += gen0();

        vector<VectorXd> xtruk;
        VectorXd w = w0;
        for (int k = 0; k < steps; k++) {
            xtruk.push_back(w);
            if (k < steps-1)
                w = f(t(k), t(k+1), xtruk.back(), T.col(k));
        }
        xtru.push_back(xtruk);

        MatrixXd Z(1, steps);
        for (int k = 0; k < steps; k++)
            Z.col(k) = h(t(k), xtruk[k]) + v.col(k);

        cout << "   HOUSE" << endl;
        timer.tick();
        house.run(t, Z);
        run_times(j, 0) = timer.tock();
        vector<VectorXd> xest_house_trial;
        for (int k = 0; k < steps; k++)
            xest_house_trial.push_back(house.distx[k].mean);
        xest_house.push_back(xest_house_trial);

        cout << "   UKF" << endl;
        timer.tick();
        ukf.run(t, Z);
        run_times(j, 1) = timer.tock();
        xest_ukf.push_back(ukf.xest);

        cout << "   CUT4" << endl;
        timer.tick();
        cut4.run(t, Z);
        run_times(j, 2) = timer.tock();
        xest_cut4.push_back(cut4.xest);

        cout << "   CUT6" << endl;
        timer.tick();
        cut6.run(t, Z);
        run_times(j, 3) = timer.tock();
        xest_cut6.push_back(cut6.xest);

        cout << "   CUT8" << endl;
        timer.tick();
        cut8.run(t, Z);
        run_times(j, 4) = timer.tock();
        xest_cut8.push_back(cut8.xest);

    }

    string dist = gauss ? "gauss" : "pearson";

    save_rmse(t, xtru, xest_house, "out/house_rmse_" + dist + ".csv");
    save_rmse(t, xtru, xest_ukf,   "out/ukf_rmse_"   + dist + ".csv");
    save_rmse(t, xtru, xest_cut4,  "out/cut4_rmse_"  + dist + ".csv");
    save_rmse(t, xtru, xest_cut6,  "out/cut6_rmse_"  + dist + ".csv");
    save_rmse(t, xtru, xest_cut8,  "out/cut8_rmse_"  + dist + ".csv");

    save_abs_err_lump(t, xtru, xest_house, "out/house_err_" + dist + ".csv");
    save_abs_err_lump(t, xtru, xest_ukf,   "out/ukf_err_"   + dist + ".csv");
    save_abs_err_lump(t, xtru, xest_cut4,  "out/cut4_err_"  + dist + ".csv");
    save_abs_err_lump(t, xtru, xest_cut6,  "out/cut6_err_"  + dist + ".csv");
    save_abs_err_lump(t, xtru, xest_cut8,  "out/cut8_err_"  + dist + ".csv");

    vector<string> filters;
    filters.push_back("house");
    filters.push_back("ukf");
    filters.push_back("cut4");
    filters.push_back("cut6");
    filters.push_back("cut8");

    string time_file = "out/run_time_";
    time_file += dist;
    time_file += ".csv";

    EigenCSV::write(run_times, filters, time_file);

}

int main() {

    cout << "--- Gaussian Distributions ---" << endl;
    run(true);

    cout << "--- Pearson Type IV Distributions ---" << endl;
    run(false);

    return 0;

}

