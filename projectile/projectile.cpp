#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "house.hpp"
#include "ukf.hpp"
#include "dyn.hpp"
#include "filter_aux.hpp"
#include "pearsonator.hpp"
#include "eigen_csv.hpp"
#include "timer.hpp"

#define G_0 9.80665

#define DEG M_PI / 180
#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)

#define BB 1E-3

#define THRESHOLD 1E-20

using namespace Eigen;
using namespace std;

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
            double d = -BB * X.tail(3).norm();
            VectorXd Xd(6);
            Xd.head(3) = X.tail(3);
            Xd.tail(3) = d * X.tail(3) + fd;
            Xd(5) -= G_0;
            return Xd;
        };

    DynamicModel f(g, 6, 1E-9, 1E-9);

    UKF::meas_model h = [] (double t, const VectorXd& X)
        -> VectorXd {
            double r = X.head(2).norm();
            Vector2d z;
            if (r > THRESHOLD) {
                z(0) = atan2(X(1), -X(0));
                z(1) = atan2(X(2), r);
            } else {
                z(0) = 0;
                z(1) = M_PI / 2;
            }
            return z;
        };

    HOUSE::meas_model hh = [h] (double t, const VectorXd& X, const VectorXd& n)
        -> VectorXd {
            return h(t, X) + n;
        };

    double stdw, stdn;
    stdw = 0.01;
    stdn = ARC_MIN;

    Matrix3d Pww = Matrix3d::Identity() * stdw * stdw;
    Matrix2d Pnn = Matrix2d::Identity() * stdn * stdn;

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

        skewn = -1;
        kurtn =  30;

        skew0 = 1;
        kurt0 = 30;

    }

    stdx0 = 250;
    stdv0 = 100;

    VectorXd X0m(6), X0std(6), X0skew(6), X0kurt(6);

    X0m  << 1000, 1000, 0, 500, 0, 500;

    X0std.head(3).setConstant(stdx0);
    X0std.tail(3).setConstant(stdv0);
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

    double dt = 0.1;
    double tmin = 1;

    int trials = 100;

    vector<string> header(7);
    header[0] = "t";
    header[1] = "x";
    header[2] = "y";
    header[3] = "z";
    header[4] = "vx";
    header[5] = "vy";
    header[6] = "vz";

    MatrixXd run_times(trials, 4);
    Timer timer;

    for (int j = 1; j <= trials; j++) {

        cout << "Running trial " << j << endl;

        house.reset(0, distx0);
        ukf.reset(0, X0m, Pxx0);
        cut4.reset(0, X0m, Pxx0);
        cut6.reset(0, X0m, Pxx0);

        vector<VectorXd> Xtru;
        double time;
        do {
            cout << "       Trying" << endl;
            time = 0;
            VectorXd X(6);
            X = X0m;
            for (int i = 0; i < 3; i++) {
                X(i)   += genx0();
                X(i+3) += genv0();
            }
            Xtru.clear();
            Xtru.push_back(X);
            do {
                Vector3d fd;
                fd << genw(), genw(), genw();
                X = f(time, time+dt, X, fd);
                time += dt;
                Xtru.push_back(X);
            } while (X(2) > 0);
         } while (time < tmin);

        int steps = Xtru.size();

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

        MatrixXd Ztru(2, steps);
        for (int k = 0; k < steps; k++) {
            Ztru.col(k) = h(t(k), Xtru[k]);
            Ztru(0,k) += genn();
            Ztru(1,k) += genn();
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

        cout << "   CUT-6" << endl;
        timer.tick();
        cut6.run(t, Ztru);
        run_times(j-1, 3) = timer.tock();
        cut6.save(filter_file(gauss, "cut6", j));

    }

    vector<string> filters;
    filters.push_back("house");
    filters.push_back("ukf");
    filters.push_back("cut4");
    filters.push_back("cut6");

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
