#include "ukf.hpp"
#include "house.hpp"
#include "eigen_csv.hpp"
#include "timer.hpp"
#include "pearsonator.hpp"

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

// One degree (in radians)
const double deg = M_PI / 180;

// Dynamical model
VectorXd ct_prop(double ti, double tf, const VectorXd& x, const VectorXd& w) {
    VectorXd xf(5);
    double W = x(4);
    double T = tf - ti;
    /*
    if (fabs(W) < 1E-9) {
        xf = x;
        xf(0) += T * xf(1);
        xf(2) += T * xf(3);
    } else {*/
        double sinWT = sin(W*T);
        double cosWT = cos(W*T);
        MatrixXd F(5,5);
        F << 1, sinWT/W,     0, -(1-cosWT)/W, 0,
             0, cosWT,       0, -sinWT,       0,
             0, (1-cosWT)/W, 1, sinWT/W,      0,
             0, sinWT,       0, cosWT,        0,
             0, 0,           0, 0,            1;
        xf = F * x;
    //}
    return xf + w;
}

// Measurement model
Vector2d ct_meas(double t, const VectorXd& x) {

    Vector2d z;

    double xi  = x(0);
    double eta = x(2);

    z << sqrt(xi*xi + eta*eta), atan2(eta, xi);

    return z;

}

// True state as function of time
VectorXd ct_tru(double t) {

    const double v  = 120;
    const double w1 =  1 * deg;
    const double w2 = -3 * deg;
    const double r1 = fabs(v / w1);
    const double r2 = fabs(v / w2);
    const double t1 = 90;
    const double t2 = 30;
    const double tl = 125;
    const double l  = v * tl;
    const double x0 = 25000;
    const double y0 = 10000;

    VectorXd x(5);
    double dt, xc, yc;

    if (t <= tl) {
        x(0) = x0 - v * t;
        x(2) = y0;
        x(1) = -v;
        x(3) = 0;
        x(4) = 0;
    } else if (t <= tl + t1) {
        dt = t - tl;
        xc = x0 - l;
        yc = y0 - r1;
        x(0) =  xc + r1 * cos(w1 * dt + M_PI/2);
        x(2) =  yc + r1 * sin(w1 * dt + M_PI/2);
        x(1) = -w1 * r1 * sin(w1 * dt + M_PI/2);
        x(3) =  w1 * r1 * cos(w1 * dt + M_PI/2);
        x(4) =  w1;
    } else if (t <= 2*tl + t1) {
        dt = t - tl - t1;
        x(0) = x0 - l - r1;
        x(2) = y0 - r1 - v * dt;
        x(1) = 0;
        x(3) = -v;
        x(4) = 0;
    } else if (t <= 2*tl + t1 + t2) {
        dt = t - 2*tl - t1;
        xc = x0 - l - r1 - r2;
        yc = y0 - l - r1;
        x(0) =  xc + r2 * cos(w2 * dt);
        x(2) =  yc + r2 * sin(w2 * dt);
        x(1) = -w2 * r2 * sin(w2 * dt);
        x(3) =  w2 * r2 * cos(w2 * dt);
        x(4) =  w2;
    } else {
        dt = t - 2*tl - t1 - t2;
        x(0) = x0 - l - r1 - r2 - v*dt;
        x(2) = y0     - r1 - r2 - l;
        x(1) = -v;
        x(3) = 0;
        x(4) = 0;
    }

    return x;

}

// Generate filename for state estimates
std::string est_filename(const std::string& filter, int i, int j) {
    std::string f = "out/ct_est_";
    f += filter;
    f += "_";
    f += std::to_string(i+1);
    f += "_";
    f += std::to_string(j+1);
    f += ".csv";
    return f;
}

int main() {

    // CUT directory
    UKF::cut_dir = "../CUT/";

    // UKF state & measurement models
    UKF::dyn_model f = ct_prop;
    UKF::meas_model h = ct_meas;

    // HOUSE measurement model
    HOUSE::meas_model hh = [] (double t, const VectorXd& x, const VectorXd& n)
        -> VectorXd {
            return ct_meas(t, x) + n;
        };

    // Constants
    double sr, sth, l1, l2;
    sr  = 100;
    sth = 1 * deg;
    l1  = 0.16;
    l2  = 0.01;

    // Skewness & kurtosis
    double skewnth, skewnr, kurtnth, kurtnr, kurt0, kurtw;
    skewnth = -1;
    kurtnth =  20;
    skewnr  = -1;
    kurtnr  = 20;
    kurt0 = 10;
    kurtw = 10;

    // Measurement noise covariance
    Matrix2d R;
    R << sr*sr, 0,
         0,     sth*sth;

    // Prior mean & covariance
    VectorXd mxi(5), cxx(5);
    cxx << 1000, 10, 1000, 10, 1*deg;
    mxi << 25000, -120, 10000, 0, 0.000001;
    MatrixXd Pxxi = cxx.cwiseAbs2().asDiagonal();

    // HOUSE distributions
    HOUSE::Dist distXi(Pxxi), distn(R);
    distXi.mean = mxi;
    distn.skew << skewnr, skewnth;
    distn.kurt << kurtnr, kurtnth;
    distXi.kurt.setConstant(kurt0);

    // Pearson noise generator
    mt19937_64 mt;
    Pearsonator::TypeIV gennr(0, sr, skewnr, kurtnr),
        gennth(0, sth, skewnth, kurtnth);

    // Time step lengths
    const int Ntimes = 4;
    VectorXd Ts(Ntimes);
    Ts << 1, 2.5, 3.9, 5;
    //Ts.setLinSpaced(Ntimes, 1, 5);

    // Save time step lengths
    EigenCSV::write(Ts, "out/times.csv");

    // Number of tests
    const int Ntest = 100;

    // Timer
    Timer timer;

    // Average runtime table
    MatrixXd runtime_table(Ntimes, 6);
    runtime_table.col(0) = Ts;

    // Run tests for various time step lengths
    for (int i = 0; i < Ntimes; i++) {

        // Time step
        double T = Ts(i);

        //cout << "-------------------------------------------" << endl;
        cout << "Running T = " << T << " s" << endl;
        //cout << "-------------------------------------------" << endl;

        // Times
        double Ttot = 125 + 90 + 125 + 30 + 125;
        int steps = (int) ceil(Ttot / T);
        VectorXd t(steps);
        t.setLinSpaced(steps, 0, T*(steps-1));

        // True states
        vector<VectorXd> Xtru;
        for (int k = 0; k < steps; k++)
            Xtru.push_back(ct_tru(t(k)));

        // Save true states
        MatrixXd xtru_table(steps, 6);
        xtru_table.col(0) = t;
        for (int k = 0; k < steps; k++)
            xtru_table.row(k).tail(5) = Xtru[k];
        string filename = "out/ct_true_";
        filename += to_string(i+1);
        filename += ".csv";
        vector<string> header(6);
        header[0] = "Time (s)";
        header[1] = "x (m)";
        header[2] = "xd (m/s)";
        header[3] = "y (m)";
        header[4] = "yd (m/s)";
        header[5] = "Omega (rad/s)";
        EigenCSV::write(xtru_table, header, filename);

        // Process noise covariance
        double t22 = T * T / 2;
        double t33 = T * T * T / 3;
        MatrixXd Q(5,5);
        Q << t33, t22, 0,   0,   0,
             t22, T,   0,   0,   0,
             0,   0,   t33, t22, 0,
             0,   0,   t22, T,   0,
             0,   0,   0,   0,   (l2/l1)*T;
        Q *= l1;

        // HOUSE process noise distributions
        HOUSE::Dist distw(Q);
        distw.kurt.setConstant(kurtw);

        // Initialize UKF & CUT filters
        UKF ukf (f, h, true, 0, mxi, Pxxi, Q, R, UKF::sig_type::JU,   1);
        UKF cut4(f, h, true, 0, mxi, Pxxi, Q, R, UKF::sig_type::CUT4, 1);
        UKF cut6(f, h, true, 0, mxi, Pxxi, Q, R, UKF::sig_type::CUT6, 1);
        UKF cut8(f, h, true, 0, mxi, Pxxi, Q, R, UKF::sig_type::CUT8, 1);

        // Initialize HOUSE
        HOUSE house(f, hh, 2, 0, distXi, distw, distn, 0.1);

        // Runtimes
        MatrixXd runtime(5, Ntest);

        // Run multiple tests
        for (int j = 0; j < Ntest; j++) {

            // cout << "Running Trial " << j+1 << endl;

            // Reset filters
            ukf.reset(0, mxi, Pxxi);
            cut4.reset(0, mxi, Pxxi);
            cut6.reset(0, mxi, Pxxi);
            cut8.reset(0, mxi, Pxxi);
            house.reset(0, distXi);

            // Measurements
            MatrixXd Ztru(2,steps);
            for (int k = 0; k < steps; k++) {
                Ztru.col(k) = ct_meas(t(k), Xtru[k]);
                Ztru(0,k) += gennr (mt);
                Ztru(1,k) += gennth(mt);
            }

            // Run UKF
            // cout << "    UKF" << endl;
            timer.tick();
            ukf.run(t, Ztru);
            runtime(0, j) = timer.tock();
            ukf.save(est_filename("ukf", i, j));

            // Run CUT-4
            //cout << "    CUT-4" << endl;
            timer.tick();
            cut4.run(t, Ztru);
            runtime(1, j) = timer.tock();
            cut4.save(est_filename("cut4", i, j));

            // Run CUT-6
            //cout << "    CUT-6" << endl;
            timer.tick();
            cut6.run(t, Ztru);
            runtime(2, j) = timer.tock();
            cut6.save(est_filename("cut6", i, j));

            // Run CUT-8
            //cout << "    CUT-8" << endl;
            timer.tick();
            cut8.run(t, Ztru);
            runtime(3, j) = timer.tock();
            cut8.save(est_filename("cut8", i, j));

            // Run HOUSE
            //cout << "    HOUSE" << endl;
            timer.tick();
            house.run(t, Ztru);
            runtime(4, j) = timer.tock();
            house.save(est_filename("house", i, j));

        }

        // Average runtimes
        runtime_table.row(i).tail(5) = runtime.rowwise().mean();

    }

    // Save average runtimes
    vector<string> rth(6);
    rth[0] = "STEP";
    rth[1] = "UKF";
    rth[2] = "CUT-4";
    rth[3] = "CUT-6";
    rth[4] = "CUT-8";
    rth[5] = "HOUSE";
    EigenCSV::write(runtime_table, rth, "out/runtimes.csv");

    return 0;

}
