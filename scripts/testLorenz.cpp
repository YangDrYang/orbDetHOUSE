#include "testOrbDet.hpp"
#include "pearsonator.hpp"

#define R_EARTH 6371E3
#define DEG M_PI / 180
#define ARC_MIN M_PI / (180 * 60)
#define ARC_SEC M_PI / (180 * 60 * 60)
#define MAXLEAPS 18

string filter_file(bool gauss, const string &filter, int trial)
{
    string filename = "out/out_lorenz";
    filename += (gauss ? "_gauss/" : "_pearson/");
    filename += filter;
    filename += "_est_";
    filename += (gauss ? "gauss_" : "pearson_");
    filename += to_string(trial);
    filename += ".csv";
    return filename;
}

// VectorXd lorenz96Model(double t, const VectorXd &X, const VectorXd &fd)
// {
//     // Setting up vector
//     int dimState = X.size();
//     VectorXd dX(dimState);
//     const double F = 8.0;

//     // % first the 3 problematic cases : 1, 2, J
//     // k(1)=(X(2)-X(J-1))*X(J)-X(1);
//     // k(2)=(X(3)-X(J))*X(1)-X(2);
//     // k(J)=(X(1)-X(J-2))*X(J-1)-X(J);
//     // %then the general case
//     // for j=3:J-1
//     //  k(j)=(X(j+1)-X(j-2)).*X(j-1)-X(j);
//     // end
//     // %add the F
//     // k=k+F;
//     //     dX(0)=(X(1)-X(dimState-1))*X(J)-X(1);

//     // Loops over indices (with operations and C++ underflow indexing handling edge cases)
//     for (int i = 0; i < dimState; ++i)
//     {
//         dX[i] = ((X[(i + 1) % dimState] - X[(i - 2 + dimState) % dimState]) * X[(i - 1 + dimState) % dimState]) - X[i] + F;
//     }

//     return dX;
// }

VectorXd lorenz96Model(double t, const VectorXd &state, const VectorXd &fd)
{
    int J = state.size();
    VectorXd k(J);
    double F = 8.0;

    k(0) = (state(1) - state(J - 2)) * state(J - 1) - state(0);
    k(1) = (state(2) - state(J - 1)) * state(0) - state(1);
    k(J - 1) = (state(0) - state(J - 3)) * state(J - 2) - state(J - 1);

    for (int j = 2; j < J - 1; ++j)
    {
        k(j) = (state(j + 1) - state(j - 2)) * state(j - 1) - state(j);
    }

    return k.array() + F;
}

// VectorXd lorenz96Model(const VectorXd &state, double F)
// {
//     int J = state.size();
//     VectorXd k(J);

//     k(0) = (state(1) - state(J - 2)) * state(J - 1) - state(0);
//     k(1) = (state(2) - state(J - 1)) * state(0) - state(1);
//     k(J - 1) = (state(0) - state(J - 3)) * state(J - 2) - state(J - 1);

//     for (int j = 2; j < J - 1; ++j)
//     {
//         k(j) = (state(j + 1) - state(j - 2)) * state(j - 1) - state(j);
//     }
//     // cout << k.transpose().array() << endl;
//     // cout << k.array() + F << endl;
//     return k.array() + F;
// }

// VectorXd rk4(VectorXd state, double dt, double F)
// {
//     int J = state.size();
//     VectorXd k1(J), k2(J), k3(J), k4(J);

//     k1 = lorenz96Model(state, F);
//     k2 = lorenz96Model(state + 0.5 * dt * k1, F);
//     k3 = lorenz96Model(state + 0.5 * dt * k2, F);
//     k4 = lorenz96Model(state + dt * k3, F);

//     return state + (1.0 / 6.0) * dt * (k1 + 2 * k2 + 2 * k3 + k4);
// }

int main(int argc, char *argv[])
{
    bool gauss = true;
    // bool gauss = false;
    int numTrials = 100;

    const int dimState = 5;
    VectorXd initialStateVec(dimState);
    initialStateVec << 8, 8, 8, 8, 8;
    MatrixXd initialCov(dimState, dimState);

    double absErr = 1.0e-8;
    double relErr = 1.0e-8;
    DynamicModel::stf accMdl = lorenz96Model;
    DynamicModel orbFun(accMdl, dimState, absErr, relErr);

    double time = 0, dt = 0.5;
    int nTotalSteps = 100;
    cout << "total steps:\t" << nTotalSteps << endl;
    // linear spaced times
    VectorXd tSec;
    tSec.setLinSpaced(nTotalSteps, 0, (nTotalSteps - 1) * dt);
    cout << "tSec\t" << tSec << endl;
    MatrixXd tableTrajTruth(nTotalSteps, dimState + 1);
    tableTrajTruth.col(0) = tSec;
    tableTrajTruth.row(0).tail(dimState) = initialStateVec;

    VectorXd propStateVec = initialStateVec;
    Timer timer;
    timer.tick();
    for (int k = 1; k < nTotalSteps; k++)
    {
        // propStateVec = rk4(propStateVec, dt, 8.0);
        propStateVec = orbFun(time, time + dt, propStateVec, VectorXd::Zero(dimState));
        time += dt;
        tableTrajTruth.row(k).tail(dimState) = propStateVec;
        cout << time << " " << propStateVec.transpose() << endl;

        // cout << "The " << k + 1 << "th time step" << endl;
    }
    cout << "The total time consumption is:\t" << timer.tock() << endl;
    // header for the saved file
    vector<string> headerTraj({"tSec", "x1", "x2", "x3", "x4", "x5"});
    string propFile = "out/out_lorenz";
    propFile += (gauss ? "_gauss/trajectory_truth.csv" : "_pearson/trajectory_truth.csv");
    EigenCSV::write(tableTrajTruth, headerTraj, propFile);

    const int dimMeas = 1;
    UKF::meas_model h = [](double t, const VectorXd &x)
        -> VectorXd
    {
        return x.head(1);
    };
    HOUSE::meas_model hh = [](double t, const VectorXd &x, const VectorXd &n)
        -> VectorXd
    {
        return x.head(1) + n;
    };

    double stdx0, stdw, stdn;
    stdx0 = 0.5;
    stdw = 1e-6;
    stdn = 0.2;

    double skeww, skewn, kurtw, kurtn, skew0, kurt0;

    if (gauss)
    {

        skeww = 0;
        kurtw = 3;

        skewn = 0;
        kurtn = 3;

        skew0 = 0;
        kurt0 = 3;
    }
    else
    {

        skeww = 1;
        kurtw = 30;

        skewn = -1;
        kurtn = 30;

        skew0 = 1;
        kurt0 = 30;
    }

    Pearsonator::TypeIV genw_p(0, stdw, skeww, kurtw),
        genn_p(0, stdn, skewn, kurtn),
        genx0_p(0, stdx0, skew0, kurt0);

    normal_distribution<double> genw_g(0, stdw),
        genn_g(0, stdn),
        genx0_g(0, stdx0);

    mt19937_64 mt(0);

    typedef function<double(void)> noisemaker;

    noisemaker genw = [&]() -> double
    { return gauss ? genw_g(mt) : genw_p(mt); };
    noisemaker genn = [&]() -> double
    { return gauss ? genn_g(mt) : genn_p(mt); };
    noisemaker genx0 = [&]() -> double
    { return gauss ? genx0_g(mt) : genx0_p(mt); };

    VectorXd x0stdVec(dimState), wstdVec(dimState);
    x0stdVec.setConstant(stdx0);
    MatrixXd Pxx0 = x0stdVec.array().square().matrix().asDiagonal();
    wstdVec.setConstant(stdw);
    MatrixXd Pww = wstdVec.array().square().matrix().asDiagonal();
    MatrixXd Pnn(1, 1);
    Pnn(0, 0) = stdn * stdn;

    MatrixXd Ztru(1, nTotalSteps);
    for (int k = 0; k < nTotalSteps; k++)
    {
        const VectorXd stateVec = tableTrajTruth.row(k).tail(dimState).transpose();
        Ztru.col(k) = h(tSec(k), stateVec);
        Ztru(0, k) += genn();
    }

    UKF ukf(orbFun, h, true, 0, 1e6, initialStateVec, Pxx0, Pww, Pnn, UKF::sig_type::JU, 1);
    UKF cut4(orbFun, h, true, 0, 1e6, initialStateVec, Pxx0, Pww, Pnn, UKF::sig_type::CUT4, 1);
    UKF cut6(orbFun, h, true, 0, 1e6, initialStateVec, Pxx0, Pww, Pnn, UKF::sig_type::CUT6, 1);

    // HOUSE distributions for state
    Dist distXi(Pxx0);
    distXi.mean = initialStateVec;
    // HOUSE distributions for state noise
    Dist distw(Pww);
    // HOUSE distributions for measurement noise
    Dist distn(Pnn);
    HOUSE house(orbFun, hh, dimMeas, 0, 1e6, distXi, distw, distn, 0);
    SRHOUSE srhouse(orbFun, hh, dimMeas, 0, 1e6, distXi, distw, distn, -0.5);

    MatrixXd run_times(numTrials, 5);
    for (int j = 1; j <= numTrials; j++)
    {
        cout << "Trial " << j << endl;

        for (int k = 0; k < dimState; k++)
            propStateVec(k) = initialStateVec(k) + genx0();
        cout << "propStateVec" << propStateVec.transpose() << endl;
        distXi.mean = propStateVec;

        house.reset(0, distXi);
        srhouse.reset(0, distXi);
        ukf.reset(0, propStateVec, Pxx0);
        cut4.reset(0, propStateVec, Pxx0);
        cut6.reset(0, propStateVec, Pxx0);

        cout << "   SRHOUSE" << endl;
        timer.tick();
        srhouse.run(tSec, Ztru);
        run_times(j - 1, 0) = timer.tock();
        srhouse.save(filter_file(gauss, "srhouse", j));

        cout << "   HOUSE" << endl;
        timer.tick();
        house.run(tSec, Ztru);
        run_times(j - 1, 1) = timer.tock();
        house.save(filter_file(gauss, "house", j));

        cout << "   UKF" << endl;
        timer.tick();
        ukf.run(tSec, Ztru);
        run_times(j - 1, 2) = timer.tock();
        ukf.save(filter_file(gauss, "ukf", j));

        cout << "   CUT4" << endl;
        timer.tick();
        cut4.run(tSec, Ztru);
        run_times(j - 1, 3) = timer.tock();
        cut4.save(filter_file(gauss, "cut4", j));

        cout << "   CUT6" << endl;
        timer.tick();
        cut6.run(tSec, Ztru);
        run_times(j - 1, 4) = timer.tock();
        cut6.save(filter_file(gauss, "cut6", j));
    }
}
