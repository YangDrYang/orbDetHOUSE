// library headers
#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>

// filter headers
#include "pearsonator.hpp"
#include "srhouse.hpp"
#include "srukf.hpp"
#include "house.hpp"
#include "ukf.hpp"
#include "ut.hpp"
#include "dyn.hpp"
#include "eigen_csv.hpp"
#include "timer.hpp"

#define DEFAULT_CONFIG_FILENAME "yamls/config_lorenz.yml"

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
    // cout << "fd:\t" << fd.transpose();
    k += fd;

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
struct FilterOpt
{
    bool srhouse;
    bool house;
    bool ukf;
    bool cut4;
    bool cut6;
    double delta;
    double weight;
    int numTrials;
};
struct SimulationInfo
{
    double tStep;
    int nTotalSteps;
    string outDir;
};
struct InitialState
{
    int dimState;
    double fValue;
    double stdInitialNoise;
    double skewInitialNoise;
    double kurtInitialNoise;
};
struct ProcessNoise
{
    double stdProcessNoise;
    double skewProcessNoise;
    double kurtProcessNoise;
};
struct MeasurementNoise
{
    int dimMeasurement;
    string typeMeasurement;
    double stdMeasurementNoise;
    double skewMeasurementNoise;
    double kurtMeasurementNoise;
};
void readConfigFile(string fileName, FilterOpt &optFilter, SimulationInfo &snrInfo,
                    InitialState &iniState, ProcessNoise &proNoise, MeasurementNoise &measNoise)
{
    // load file
    YAML::Node config = YAML::LoadFile(fileName);

    // read filter options (required)
    YAML::Node filterOpts = config["filter_options"];
    optFilter.srhouse = filterOpts["SRHOUSE"].as<bool>();
    optFilter.house = filterOpts["HOUSE"].as<bool>();
    optFilter.ukf = filterOpts["UKF"].as<bool>();
    optFilter.cut4 = filterOpts["CUT4"].as<bool>();
    optFilter.cut6 = filterOpts["CUT6"].as<bool>();
    optFilter.numTrials = filterOpts["num_trials"].as<int>();
    optFilter.delta = filterOpts["delta"].as<double>();
    optFilter.weight = filterOpts["weight"].as<double>();

    // read simulation parameters (required)
    YAML::Node simParams = config["simulation_parameters"];
    snrInfo.tStep = simParams["time_step"].as<double>();
    snrInfo.nTotalSteps = simParams["total_step"].as<int>();
    snrInfo.outDir = simParams["output_directory"].as<string>();

    // read initial state (required)
    YAML::Node initialState = config["initial_state"];
    iniState.dimState = initialState["state_dim"].as<int>();
    iniState.fValue = initialState["fvalue"].as<double>();
    iniState.stdInitialNoise = initialState["initial_noise_std"].as<double>();
    iniState.skewInitialNoise = initialState["initial_noise_skew"].as<double>();
    iniState.kurtInitialNoise = initialState["initial_noise_kurt"].as<double>();

    // read process noise (required)
    YAML::Node processNoise = config["process_noise"];
    proNoise.stdProcessNoise = processNoise["process_noise_std"].as<double>();
    proNoise.skewProcessNoise = processNoise["process_noise_skew"].as<double>();
    proNoise.kurtProcessNoise = processNoise["process_noise_kurt"].as<double>();

    // read measurement noise (required)
    YAML::Node measurementNoise = config["measurement_noise"];
    measNoise.dimMeasurement = measurementNoise["measurement_dim"].as<int>();
    measNoise.typeMeasurement = measurementNoise["measurement_type"].as<string>();
    measNoise.stdMeasurementNoise = measurementNoise["measurement_std"].as<double>();
    measNoise.skewMeasurementNoise = measurementNoise["measurement_skew"].as<double>();
    measNoise.kurtMeasurementNoise = measurementNoise["measurement_kurt"].as<double>();
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

    struct FilterOpt optFilter;
    struct SimulationInfo snrInfo;
    struct InitialState iniState;
    struct ProcessNoise proNoise;
    struct MeasurementNoise measNoise;
    readConfigFile(configFilename, optFilter, snrInfo, iniState, proNoise, measNoise);

    bool gauss = false;
    if (measNoise.typeMeasurement == "gauss")
    {
        gauss = true;
    }

    int numTrials = optFilter.numTrials;

    double stdx0, stdw, stdn;
    stdx0 = iniState.stdInitialNoise;
    stdw = proNoise.stdProcessNoise;
    stdn = measNoise.stdMeasurementNoise;

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

        skeww = proNoise.skewProcessNoise;
        kurtw = proNoise.kurtProcessNoise;

        skewn = measNoise.skewMeasurementNoise;
        kurtn = measNoise.kurtMeasurementNoise;

        skew0 = iniState.skewInitialNoise;
        kurt0 = iniState.kurtInitialNoise;
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

    const int dimState = iniState.dimState;
    VectorXd initialStateVec(dimState);
    initialStateVec.setConstant(iniState.fValue);

    double absErr = 1.0e-8;
    double relErr = 1.0e-8;
    DynamicModel::stf accMdl = lorenz96Model;
    DynamicModel stateFun(accMdl, dimState, absErr, relErr);

    double time = 0, dt = snrInfo.tStep;
    int nTotalSteps = snrInfo.nTotalSteps;
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
        VectorXd fd(dimState);
        for (int j = 0; j < dimState; j++)
            fd(j) = genw();
        // propStateVec = rk4(propStateVec, dt, 8.0);
        // propStateVec = stateFun(time, time + dt, propStateVec, VectorXd::Zero(dimState));
        propStateVec = stateFun(time, time + dt, propStateVec, fd);
        time += dt;
        tableTrajTruth.row(k).tail(dimState) = propStateVec;
        cout << time << " " << propStateVec.transpose() << endl;

        // cout << "The " << k + 1 << "th time step" << endl;
    }
    cout << "The total time consumption is:\t" << timer.tock() << endl;
    // header for the saved file

    vector<string> headerTraj(dimState + 1);
    headerTraj[0] = "tSec";
    for (int k = 1; k <= dimState; k++)
    {
        headerTraj[k] = "x";
        headerTraj[k] += to_string(k);
    }
    string propFile = snrInfo.outDir;
    propFile += (gauss ? "_gauss/trajectory_truth.csv" : "_pearson/trajectory_truth.csv");
    EigenCSV::write(tableTrajTruth, headerTraj, propFile);

    const int dimMeas = measNoise.dimMeasurement;
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

    VectorXd x0stdVec(dimState), wstdVec(dimState);
    x0stdVec.setConstant(stdx0);
    MatrixXd Pxx0 = x0stdVec.array().square().matrix().asDiagonal();
    wstdVec.setConstant(stdw);
    MatrixXd Pww = wstdVec.array().square().matrix().asDiagonal();
    MatrixXd Pnn(1, 1);
    Pnn(0, 0) = stdn * stdn;

    MatrixXd Ztru(dimMeas, nTotalSteps);
    for (int k = 0; k < nTotalSteps; k++)
    {
        const VectorXd stateVec = tableTrajTruth.row(k).tail(dimState).transpose();
        Ztru.col(k) = h(tSec(k), stateVec);
        Ztru(0, k) += genn();
    }

    UKF ukf(stateFun, h, true, 0, 1e6, initialStateVec, Pxx0, Pww, Pnn, UKF::sig_type::JU, 1);
    UKF cut4(stateFun, h, true, 0, 1e6, initialStateVec, Pxx0, Pww, Pnn, UKF::sig_type::CUT4, 1);
    UKF cut6(stateFun, h, true, 0, 1e6, initialStateVec, Pxx0, Pww, Pnn, UKF::sig_type::CUT6, 1);

    // HOUSE distributions for state
    Dist distXi(Pxx0);
    distXi.mean = initialStateVec;
    // HOUSE distributions for state noise
    Dist distw(Pww);
    // HOUSE distributions for measurement noise
    Dist distn(Pnn);

    HOUSE house(stateFun, hh, dimMeas, 0, 1e6, distXi, distw, distn, optFilter.delta);
    SRHOUSE srhouse(stateFun, hh, dimMeas, 0, 1e6, distXi, distw, distn, optFilter.weight);

    MatrixXd run_times(numTrials, 5);
    for (int j = 1; j <= numTrials; j++)
    {
        cout << "Trial " << j << endl;

        for (int k = 0; k < dimState; k++)
            propStateVec(k) = initialStateVec(k) + genx0();
        // cout << "propStateVec" << propStateVec.transpose() << endl;
        distXi.mean = propStateVec;

        if (optFilter.srhouse)
        {
            cout << "   SRHOUSE" << endl;
            srhouse.reset(0, distXi);
            timer.tick();
            srhouse.run(tSec, Ztru);
            run_times(j - 1, 0) = timer.tock();
            srhouse.save(filter_file(gauss, "srhouse", j));
        }

        if (optFilter.house)
        {
            cout << "   HOUSE" << endl;
            house.reset(0, distXi);
            timer.tick();
            house.run(tSec, Ztru);
            run_times(j - 1, 1) = timer.tock();
            house.save(filter_file(gauss, "house", j));
        }

        if (optFilter.ukf)
        {
            cout << "   UKF" << endl;
            ukf.reset(0, propStateVec, Pxx0);
            timer.tick();
            ukf.run(tSec, Ztru);
            run_times(j - 1, 2) = timer.tock();
            ukf.save(filter_file(gauss, "ukf", j));
        }

        if (optFilter.cut4)
        {
            cout << "   CUT4" << endl;
            cut4.reset(0, propStateVec, Pxx0);
            timer.tick();
            cut4.run(tSec, Ztru);
            run_times(j - 1, 3) = timer.tock();
            cut4.save(filter_file(gauss, "cut4", j));
        }

        if (optFilter.cut6)
        {
            cout << "   CUT6" << endl;
            cut6.reset(0, propStateVec, Pxx0);
            timer.tick();
            cut6.run(tSec, Ztru);
            run_times(j - 1, 4) = timer.tock();
            cut6.save(filter_file(gauss, "cut6", j));
        }
    }
}
