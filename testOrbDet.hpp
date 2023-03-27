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

// model headers
#include "auxillaryData.hpp"
#include "forceModels.hpp"
#include "coordTrans.hpp"
#include "satRefSys.hpp"
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
// #inlcude "multiNormal.hpp"

using namespace Eigen;
using namespace std;

#define DEFAULT_CONFIG_FILENAME "config.yaml"
#define JPL_EPHEMERIS_FILENAME "./auxdata/unxp2000.405"

struct EpochInfo
{
    double startMJD, endMJD; // MJD in days
    double timeStep;         // time step in seconds
    double timePass;         // time duration in each pass in seconds
};
struct FileInfo
{
    string outDir;
};
struct SimInfo
{
    EpochInfo epoch;
    FileInfo file;
};

struct Filters
{
    bool house;
    bool ukf;
    bool cut4;
    bool cut6;
    int numTrials;
};

struct InitialState
{
    int dimState;
    string initialStateType;
    VectorXd initialStateVec;
    MatrixXd initialCovarianceMat;
    MatrixXd processNoiseCovarianceMat;
};

struct Errors
{
    double elevationErr, azimuthErr, rangeErr, rangeRateErr;
};

struct MeasModel
{
    VectorXd groundStation;
    int dimMeas;
    Errors errorStd;
};

// prototypes
VectorXd
orbitModel(double t, const VectorXd &X);
void initEGMCoef(string filename);
void initGlobalVariables(VectorXd &initialState, string stateType);
VectorXd stdVec2EigenVec(const vector<double> &stdVec);
void readConfigFile(string fileName, ForceModels &optTruth, ForceModels &optFilter, struct SimInfo &simInfo, struct InitialState &initialState,
                    MeasModel &measMdl, struct Filters &filters);
MatrixXd generateTrueResults(DynamicModel &f, struct EpochInfo epoch, VectorXd initialState);