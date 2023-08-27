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
#include "configDefault.hpp"
#include "auxillaryData.hpp"
#include "forceModels.hpp"
#include "coordTrans.hpp"
#include "satRefSys.hpp"
#include "jplEph.hpp"
#include "jpl_eph.hpp"

// filter headers
#include "srhouse.hpp"
#include "srukf.hpp"
#include "house.hpp"
#include "ukf.hpp"
#include "ut.hpp"
#include "dyn.hpp"
#include "eigen_csv.hpp"
#include "timer.hpp"
// #inlcude "multiNormal.hpp"

using namespace Eigen;
using namespace std;

#define DEFAULT_CONFIG_FILENAME "config.yaml"

struct EpochInfo
{
    double startMJD, endMJD; // MJD in days
    double timeStep;         // time step in seconds
    double maxTimeStep;      // max time step in seconds
    double timePass;         // time duration in each pass in seconds
};
struct FileInfo
{
    string grvFile;
    string ephFile;
    string erpFile;
};
struct ScenarioInfo
{
    EpochInfo epoch;
    string outDir;
};

struct Filters
{
    bool squareRoot;
    bool house;
    bool ukf;
    bool cut4;
    bool cut6;
    bool initNoise;
    int numTrials;
};

struct InitialState
{
    int dimState;
    string initialStateType;
    VectorXd initialStateVec;
    MatrixXd initialCovarianceMat;
    MatrixXd processNoiseCovarianceMat;
    double initialSkewness;
    double initialKurtosis;
};

struct Errors
{
    double rightAscensionErr, declinationErr;
    double elevationErr, azimuthErr, rangeErr, rangeRateErr;
    VectorXd stdVec, skewVec, kurtVec;
};

struct MeasModel
{
    string measFile;
    VectorXd groundStation;
    int dimMeas;
    Errors errorStatistics;
};

// // prototypes
// VectorXd orbitModel(double t, const VectorXd &X);
void initEGMCoef(string filename);
VectorXd stdVec2EigenVec(const vector<double> &stdVec);
void readConfigFile(string fileName, ForceModels &optFilter, struct ScenarioInfo &snrInfo, struct InitialState &initialState,
                    struct MeasModel &measMdl, struct Filters &filters, struct FileInfo &suppFiles);
