// library headers
#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>

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

#define DEFAULT_CONFIG_FILENAME "config.yaml"
#define JPL_EPHEMERIS_FILENAME "unxp2000.405"



struct EpochInfo {
    double startMJD, endMJD;    // MJD in days
    double timeStep;            // time step in seconds
};

struct Errors {
    double elevationErr, azimuthErr, rangeErr, rangeRateErr;
};

struct Filters {
    bool house;
    bool ukf;
    bool cut4;
    bool cut6;    
};

// prototypes
VectorXd orbitModel(double t, const VectorXd& X);
void initEGMCoef(string filename);
void initGlobalVariables(VectorXd& initialState, string stateType);
Eigen::VectorXd stdVec2EigenVec(const std::vector<double>& stdVec);
void readConfigFile(string fileName, ForceModels& options, struct EpochInfo& epoch, Eigen::VectorXd& initialState,
                    Eigen::VectorXd& groundStation, struct Filters& filters, int& numTrials, struct Errors& errorStd,
                    string& iniatialStateType);
Eigen::MatrixXd generateTrueResults(DynamicModel& f, struct EpochInfo epoch, VectorXd initialState);