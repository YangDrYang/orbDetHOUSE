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



// prototypes
VectorXd orbitModel(double t, const VectorXd& X);
void initEGMCoef(string filename);
void initGlobalVariables(VectorXd rvECI);
