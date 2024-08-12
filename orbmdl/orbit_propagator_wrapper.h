#ifndef ORBIT_PROPAGATOR_WRAPPER_H
#define ORBIT_PROPAGATOR_WRAPPER_H

#include "testOrbDet.hpp"
using namespace std;

class OrbitPropagatorWapper
{
public:
    DynamicModel::stf accelerationModel;
    // OrbitPropagatorWapper(const string &configFilename, double param1, int param2);
    OrbitPropagatorWapper(const string &configFilename);
    vector<double> propagate();

private:
    void readConfigFile(const string &fileName);
    void initGlobalVariables();
    VectorXd stdVec2EigenVec(const vector<double> &stdVec);
    bool initEGMCoef(const string &filename);

    // Member variables
    double param1_;
    int param2_;
    struct ForceModels forceModelsPropOpt_;
    struct ScenarioInfo snrInfo_;
    struct MeasModel measMdl_;
    struct Filters filters_;
    struct InitialState initialState_;
    struct FileInfo suppFiles_;
    double leapSec_;
    erp_t erpt_;
    IERS iersInstance_;
    EGMCoef egm_;
    void *pJPLEph_;
    struct EpochInfo epoch_;
};

#endif // ORBIT_PROPAGATOR_H