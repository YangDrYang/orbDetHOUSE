#ifndef ORBIT_PROPAGATOR_WRAPPER_H
#define ORBIT_PROPAGATOR_WRAPPER_H

#include "testOrbDet.hpp"
using namespace std;

class OrbitPropagatorWrapper
{
public:
    // Use a lambda function to capture the this pointer and call the member function
    OrbitPropagatorWrapper(const string &configFilename);
    MatrixXd propagateOrbit();
    void saveResults(const MatrixXd &results, const vector<string> &header, const string &filePath);

    void readConfigFile(const string &fileName);
    void initGlobalVariables();
    VectorXd stdVec2EigenVec(const vector<double> &stdVec);
    bool initEGMCoef(const string &filename);

private:
    VectorXd accelerationModel(double t, const VectorXd &x, const VectorXd &w);

private:
    // Member variables
    double param1_;
    int param2_;
    DynamicModel orbFun;
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
    Propagator orbitProp;
};

#endif // ORBIT_PROPAGATOR_H