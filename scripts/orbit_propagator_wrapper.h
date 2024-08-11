#ifndef ORBIT_PROPAGATOR_WRAPPER_H
#define ORBIT_PROPAGATOR_WRAPPER_H

#include "testOrbDet.hpp"

class OrbitPropagatorWapper
{
public:
    DynamicModel::stf accelerationModel;
    OrbitPropagatorWapper(const std::string &configFilename, double param1, int param2);
    std::vector<double> propagate(double initial_position, double initial_velocity, double time);

private:
    void readConfigFile(const std::string &fileName);
    void initGlobalVariables();
    Eigen::VectorXd stdVec2EigenVec(const std::vector<double> &stdVec);
    bool initEGMCoef(const std::string &filename);

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