#ifndef UT_HPP
#define UT_HPP

#include <functional>
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class UT
{

public:
    // Available sigma point types
    enum sig_type
    {
        JU,
        CUT4,
        CUT6,
        CUT8
    };

    // Coordinate transformation model type
    typedef function<VectorXd(const VectorXd &)> trans_model;

    // Dimensions of state, & process noise
    const int nx, nw;

    // Indicates whether process noise is additive
    const bool addw;

    // Coordinate transformation model
    const trans_model gt;

    // Process noise covariance
    const MatrixXd Pww;

    // Cholesky decomposition of process noise covariance
    const MatrixXd Cww;

    // Standardized sigma points & weights
    const MatrixXd Sp, Su;
    const VectorXd wp, wu;

    // Number of sigma points in prediction step
    const int nsp;

    // Constructor
    UT(
        const trans_model &gt_,
        bool addw_,
        double t0,
        const VectorXd &xm0,
        const MatrixXd &Pxx0,
        const MatrixXd &Pww_,
        sig_type stype,
        double k);

    // Unscented transformation from ECI xi to MEE xf for state and covariance
    void operator()(VectorXd &xx, MatrixXd &Pxx);

    // Unscented transformation from ECI xi to MEE xf for process noise covariance
    void operator()(MatrixXd &Pww);

    // Generate standardized sigma points
    static MatrixXd sigmaSt(sig_type stype, int n, double k);

    // Generate sigma point weights
    static VectorXd sigmaWt(sig_type stype, int n, double k);

    // State estimates
    VectorXd xi, xf;

    // State estimate covariances
    MatrixXd Pi, Pf;

    // Directory for CUT files
    static string cut_dir;
};

#endif
