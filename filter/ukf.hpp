#ifndef UKF_H
#define UKF_H

#include <functional>
#include <string>
#include <vector>

#include <Eigen/Dense>

class UKF {

    public:

        // Available sigma point types
        enum sig_type {JU, CUT4, CUT6, CUT8};

        // Dynamical model type
        typedef std::function <Eigen::VectorXd (double, double,
            const Eigen::VectorXd&, const Eigen::VectorXd&)> dyn_model;

        // Measurement model type
        typedef std::function <Eigen::VectorXd (double,
            const Eigen::VectorXd&)> meas_model;

        // Dimensions of state, measurement, & process noise
        const int nx, nz, nw;

        // Indicates whether process noise is additive
        const bool addw;

        // Dynamical model: x(tf) = f(ti, tf, x(ti), w)
        const dyn_model f;

        // Measurement model: z = h(t, x) + v
        const meas_model h;

        // Process & measurement noise covariance
        const Eigen::MatrixXd Pww, Pnn;

        // Cholesky decomposition of process noise covariance
        const Eigen::MatrixXd Cww;

        // Standardized sigma points & weights
        const Eigen::MatrixXd Sp, Su;
        const Eigen::VectorXd wp, wu;

        // Number of sigma points in prediction & update step
        const int nsp, nsu;

        // Constructor
        UKF (
            const dyn_model& f_,
            const meas_model& h_,
            bool addw_,
            double t0,
            const Eigen::VectorXd& xm0,
            const Eigen::MatrixXd& Pxx0,
            const Eigen::MatrixXd& Pww_,
            const Eigen::MatrixXd& Pnn_,
            sig_type stype,
            double k);

        // Prediction step
        void predict(double tp);

        // Update step with one measurement
        void update(const Eigen::VectorXd& z);

        // Run filter for sequence of measurements
        void run(const Eigen::VectorXd& tz, const Eigen::MatrixXd& Z);

        // Generate standardized sigma points
        static Eigen::MatrixXd sigmaSt(sig_type stype, int n, double k);

        // Generate sigma point weights
        static Eigen::VectorXd sigmaWt(sig_type stype, int n, double k);

        // Times
        std::vector<double> t;

        // State estimates
        std::vector<Eigen::VectorXd> xest;

        // State estimate covariances
        std::vector<Eigen::MatrixXd> Pxx;

        // Reset filter
        void reset(
            double t0,
            const Eigen::VectorXd& xm0,
            const Eigen::MatrixXd& Pxx0
        );

        // Save results
        void save(const std::string& filename);

        // Directory for CUT files
        static std::string cut_dir;

};

#endif
