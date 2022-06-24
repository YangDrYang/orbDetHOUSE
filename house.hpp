#ifndef HOUSE_H
#define HOUSE_H

#include <Eigen/Dense>

#include <vector>

class HOUSE {

    public:

        // Dynamical model type
        typedef std::function <Eigen::VectorXd (double, double,
            const Eigen::VectorXd&, const Eigen::VectorXd&)> dyn_model;

        // Measurement model type
        typedef std::function <Eigen::VectorXd (double,
            const Eigen::VectorXd&, const Eigen::VectorXd&)> meas_model;

        // Distribution type
        class Dist {

            public:

                int n;

                Eigen::VectorXd mean, skew, kurt;
                Eigen::MatrixXd cov, covL;

                // Generate distribution from sigma points & weights
                Dist(const Eigen::MatrixXd& X, const Eigen::VectorXd& w);

                // Generate zero-mean Gaussian distribution
                Dist(const Eigen::MatrixXd& S);

        };

        // Augmented state sigma point type
        class Sigma {

            public:

                int n_state, n_noise, n_pts;

                Eigen::MatrixXd state, noise;
                Eigen::VectorXd wgt;

                Sigma(const Dist& distX, const Dist& distW, double delta);

        };

        // Constructor
        HOUSE(
            const dyn_model& f_,
            const meas_model& h_,
            int nz_,
            double t0,
            const Dist& distx0,
            const Dist& distw_,
            const Dist& distv_,
            double delta_
        );

        // Prediction step
        void predict(double tp);

        // Update step with one measurement
        void update(const Eigen::VectorXd& z);

        // Run filter for sequence of measurements
        void run(const Eigen::VectorXd& tz, const Eigen::MatrixXd& Z);

        // Dynamical model
        const dyn_model f;

        // Measurement model type
        const meas_model h;

        // Dimensions of state & measurement
        const int nx, nz;

        // Dimensions of process & measurement noise
        const int nw, nv;

        // Distribution of process noise
        const Dist distw;

        // Distribution of measurement noise
        const Dist distv;

        // Minimal weight at mean
        const double delta;

        // History of state estimate distributions
        std::vector<Dist> distx;

        // Times
        std::vector<double> t;

        // Reset filter
        void reset(double t0, const Dist& distx0);

        // Save results
        void save(const std::string& filename);

};

#endif
