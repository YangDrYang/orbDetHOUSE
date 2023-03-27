#ifndef FILTER_AUX_H
#define FILTER_AUX_H

#include <Eigen/Dense>
#include <string>
#include <vector>

void save_rmse(
        const Eigen::VectorXd& t,
        const std::vector<std::vector<Eigen::VectorXd>>& xtru,
        const std::vector<std::vector<Eigen::VectorXd>>& xest,
        const std::string& filename
    );

void save_abs_err(
        const Eigen::VectorXd& t,
        const std::vector<std::vector<Eigen::VectorXd>>& xtru,
        const std::vector<std::vector<Eigen::VectorXd>>& xest,
        const std::string& filename
    );

void save_abs_err_lump(
        const Eigen::VectorXd& t,
        const std::vector<std::vector<Eigen::VectorXd>>& xtru,
        const std::vector<std::vector<Eigen::VectorXd>>& xest,
        const std::string& filename
    );

#endif
