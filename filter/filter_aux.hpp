#ifndef FILTER_AUX_H
#define FILTER_AUX_H

#include "Eigen/Dense"
#include <string>
#include <vector>
using namespace Eigen;
using namespace std;

void save_rmse(
    const VectorXd &t,
    const vector<vector<VectorXd>> &xtru,
    const vector<vector<VectorXd>> &xest,
    const std::string &filename);

void save_abs_err(
    const VectorXd &t,
    const vector<vector<VectorXd>> &xtru,
    const vector<vector<VectorXd>> &xest,
    const std::string &filename);

void save_abs_err_lump(
    const VectorXd &t,
    const vector<vector<VectorXd>> &xtru,
    const vector<vector<VectorXd>> &xest,
    const std::string &filename);

#endif
