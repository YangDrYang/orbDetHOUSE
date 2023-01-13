#include "filter_aux.hpp"

#include "eigen_csv.hpp"

using namespace std;
using namespace Eigen;

void save_rmse(
        const VectorXd& t,
        const vector<vector<VectorXd>>& xtru,
        const vector<vector<VectorXd>>& xest,
        const string& filename
    ) {

    int steps, trials, nx;

    steps = t.size();
    trials = xest.size();
    nx = xtru.front().front().size();

    MatrixXd err(nx, trials), table(steps, nx+1);

    table.col(0) = t;

    for (int k = 0; k < steps; k++) {

        for (int j = 0; j < trials; j++)
            err.col(j) = xest[j][k] - xtru[j][k];

        table.row(k).tail(nx) = err.array().square().rowwise().mean().sqrt();

    }

    vector<string> header(nx+1);
    header[0] = "TIME";
    for (int i = 1; i <= nx; i++) {
        header[i] = "RMSE X";
        header[i] += to_string(i);
    }

    EigenCSV::write(table, header, filename);

}

void save_abs_err(
        const VectorXd& t,
        const vector<vector<VectorXd>>& xtru,
        const vector<vector<VectorXd>>& xest,
        const string& filename
    ) {

    int steps, trials, nx;

    steps = t.size();
    trials = xest.size();
    nx = xtru.front().front().size();

    MatrixXd err(nx, trials), table(steps, 3*nx+1);
    VectorXd avg_err(nx), std_err(nx), lb_err(nx), ub_err(nx);

    table.col(0) = t;

    for (int k = 0; k < steps; k++) {

        for (int j = 0; j < trials; j++)
            err.col(j) = xest[j][k] - xtru[j][k];

        avg_err = err.rowwise().mean();

        std_err = (err.colwise() - avg_err).cwiseAbs2().rowwise().mean()
            .cwiseSqrt();

        ub_err = avg_err + std_err;
        lb_err = avg_err - std_err;

        table.row(k).segment(1, nx) = avg_err;
        table.row(k).segment(1+nx, nx) = lb_err;
        table.row(k).segment(1+2*nx, nx) = ub_err;

    }

    vector<string> header(3*nx+1);
    header[0] = "TIME";
    for (int i = 1; i <= nx; i++) {
        header[i] = "AVG ERR X";
        header[i] += to_string(i);
    }
    for (int i = 1; i <= nx; i++) {
        header[nx+i] = "LB ERR X";
        header[nx+i] += to_string(i);
    }
    for (int i = 1; i <= nx; i++) {
        header[2*nx+i] = "UB ERR X";
        header[2*nx+i] += to_string(i);
    }

    EigenCSV::write(table, header, filename);

}

void save_abs_err_lump(
        const Eigen::VectorXd& t,
        const std::vector<std::vector<Eigen::VectorXd>>& xtru,
        const std::vector<std::vector<Eigen::VectorXd>>& xest,
        const std::string& filename
    ) {

    int steps, trials, nx;

    steps = t.size();
    trials = xest.size();
    nx = xtru.front().front().size();

    MatrixXd table(steps*trials, nx);

    for (int i = 0; i < trials; i++)
        for (int j = 0; j < steps; j++)
            table.row(i*steps+j) = (xest[i][j] - xtru[i][j]).cwiseAbs();

    vector<string> header(nx);
    for (int i = 0; i < nx; i++) {
        header[i] = "DX";
        header[i] += to_string(i+1);
    }

    EigenCSV::write(table, header, filename);

}

