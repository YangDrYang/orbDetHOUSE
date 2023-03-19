#include "eigen_csv.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

// Generate filename for state estimates
std::string est_filename(const std::string& filter, int i, int j) {
    std::string f = "out/ct_est_";
    f += filter;
    f += "_";
    f += std::to_string(i+1);
    f += "_";
    f += std::to_string(j+1);
    f += ".csv";
    return f;
}

int main() {

    using namespace std;
    using namespace Eigen;

    string filename;

    // Number of tests
    int Ntest = 100;

    // Number of filters & filter names
    int Nfilters = 5;
    vector<string> filter_names(Nfilters);
    filter_names[0] = "ukf";
    filter_names[1] = "cut4";
    filter_names[2] = "cut6";
    filter_names[3] = "cut8";
    filter_names[4] = "house";

    // Read times
    VectorXd T;
    EigenCSV::read("out/times.csv", false, true, T);
    int Ntimes = T.size();

    // Results table
    MatrixXd table(Ntimes, 1+3*Nfilters);

    for (int l = 0; l < Nfilters; l++) {

        VectorXd
        pos_rmse(Ntimes), vel_rmse(Ntimes), ome_rmse(Ntimes),
        pos_avge(Ntimes), vel_avge(Ntimes), ome_avge(Ntimes),
        pos_stde(Ntimes), vel_stde(Ntimes), ome_stde(Ntimes);

        for (int i = 0; i < Ntimes; i++) {

            MatrixXd tabtru, pos_tru, vel_tru, ome_tru;

            filename = "out/ct_true_";
            filename += to_string(i+1);
            filename += ".csv";

            EigenCSV::read(filename, true, true, tabtru);

            int nt = tabtru.rows();

            pos_tru.resize(nt, 2);
            vel_tru.resize(nt, 2);
            ome_tru.resize(nt, 1);

            pos_tru << tabtru.col(1), tabtru.col(3);
            vel_tru << tabtru.col(2), tabtru.col(4);
            ome_tru = tabtru.col(5);

            MatrixXd pos_err(nt*Ntest, 2), vel_err(nt*Ntest, 2), ome_err(nt*Ntest, 1);

            for (int j = 0; j < Ntest; j++) {

                MatrixXd tabest, pos_est(nt,2), vel_est(nt,2), ome_est(nt,1);

                filename = est_filename(filter_names[l], i, j);

                EigenCSV::read(filename, true, true, tabest);

                pos_est << tabest.col(1), tabest.col(3);
                vel_est << tabest.col(2), tabest.col(4);
                ome_est = tabest.col(5);

                pos_err.block(j*nt, 0, nt, 2) = pos_est - pos_tru;
                vel_err.block(j*nt, 0, nt, 2) = vel_est - vel_tru;
                ome_err.block(j*nt, 0, nt, 1) = ome_est - ome_tru;

            }

            pos_rmse(i) = sqrt(pos_err.array().square().mean() * 2);
            vel_rmse(i) = sqrt(vel_err.array().square().mean() * 2);
            ome_rmse(i) = sqrt(ome_err.array().square().mean());

            pos_avge(i) = pos_err.mean();
            vel_avge(i) = vel_err.mean();
            ome_avge(i) = ome_err.mean();

            pos_stde(i) = sqrt((pos_err.array() - pos_avge(i)).square().mean());
            vel_stde(i) = sqrt((vel_err.array() - pos_avge(i)).square().mean());
            ome_stde(i) = sqrt((ome_err.array() - ome_avge(i)).square().mean());

        }

        MatrixXd tabfil(Ntimes, 10);
        tabfil << T, pos_rmse, pos_avge, pos_stde,
                     vel_rmse, vel_avge, vel_stde,
                     ome_rmse, ome_avge, ome_stde;

        vector<string> header(10);
        header[0] = "TIME";
        header[1] = "POS RMSE"; header[2] = "POS AVGE"; header[3] = "POS STDE";
        header[4] = "VEL RMSE"; header[5] = "VEL AVGE"; header[6] = "VEL STDE";
        header[7] = "OME RMSE"; header[8] = "OME AVGE"; header[9] = "OME STDE";

        filename = "out/ct_err_";
        filename += filter_names[l];
        filename += ".csv";

        EigenCSV::write(tabfil, header, filename);

    }

    return 0;

}

