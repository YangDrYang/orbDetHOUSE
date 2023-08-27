
#include "eigen_csv.hpp"
#include "srukf.hpp"
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
using namespace Eigen;

MatrixXd cholupdate(MatrixXd matL, MatrixXd matW, double alpha)
{
    int nRows = matL.rows();
    int nColumns = matW.cols();

    MatrixXd matLo;

    for (int i = 0; i < nColumns; ++i)
    {
        matLo = matL;

        double b = 1.0;

        for (int j = 0; j < nRows; ++j)
        {
            const double ljjPow2 = matL(j, j) * matL(j, j);
            const double wjiPow2 = matW(j, i) * matW(j, i);
            const double upsilon = (ljjPow2 * b) + (alpha * wjiPow2);

            matLo(j, j) = std::sqrt(ljjPow2 + ((alpha / b) * wjiPow2));

            for (int k = j + 1; k < nRows; ++k)
            {
                matW(k, i) -= (matW(j, i) / matL(j, j)) * matL(k, j);
                matLo(k, j) = ((matLo(j, j) / matL(j, j)) * matL(k, j)) + (matLo(j, j) * alpha * matW(j, i) * matW(k, i) / upsilon);
            }

            b += alpha * (wjiPow2 / ljjPow2);
        }

        matL = matLo;
    }
    return matLo;
}

// Prediction step
void SRUKF::predict(double tp)
{

    MatrixXd matSxxi, matSxxp, xi(nx, nx), matXi(nx, nsp), matXp(nx, nsp), matW(nw, nsp);
    VectorXd xesti, xestp;

    xesti = xest.back();
    matSxxi = Sxx.back();
    cout << "xesti value: \n"
         << xesti << endl;
    cout << "matSxxi size:\t" << matSxxi.rows() << "\t" << matSxxi.cols() << endl;
    cout << "matSxxi value:\n"
         << matSxxi << endl;

    double ti = t.back();

    // cout << tp << endl
    //      << ti << endl;

    if (tp > ti)
    {
        // control the prediction step, no larger than dtMax
        while (tp > ti + dtMax)
        {
            matXi = xesti.rowwise().replicate(nsp) + matSxxi * Sp.topRows(nx);

            matW.setZero();
            for (int i = 0; i < nsp; i++)
                matXp.col(i) = f(ti, ti + dtMax, matXi.col(i), matW.col(i));
            xestp = matXp * wp;

            MatrixXd matRes(nx, 2 * nx + nx);
            MatrixXd matxestp = xestp.replicate(1, nx);
            VectorXd sqrtWeight1 = wp.segment(0, nx).array().sqrt();
            VectorXd sqrtWeight2 = wp.segment(nx + 1, nx).array().sqrt();
            matRes << (matXp.leftCols(nx) - matxestp) * sqrtWeight1.asDiagonal(),
                (matXp.rightCols(nx) - matxestp) * sqrtWeight2.asDiagonal(), Sww;
            HouseholderQR<MatrixXd> qr(matRes.transpose());
            // Get the upper triangular matrix from the factorization
            MatrixXd matS2 = qr.matrixQR().triangularView<Upper>();
            MatrixXd matS = matS2.block(0, 0, nx, nx).transpose();
            // cout << "matS:\n"
            //      << matS << endl;
            // cout << "matS size:\t" << matS.rows() << "\t" << matS.cols() << endl;
            // cout << "matS2:\n"
            //      << matS2 << endl;
            // cout << "matS2 size:\t" << matS2.rows() << "\t" << matS2.cols() << endl;
            // Update the Cholesky decomposition with the nth vector
            MatrixXd matSxxp = cholupdate(matS, matXp.col(nx) - xestp, wp(nx));
            // LLT<MatrixXd, Lower> lltOfS(matS);
            // lltOfS.rankUpdate(matXp.col(nx) - xestp, wp(nx));
            // MatrixXd matSxxp = matS;
            // cout << "matSxxp:\n"
            //      << matSxxp << endl;
            xesti = xestp;
            ti += dtMax;
            matSxxi = matSxxp;
            // cout << "matSxxp value:\t" << matSxxp << endl;
        }

        matXi = xesti.rowwise().replicate(nsp) + matSxxi * Sp.topRows(nx);
        matW.setZero();
        for (int i = 0; i < nsp; i++)
            matXp.col(i) = f(ti, tp, matXi.col(i), matW.col(i));

        xestp = matXp * wp;
        cout << "mean in SRUKF prediction:\n"
             << xestp << endl;

        MatrixXd matRes(nx, 2 * nx + nx);
        MatrixXd matxestp = xestp.replicate(1, nx);
        VectorXd sqrtWeight1 = wp.segment(0, nx).array().sqrt();
        VectorXd sqrtWeight2 = wp.segment(nx + 1, nx).array().sqrt();
        matRes << (matXp.leftCols(nx) - matxestp) * sqrtWeight1.asDiagonal(),
            (matXp.rightCols(nx) - matxestp) * sqrtWeight2.asDiagonal(), Sww;
        HouseholderQR<MatrixXd> qr(matRes.transpose());
        // Get the upper triangular matrix from the factorization
        MatrixXd matS2 = qr.matrixQR().triangularView<Upper>();
        MatrixXd matS = matS2.block(0, 0, nx, nx).transpose();
        // Update the Cholesky decomposition with the nth vector
        MatrixXd matSxxp = cholupdate(matS, matXp.col(nx) - xestp, wp(nx));
        // LLT<MatrixXd, Lower> lltOfS(matS);
        // lltOfS.rankUpdate(matXp.col(nx) - xestp, wp(nx));
        // MatrixXd matSxxp = matS;
        // cout << "matSxxp size:\t" << matSxxp.rows() << "\t" << matSxxp.cols() << endl;
        cout << "matSxxp value:\n"
             << matSxxp << endl;
        cout << "predicted covariance:\n"
             << matSxxp * matSxxp.transpose() << endl;

        t.push_back(tp);
        xest.push_back(xestp);
        // Saved as a lower trangular matrix
        Sxx.push_back(matSxxp);
        Pxx.push_back(matSxxp * matSxxp.transpose());
    }
}

// Update step with one measurement
void SRUKF::update(const VectorXd &z)
{
    MatrixXd matSxxp(nx, nx), matSxxu(nx, nx), matX(nx, nsu), matZ(nz, nsu), matPzx(nx, nz), matK(nx, nz);
    VectorXd xestp(nx), xestu(nx), zm(nz);

    xestp = xest.back();
    matSxxp = Sxx.back();

    cout << "dimensions of meas: \t" << nz << endl;
    cout << "predicted xest: \n"
         << xestp << endl;
    cout << "predicted Pxx: \n"
         << matSxxp * matSxxp.transpose() << endl;

    double tz = t.back();

    // cout << "Su" << Su << endl;
    matX = xestp.rowwise().replicate(nsu) + matSxxp * Su;

    for (int i = 0; i < nsu; i++)
        matZ.col(i) = h(tz, matX.col(i));

    // VectorXd res = z - matZ.col(0);
    // cout << "residuals: " << res(0) << "\t" << res(1) << endl;

    zm = matZ * wu;

    MatrixXd matzm = zm.replicate(1, nx); // Repeat zm horizontally to match the number of columns of Z
    MatrixXd matRes(nz, 2 * nx + nz);
    VectorXd sqrtWeight1 = wu.segment(0, nx).array().sqrt();
    VectorXd sqrtWeight2 = wu.segment(nx + 1, nx).array().sqrt();
    matRes << (matZ.leftCols(nx) - matzm) * sqrtWeight1.asDiagonal(),
        (matZ.rightCols(nx) - matzm) * sqrtWeight2.asDiagonal(), Snn;
    HouseholderQR<MatrixXd> qr(matRes.transpose());
    // Get the upper triangular matrix from the factorization
    MatrixXd matS2 = qr.matrixQR().triangularView<Upper>();
    // Get the lower triangular matrix
    MatrixXd matS = matS2.block(0, 0, nz, nz).transpose();
    // cout << "matS:\n"
    //      << matS << endl;
    // cout << "matS size:\t" << matS.rows() << "\t" << matS.cols() << endl;
    // // cout << "matS2:\n"
    //      << matS2 << endl;
    // cout << "matS2 size:\t" << matS2.rows() << "\t" << matS2.cols() << endl;
    // Update the Cholesky decomposition with the nth vector
    MatrixXd matSzz = cholupdate(matS, matZ.col(nx) - zm, wu(nx));
    // LLT<MatrixXd, Lower> lltOfS(matS);
    // lltOfS.rankUpdate(matZ.col(nx) - zm, wu(nx));
    // matSzz = matS;
    // cout << matZ.col(nx) - zm << "\t" << wu(nx) << endl;
    cout << "matSzz:\n"
         << matSzz << endl;
    cout << "measurement covariance Pzz: \n"
         << matSzz * matSzz.transpose() << endl;

    matPzx = matZ * wu.asDiagonal() * matX.transpose() - zm * xestp.transpose();

    cout << "cross covariance Pzx: \n"
         << matPzx << endl;
    // Kalman Gain
    MatrixXd matK1 = matSzz.householderQr().solve(matPzx);
    MatrixXd matK2 = matSzz.transpose().householderQr().solve(matK1);
    matK = matK2.transpose();
    cout << "Kalman gain:\n"
         << matK << endl;

    // Update the state
    xestu = xestp + matK * (z - zm);
    cout << "updated xest: \n"
         << xestu << endl;

    // Create a temporary matrix to hold the result of K * Szz
    MatrixXd matU = matK * matSzz;
    cout << "matU:\n"
         << matU << endl;
    // Update the covariance
    cout << "Sxxp:\n"
         << matSxxp << endl;
    // LLT<MatrixXd, Lower> lltOfSxxp(matSxxp);
    // for (int j = 0; j < matU.cols(); j++)
    // {
    //     lltOfSxxp.rankUpdate(matU.col(j), -1.0);
    //     cout << "rankUpdate middle result:\n"
    //          << j << matSxxp << endl;
    // }
    // matSxxu = matSxxp;
    matSxxu = cholupdate(matSxxp, matU, -1.0);
    cout << "Sxxu:\n"
         << matSxxu << endl;
    cout << "updated Pxx: \n"
         << matSxxu * matSxxu.transpose() << endl;

    xest.back() = xestu;
    Sxx.back() = matSxxu;
    Pxx.back() = matSxxu * matSxxu.transpose();
    cout << Sxx[0] << endl;
}

// Run filter for sequence of measurements
void SRUKF::run(const VectorXd &tz, const MatrixXd &Z)
{
    for (int i = 0; i < tz.size(); i++)
    {
        cout << "the " << i << "th epoch" << endl;
        predict(tz(i));
        // only update if given a measurment
        if (abs(Z(1, i)) <= M_PI * 2)
        {
            update(Z.col(i));
        }
    }
}

// Constructor
SRUKF::SRUKF(
    const dyn_model &f_,
    const meas_model &h_,
    bool addw_,
    double t0,
    double dtMax_,
    const VectorXd &xm0,
    const MatrixXd &Pxx0,
    const MatrixXd &Pww_,
    const MatrixXd &Pnn_,
    sig_type stype,
    double k) : UKF(f_, h_, addw_, t0, dtMax_, xm0, Pxx0, Pww_, Pnn_, stype, k),
                nx(xm0.size()),
                nz(Pnn_.rows()),
                nw(Pww_.rows()),
                addw(addw_),
                f(f_),
                dtMax(dtMax_),
                h(h_),
                Sww(Pww_.llt().matrixL()),
                Snn(Pnn_.llt().matrixL()),
                Sp(sigmaSt(stype, addw ? nx : nx + nw, k)),
                Su(sigmaSt(stype, nx, k)),
                wp(sigmaWt(stype, addw ? nx : nx + nw, k)),
                wu(sigmaWt(stype, nx, k)),
                nsp(Sp.cols()),
                nsu(Su.cols())
{
    // // cholesky factorization to get matrix Pk square-root
    // LLT<MatrixXd> lltOfPxx(Pxx0);
    // MatrixXd tmpMat = lltOfPxx.matrixL();
    t.push_back(t0);
    xest.push_back(xm0);
    Pxx.push_back(Pxx0);
    Sxx.push_back(Pxx0.llt().matrixL());
}

// Generate standardized sigma points
MatrixXd SRUKF::sigmaSt(UKF::sig_type stype, int n, double k)
{
    return UKF::sigmaSt(stype, n, k);
}

// Generate sigma point weights
VectorXd SRUKF::sigmaWt(UKF::sig_type stype, int n, double k)
{
    return UKF::sigmaWt(stype, n, k);
}

// Reset filter
void SRUKF::reset(double t0, const VectorXd &xm0, const MatrixXd &Pxx0)
{
    UKF::reset(t0, xm0, Pxx0);
    Sxx.push_back(Pxx0.llt().matrixL());
}

// Save results for SRUKF
// void SRUKF::save(const string &filename, string stateType)
// {
//     // You can reuse the save function from the base class UKF
//     UKF::save(filename, stateType);
// }
void SRUKF::save(const string &filename, string stateType)
{

    int steps = xest.size();

    MatrixXd table(steps, nx * (nx + 1) + 1);

    for (int k = 0; k < steps; k++)
    {

        table(k, 0) = t[k];

        if (stateType == "eci")
            table.row(k).segment(1, nx) = xest[k];
        else if (stateType == "mee")
            table.row(k).segment(1, nx) = mee2eci(xest[k], GM_Earth);

        table.row(k).tail(nx * nx) = Pxx[k].reshaped(1, nx * nx);

        // // Use Eigen::Map to copy the covariance matrix elements directly
        // Map<MatrixXd> covMatrixMap(table.row(k).tail(nx * nx).data(), nx, nx);
        // covMatrixMap = Pxx[k];
    }

    vector<string> header(nx * (nx + 1) + 1);
    header[0] = "TIME";
    for (int i = 1; i <= nx; i++)
    {
        header[i] = "EST X";
        header[i + nx * i] = "STD X";
        header[i] += to_string(i);
        header[i + nx * i] += to_string(i);
    }

    EigenCSV::write(table, header, filename);
}

// Default location for CUT files
string SRUKF::cut_dir = "CUT/";
