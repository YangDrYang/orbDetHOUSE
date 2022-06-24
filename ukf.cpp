#include "ukf.hpp"

#include "eigen_csv.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace std;
using namespace Eigen;

// Prediction step
void UKF::predict(double tp) {

    MatrixXd Pxxi(nx,nx), Pxxp(nx,nx), Xi(nx,nsp), Xp(nx,nsp), W(nw,nsp);
    VectorXd xesti, xestp;

    xesti = xest.back();
    Pxxi = Pxx.back();

    double ti = t.back();

    if (tp > ti) {

        Xi = xesti.rowwise().replicate(nsp) + Pxxi.llt().matrixL() * Sp.topRows(nx);

        if (addw)
            W.setZero();
        else
            W = Cww * Sp.bottomRows(nw);

        for (int i = 0; i < nsp; i++)
            Xp.col(i) = f(ti, tp, Xi.col(i), W.col(i));

        xestp = Xp * wp;

        Pxxp = Xp * wp.asDiagonal() * Xp.transpose() - xestp * xestp.transpose();

        if (addw)
            Pxxp += Pww;

        t.push_back(tp);
        xest.push_back(xestp);
        Pxx.push_back(Pxxp);

    }

}

// Update step with one measurement
void UKF::update(const Eigen::VectorXd& z) {

    MatrixXd Pxxp(nx,nx), Pxxu(nx,nx), X(nx,nsu), Z(nz,nsu), Pzz(nz,nz),
        Pzx(nz,nx), K(nx,nz);
    VectorXd xestp(nx), xestu(nx), zm(nz);

    xestp = xest.back();
    Pxxp = Pxx.back();

    double tz = t.back();

    X = xestp.rowwise().replicate(nsu) + Pxxp.llt().matrixL() * Su;

    for (int i = 0; i < nsu; i++)
        Z.col(i) = h(tz, X.col(i));

    zm = Z * wu;

    Pzz = Z * wu.asDiagonal() * Z.transpose() - zm * zm.transpose() + Pnn;

    Pzx = Z * wu.asDiagonal() * X.transpose() - zm * xestp.transpose();

    K = Pzz.llt().solve(Pzx).transpose();

    xestu = xestp + K * (z - zm);

    Pxxu = Pxxp - K * Pzx;

    xest.back() = xestu;
    Pxx.back() = Pxxu;

}

// Run filter for sequence of measurements
void UKF::run(const Eigen::VectorXd& tz, const Eigen::MatrixXd& Z) {
    for (int i = 0; i < tz.size(); i++) {
        predict(tz(i));
        update(Z.col(i));
    }
}

// Constructor
UKF::UKF (
        const dyn_model& f_,
        const meas_model& h_,
        bool addw_,
        double t0,
        const Eigen::VectorXd& xm0,
        const Eigen::MatrixXd& Pxx0,
        const Eigen::MatrixXd& Pww_,
        const Eigen::MatrixXd& Pnn_,
        sig_type stype,
        double k
    ) :
        nx(xm0.size()),
        nz(Pnn_.rows()),
        nw(Pww_.rows()),
        addw(addw_),
        f(f_),
        h(h_),
        Pww(Pww_),
        Pnn(Pnn_),
        Cww(Pww.llt().matrixL()),
        Sp(sigmaSt(stype, addw ? nx : nx+nw, k)),
        Su(sigmaSt(stype, nx, k)),
        wp(sigmaWt(stype, addw ? nx : nx+nw, k)),
        wu(sigmaWt(stype, nx, k)),
        nsp(Sp.cols()),
        nsu(Su.cols())
    {

    t.push_back(t0);
    xest.push_back(xm0);
    Pxx.push_back(Pxx0);

}

// Generate standardized sigma points
Eigen::MatrixXd UKF::sigmaSt(UKF::sig_type stype, int n, double k) {

    MatrixXd S;

    if (stype == JU) {
        S.resize(n, 2*n+1);
        S.leftCols(n).setIdentity();
        S.rightCols(n).setIdentity();
        S.rightCols(n) *= -1;
        S.col(n).setZero();
        S *= sqrt(n + k);
    } else {
        int order;
        if (stype == CUT8)
            order = 8;
        else if (stype == CUT6)
            order = 6;
        else
            order = 4;
        string file = cut_dir;
        file += "pts_";
        file += to_string(order);
        file += "_";
        file += to_string(n);
        file += ".csv";
        EigenCSV::read(file, false, true, S);
    }

    return S;

}

// Generate sigma point weights
Eigen::VectorXd UKF::sigmaWt(UKF::sig_type stype, int n, double k) {

    VectorXd w;

    if (stype == JU) {
        w.resize(2*n+1);
        w.setConstant(0.5 / (n + k));
        w(n) = k / (n + k);
    } else {
        int order;
        if (stype == CUT8)
            order = 8;
        else if (stype == CUT6)
            order = 6;
        else
            order = 4;
        string file = cut_dir;
        file += "wts_";
        file += to_string(order);
        file += "_";
        file += to_string(n);
        file += ".csv";
        EigenCSV::read(file, false, true, w);
    }

    return w;
}

// Reset filter
void UKF::reset(
            double t0,
            const Eigen::VectorXd& xm0,
            const Eigen::MatrixXd& Pxx0
        ) {

    t.clear();
    xest.clear();
    Pxx.clear();

    t.push_back(t0);
    xest.push_back(xm0);
    Pxx.push_back(Pxx0);

}

// Save results
void UKF::save(const string& filename) {

    int steps = xest.size();

    MatrixXd table(steps, 2*nx+1);

    for (int k = 0; k < steps; k++) {

        table(k,0) = t[k];

        table.row(k).segment(1, nx) = xest[k];

        table.row(k).tail(nx) = Pxx[k].diagonal().cwiseSqrt();

    }

    vector<string> header(2*nx+1);
    header[0] = "TIME";
    for (int i = 1; i <= nx; i++) {
        header[i]    = "EST X";
        header[i+nx] = "STD X";
        header[i]    += to_string(i);
        header[i+nx] += to_string(i);
    }

    EigenCSV::write(table, header, filename);

}

// Default location for CUT files
string UKF::cut_dir = "CUT/";
