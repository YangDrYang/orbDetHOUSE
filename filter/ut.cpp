#include "ut.hpp"

#include "eigen_csv.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace std;
using namespace Eigen;

// unscented transformation
void UT::operator()(VectorXd &xx, MatrixXd &Pxx)
{

    MatrixXd Xi(nx, nsp), Xp(nx, nsp);

    Xi = xi.rowwise().replicate(nsp) + Pi.llt().matrixL() * Sp.topRows(nx);

    for (int i = 0; i < nsp; i++)
        Xp.col(i) = gt(Xi.col(i));

    xf = Xp * wp;
    // cout << "transformed state by UT:\t" << endl
    //      << xf << endl;
    // pass on the state
    xx = xf;

    Pf = Xp * wp.asDiagonal() * Xp.transpose() - xf * xf.transpose();
    // cout << "transformed covariance by UT:\t" << endl
    //      << Pf << endl;
    // pass on the covariance
    Pxx = Pf;
}

// Constructor
UT::UT(
    const trans_model &gt_,
    bool addw_,
    double t0,
    const VectorXd &xm0,
    const MatrixXd &Pxx0,
    const MatrixXd &Pww_,
    sig_type stype,
    double k) : nx(xm0.size()),
                nw(Pww_.rows()),
                addw(addw_),
                gt(gt_),
                Pww(Pww_),
                Cww(Pww.llt().matrixL()),
                Sp(sigmaSt(stype, addw ? nx : nx + nw, k)),
                wp(sigmaWt(stype, addw ? nx : nx + nw, k)),
                nsp(Sp.cols())
{
    xi = xm0;
    Pi = Pxx0;
}

// Generate standardized sigma points
MatrixXd UT::sigmaSt(UT::sig_type stype, int n, double k)
{

    MatrixXd S;

    if (stype == JU)
    {
        S.resize(n, 2 * n + 1);
        S.leftCols(n).setIdentity();
        S.rightCols(n).setIdentity();
        S.rightCols(n) *= -1;
        S.col(n).setZero();
        S *= sqrt(n + k);
    }
    else
    {
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
VectorXd UT::sigmaWt(UT::sig_type stype, int n, double k)
{

    VectorXd w;

    if (stype == JU)
    {
        w.resize(2 * n + 1);
        w.setConstant(0.5 / (n + k));
        w(n) = k / (n + k);
    }
    else
    {
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

// Default location for CUT files
string UT::cut_dir = "CUT/";
