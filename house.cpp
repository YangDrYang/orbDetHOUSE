#include "house.hpp"

#include "eigen_csv.hpp"

#include <cmath>
#include <iostream>

using namespace Eigen;

// Prediction step
void HOUSE::predict(double tp) {

    double ti = t.back();

    if (tp > ti) {

        Sigma sig(distx.back(), distw, delta);

        MatrixXd Xp(nx, sig.n_pts);

        for (int i = 0; i < sig.n_pts; i++)
            Xp.col(i) = f(ti, tp, sig.state.col(i), sig.noise.col(i));

        Dist distXp(Xp, sig.wgt);

        distx.push_back(distXp);
        t.push_back(tp);

    }

}

// Update step with one measurement
void HOUSE::update(const VectorXd& z) {

    double tz = t.back();

    Sigma sig(distx.back(), distv, delta);

    MatrixXd Z(nz, sig.n_pts), Pzz(nz, nz), Pzx(nz, nx), K(nx, nz),
        Xu(nx, sig.n_pts);
    VectorXd xm, zm;

    for (int i = 0; i < sig.n_pts; i++)
        Z.col(i) = h(tz, sig.state.col(i), sig.noise.col(i));

    xm = distx.back().mean;
    zm = Z * sig.wgt;

    Pzz = Z * sig.wgt.asDiagonal() * Z.transpose() - zm*zm.transpose();
    Pzx = Z * sig.wgt.asDiagonal() * sig.state.transpose() - zm*xm.transpose();

    K = Pzz.llt().solve(Pzx).transpose();

    Xu = sig.state - K * (Z.colwise() - zm);

    Dist distXu(Xu, sig.wgt);

    distXu.mean = xm + K * (z - zm);

    distx.back() = distXu;

}

// Run filter for sequence of measurements
void HOUSE::run(const VectorXd& tz, const MatrixXd& Z) {
    for (int i = 0; i < tz.size(); i++) {
        predict(tz(i));
        update(Z.col(i));
    }
}

// Sigma point constructor
HOUSE::Sigma::Sigma(const Dist& distX, const Dist& distW, double delta) {

    double S, m, mmin;
    int nx, nw;
    nx = distX.n;
    nw = distW.n;

    n_state = nx;
    n_noise = nw;
    n_pts = 2*(nx+nw)+1;

    VectorXd sx, kx, sw, kw, ax, bx, cx, aw, bw, cw;

    sx = distX.skew;
    sw = distW.skew;
    kx = distX.kurt;
    kw = distW.kurt;

    mmin = (nx + nw) / (1 - delta);

    for (int i = 0; i < nx; i++) {
        m = kx(i) - sx(i) * sx(i);
        if (m < mmin)
            kx(i) += mmin - m;
    }

    for (int i = 0; i < nw; i++) {
        m = kw(i) - sw(i) * sw(i);
        if (m < mmin)
            kw(i) += mmin - m;
    }

    cx = (4 * kx.array() - 3 * sx.array().square()).sqrt();
    cw = (4 * kw.array() - 3 * sw.array().square()).sqrt();

    ax = (sx + cx) / 2;
    aw = (sw + cw) / 2;

    bx = ax - sx;
    bw = aw - sw;

    wgt.resize(n_pts);

    wgt.segment(1, nx)    = ax.cwiseProduct(cx).cwiseInverse();
    wgt.segment(1+nx, nx) = bx.cwiseProduct(cx).cwiseInverse();

    wgt.segment(1+2*nx, nw)    = aw.cwiseProduct(cw).cwiseInverse();
    wgt.segment(1+2*nx+nw, nw) = bw.cwiseProduct(cw).cwiseInverse();

    S =  (kx.array() - sx.array().square()).inverse().sum()
       + (kw.array() - sw.array().square()).inverse().sum();

    wgt(0) = 1 - S;

    state = distX.mean.rowwise().replicate(n_pts);
    noise = distW.mean.rowwise().replicate(n_pts);

    state.block(0, 1,    nx, nx) += distX.covL * ax.asDiagonal();
    state.block(0, 1+nx, nx, nx) -= distX.covL * bx.asDiagonal();

    noise.block(0, 1+2*nx,    nw, nw) += distW.covL * aw.asDiagonal();
    noise.block(0, 1+2*nx+nw, nw, nw) -= distW.covL * bw.asDiagonal();

}

// Constructor
HOUSE::HOUSE(
        const dyn_model& f_,
        const meas_model& h_,
        int nz_,
        double t0,
        const Dist& distx0,
        const Dist& distw_,
        const Dist& distv_,
        double delta_
    ) :
        f(f_), h(h_),
        nx(distx0.n), nz(nz_),
        nw(distw_.n), nv(distv_.n),
        distw(distw_), distv(distv_),
        delta(delta_)
    {

    distx.push_back(distx0);
    t.push_back(t0);

}


// Generate distribution from sigma points & weights
HOUSE::Dist::Dist(const MatrixXd& X, const VectorXd& w) {

    n = X.rows();

    mean = X * w;

    cov = X * w.asDiagonal() * X.transpose() - mean * mean.transpose();

    covL = cov.llt().matrixL();

    MatrixXd Xstd = covL.triangularView<Lower>().solve(X.colwise() - mean);

    skew = Xstd.array().pow(3).matrix() * w;
    kurt = Xstd.array().pow(4).matrix() * w;

}

// Generate zero-mean Gaussian distribution
HOUSE::Dist::Dist(const MatrixXd& S) {

    n = S.rows();

    mean.setZero(n);

    cov = S;

    covL = S.llt().matrixL();

    skew.setZero(n);
    kurt.setConstant(n, 3);

}

// Reset filter
void HOUSE::reset(double t0, const Dist& distx0) {

    t.clear();
    distx.clear();

    t.push_back(t0);
    distx.push_back(distx0);

}

// Save results
void HOUSE::save(const std::string& filename) {

    using namespace std;

    int steps = distx.size();

    MatrixXd table(steps, 2*nx+1);

    for (int k = 0; k < steps; k++) {

        table(k,0) = t[k];

        table.row(k).segment(1, nx) = distx[k].mean;

        table.row(k).tail(nx) = distx[k].cov.diagonal().cwiseSqrt();

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
