#include "srhouse.hpp"
#include "house.hpp"
#include "srukf.hpp"
#include "eigen_csv.hpp"

#include <cmath>
#include <iostream>

using namespace Eigen;

Dist::Dist()
{
}

// Sigma point constructor
Sigma::Sigma(const Dist &distX, const Dist &distW)
{

    double S;
    int nx, nw;
    nx = distX.n;
    nw = distW.n;

    n_state = nx;
    n_noise = nw;
    n_pts = 2 * (nx + nw) + 1;

    VectorXd sx, kx, sw, kw, ax, bx, cx, aw, bw, cw;

    sx = distX.skew;
    sw = distW.skew;
    kx = distX.kurt;
    kw = distW.kurt;

    // mmin = (nx + nw) / (1 - delta);

    // for (int i = 0; i < nx; i++)
    // {
    //     m = kx(i) - sx(i) * sx(i);
    //     if (m < mmin)
    //         kx(i) += mmin - m;
    // }

    // for (int i = 0; i < nw; i++)
    // {
    //     m = kw(i) - sw(i) * sw(i);
    //     if (m < mmin)
    //         kw(i) += mmin - m;
    // }

    cx = (4 * kx.array() - 3 * sx.array().square()).sqrt();
    cw = (4 * kw.array() - 3 * sw.array().square()).sqrt();

    ax = (sx + cx) / 2;
    aw = (sw + cw) / 2;

    bx = ax - sx;
    bw = aw - sw;

    wgt.resize(n_pts);

    wgt.segment(1, nx) = ax.cwiseProduct(cx).cwiseInverse();
    wgt.segment(1 + nx, nx) = bx.cwiseProduct(cx).cwiseInverse();

    wgt.segment(1 + 2 * nx, nw) = aw.cwiseProduct(cw).cwiseInverse();
    wgt.segment(1 + 2 * nx + nw, nw) = bw.cwiseProduct(cw).cwiseInverse();

    S = (kx.array() - sx.array().square()).inverse().sum() + (kw.array() - sw.array().square()).inverse().sum();

    wgt(0) = 1 - S;

    state = distX.mean.rowwise().replicate(n_pts);
    noise = distW.mean.rowwise().replicate(n_pts);

    state.block(0, 1, nx, nx) += distX.covL * ax.asDiagonal();
    state.block(0, 1 + nx, nx, nx) -= distX.covL * bx.asDiagonal();

    noise.block(0, 1 + 2 * nx, nw, nw) += distW.covL * aw.asDiagonal();
    noise.block(0, 1 + 2 * nx + nw, nw, nw) -= distW.covL * bw.asDiagonal();
}

// Prediction step
void SRHOUSE::predict(double tp)
{

    double ti = t.back();

    if (tp > ti)
    {
        Dist distxi = distx.back();
        while (tp > ti + dtMax)
        {
            // Sigma sig(distxi, distw);
            Sigma sig(distxi, distw, 0);

            MatrixXd Xp(nx, sig.n_pts);

            for (int i = 0; i < sig.n_pts; i++)
                Xp.col(i) = f(ti, ti + dtMax, sig.state.col(i), sig.noise.col(i));

            // Dist distXp(Xp, sig.wgt);
            Dist distXp;
            distXp.n = nx;
            distXp.mean = Xp * sig.wgt;
            MatrixXd matRes(nx, 4 * nx);
            MatrixXd matxestp = distXp.mean.replicate(1, nx);
            matRes << sqrt(sig.wgt(1)) * (Xp.block(0, 1, nx, nx) - matxestp), sqrt(sig.wgt(nx + 1)) * (Xp.block(0, nx + 1, nx, nx) - matxestp), sqrt(sig.wgt(2 * nx + 1)) * (Xp.block(0, 2 * nx + 1, nx, 2 * nx) - distXp.mean.replicate(1, 2 * nx));
            HouseholderQR<MatrixXd> qr(matRes.transpose());
            MatrixXd matS2 = qr.matrixQR().triangularView<Upper>();
            MatrixXd matS = matS2.block(0, 0, nx, nx).transpose();
            distXp.covL = cholupdate(matS, Xp.col(0) - distXp.mean, sig.wgt(0));
            distXp.cov = distXp.covL * distXp.covL.transpose();
            MatrixXd Xstd = distXp.covL.triangularView<Lower>().solve(Xp.colwise() - Xp * sig.wgt); // standardised states at the sigma points,  covariance, A7/B10
            distXp.skew = Xstd.array().pow(3).matrix() * sig.wgt;                                   // skewness of the standardised state, Eq. A8/B11
            distXp.kurt = Xstd.array().pow(4).matrix() * sig.wgt;                                   // kurtosis of the standardised state, Eq. A9/B12

            distxi = distXp;
            ti += dtMax;
        }

        // Sigma sig(distxi, distw);
        Sigma sig(distxi, distw, 0);

        MatrixXd Xp(nx, sig.n_pts);

        for (int i = 0; i < sig.n_pts; i++)
            Xp.col(i) = f(ti, tp, sig.state.col(i), sig.noise.col(i));

        // Dist distXp(Xp, sig.wgt);
        Dist distXp;
        distXp.n = nx;
        distXp.mean = Xp * sig.wgt;
        MatrixXd matRes(nx, 4 * nx);
        MatrixXd matxestp = distXp.mean.replicate(1, nx);
        matRes << sqrt(sig.wgt(1)) * (Xp.block(0, 1, nx, nx) - matxestp), sqrt(sig.wgt(nx + 1)) * (Xp.block(0, nx + 1, nx, nx) - matxestp), sqrt(sig.wgt(2 * nx + 1)) * (Xp.block(0, 2 * nx + 1, nx, 2 * nx) - distXp.mean.replicate(1, 2 * nx));
        HouseholderQR<MatrixXd> qr(matRes.transpose());
        MatrixXd matS2 = qr.matrixQR().triangularView<Upper>();
        MatrixXd matS = matS2.block(0, 0, nx, nx).transpose();
        distXp.covL = cholupdate(matS, Xp.col(0) - distXp.mean, sig.wgt(0));
        distXp.cov = distXp.covL * distXp.covL.transpose();
        MatrixXd Xstd = distXp.covL.triangularView<Lower>().solve(Xp.colwise() - Xp * sig.wgt); // standardised states at the sigma points,  covariance, A7/B10
        distXp.skew = Xstd.array().pow(3).matrix() * sig.wgt;                                   // skewness of the standardised state, Eq. A8/B11
        distXp.kurt = Xstd.array().pow(4).matrix() * sig.wgt;                                   // kurtosis of the standardised state, Eq. A9/B12

        cout << "mean in SRHOUSE prediction:\t" << endl
             << distXp.mean << endl;
        cout << "covariance in SRHOUSE prediction:\t" << endl
             << distXp.cov << endl;
        cout << "covariance lower triangle in SRHOUSE prediction:\t" << endl
             << distXp.covL << endl;

        distx.push_back(distXp);
        t.push_back(tp);
    }
}

// Update step with one measurement
void SRHOUSE::update(const VectorXd &z)
{
    double tz = t.back();

    // Sigma sig(distx.back(), distv);
    Sigma sig(distx.back(), distv, 0);

    MatrixXd Z(nz, sig.n_pts), Pzz(nz, nz), Pzx(nz, nx), K(nx, nz),
        Xu(nx, sig.n_pts);
    VectorXd xm, zm;

    for (int i = 0; i < sig.n_pts; i++)
        Z.col(i) = h(tz, sig.state.col(i), sig.noise.col(i)); // Eq. B3

    VectorXd res = z - Z.col(0);
    cout << "residuals: " << res(0) << "\t" << res(1) << endl;

    cout << "weight:\t" << sig.wgt.transpose() << endl;
    xm = distx.back().mean;
    zm = Z * sig.wgt;

    MatrixXd matRes(nz, 2 * nx + 2 * nz);
    matRes << sqrt(sig.wgt(1)) * (Z.block(0, 1, nz, nx) - zm.replicate(1, nx)), sqrt(sig.wgt(nx + 1)) * (Z.block(0, nx + 1, nz, nx) - zm.replicate(1, nx)), sqrt(sig.wgt(2 * nx + 1)) * (Z.block(0, 2 * nx + 1, nz, 2 * nz) - zm.replicate(1, 2 * nz));
    HouseholderQR<MatrixXd> qr(matRes.transpose());

    MatrixXd matS2 = qr.matrixQR().triangularView<Upper>();
    MatrixXd matS = matS2.block(0, 0, nz, nz).transpose();
    MatrixXd matSzz = cholupdate(matS, Z.col(0) - zm, sig.wgt(0));

    // Pzz = Z * sig.wgt.asDiagonal() * Z.transpose() - zm * zm.transpose();         // Eq. B5
    Pzx = Z * sig.wgt.asDiagonal() * sig.state.transpose() - zm * xm.transpose(); // Eq. B6

    // K = Pzz.llt().solve(Pzx).transpose();
    // Kalman Gain
    MatrixXd matK = matSzz.transpose().householderQr().solve(matSzz.householderQr().solve(Pzx));
    K = matK.transpose();
    cout << "Kalman gain: \n"
         << K << endl;

    Dist distXu;
    distXu.n = nx;
    distXu.mean = xm + K * (z - zm); // Eq. B7
    MatrixXd matU = K * matSzz;
    distXu.covL = cholupdate(distx.back().covL, matU, -1.0);
    distXu.cov = distXu.covL * distXu.covL.transpose();
    Xu = sig.state - K * (Z.colwise() - zm);
    MatrixXd Xstd = distx.back().covL.triangularView<Lower>().solve(Xu.colwise() - Xu * sig.wgt); // standardised states at the sigma points,  covariance, A7/B10
    distXu.skew = Xstd.array().pow(3).matrix() * sig.wgt;                                         // skewness of the standardised state, Eq. A8/B11
    distXu.kurt = Xstd.array().pow(4).matrix() * sig.wgt;                                         // kurtosis of the standardised state, Eq. A9/B12

    cout << "updated mean: \n"
         << distXu.mean << endl;
    cout << "updated covariance: \n"
         << distXu.cov << endl;
    cout << "updated covariance lower triangle: \n"
         << distXu.covL << endl;
    cout << "updated skewness: \n"
         << distXu.skew << endl;
    cout << "updated kurtosis: \n"
         << distXu.kurt << endl;
    // // doesn't use Eq. 8 for covariance. instead, Pu = Pp - K*Pzz*K^T is used.
    // Xu = sig.state - K * (Z.colwise() - zm);

    // // // test Eq. 8
    // // MatrixXd Pxx = distx.back().cov - Pzz.llt().solve(Pzx).transpose() * Pzx;
    // // MatrixXd PxxL = Pxx.llt().matrixL();
    // // MatrixXd Xstd = PxxL.triangularView<Lower>().solve(Xu.colwise() - Xu * sig.wgt); // standardised states at the sigma points,  covariance, A7/B10
    // // cout << "skewness 1:\t" << Xstd.array().pow(3).matrix() * sig.wgt << endl
    // //      << "kurtosis 1:\t" << Xstd.array().pow(4).matrix() * sig.wgt << endl;

    // Dist distXu(Xu, sig.wgt);
    // cout << "skewness 2:\t" << distXu.skew << endl
    //      << "kurtosis 2:\t" << distXu.kurt << endl;
    // distXu.mean = xm + K * (z - zm); // Eq. B7

    distx.back() = distXu;
}

// Run filter for sequence of measurements
void SRHOUSE::run(const VectorXd &tz, const MatrixXd &Z)
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

SRHOUSE::SRHOUSE(
    const dyn_model &f_,
    const meas_model &h_,
    int nz_,
    double t0,
    double dtMax_,
    const Dist &distx0,
    const Dist &distw_,
    const Dist &distv_,
    double delta_) : HOUSE(f_,
                           h_,
                           nz_,
                           t0,
                           dtMax_,
                           distx0,
                           distw_,
                           distv_,
                           delta_),
                     f(f_),
                     h(h_),
                     nx(distx0.n),
                     nz(nz_),
                     nw(distw_.n),
                     nv(distv_.n),
                     distw(distw_),
                     distv(distv_),
                     delta(delta_),
                     dtMax(dtMax_)
{

    distx.push_back(distx0);
    t.push_back(t0);
}

// Reset filter
void SRHOUSE::reset(double t0, const Dist &distx0)
{

    t.clear();
    distx.clear();

    t.push_back(t0);
    distx.push_back(distx0);
}

// Save results
void SRHOUSE::save(const string &filename, string stateType)
{

    using namespace std;

    int steps = distx.size();

    MatrixXd table(steps, nx * (nx + 1) + 1);

    for (int k = 0; k < steps; k++)
    {

        table(k, 0) = t[k];
        if (stateType == "eci")
            table.row(k).segment(1, nx) = distx[k].mean;
        else if (stateType == "mee")
            table.row(k).segment(1, nx) = mee2eci(distx[k].mean, GM_Earth);

        table.row(k).tail(nx * nx) = distx[k].cov.reshaped(1, nx * nx);
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