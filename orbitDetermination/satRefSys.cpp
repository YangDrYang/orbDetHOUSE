
#include "satRefSys.hpp"

/* Class constants
 */
const double IERS::TT_TAI = +32.184; //  TT-TAI time difference [s]

const double IERS::GPS_TAI = -19; // GPS-TAI time difference [s]
const double IERS::TAI_GPS = +19; // TAI-GPS time difference [s]

const double IERS::TT_GPS = IERS::TT_TAI //  TT-GPS time difference [s]
                              - IERS::GPS_TAI;

/* Default values of Earth rotation parameters
 */
double IERS::UT1_TAI_ = 0; // UT1-TAI time difference [s]
double IERS::UTC_TAI_ = 0; // UTC-TAI time difference [s]
double IERS::x_pole_ = 0;  // Pole coordinate [rad]
double IERS::y_pole_ = 0;  // Pole coordinate [rad]
double IERS::lod_ = 0;     // length of day [s/d]

/* Setting of IERS Earth rotation parameters
 * (UT1-UTC [s], UTC-TAI [s], x [radian], y [radian])
 */
void IERS::Set(
    double UT1_UTC,
    double UTC_TAI,
    double x_pole,
    double y_pole,
    double lod)
{
    UT1_TAI_ = UT1_UTC + UTC_TAI;
    UTC_TAI_ = UTC_TAI;
    x_pole_ = x_pole;
    y_pole_ = y_pole;
    lod_ = lod;
};

/* Time differences [s]
 */
double IERS::UTC_TAI(double Mjd_UTC) { return UTC_TAI_; };
double IERS::UT1_TAI(double Mjd_UTC) { return UT1_TAI_; };

double IERS::UTC_GPS(double Mjd_UTC) { return UTC_TAI(Mjd_UTC) - GPS_TAI; };
double IERS::UT1_GPS(double Mjd_UTC) { return UT1_TAI(Mjd_UTC) - GPS_TAI; };

double IERS::TT_UTC(double Mjd_UTC) { return TT_TAI - UTC_TAI(Mjd_UTC); };
double IERS::GPS_UTC(double Mjd_UTC) { return GPS_TAI - UTC_TAI(Mjd_UTC); };
double IERS::UT1_UTC(double Mjd_UTC) { return UT1_TAI(Mjd_UTC) - UTC_TAI(Mjd_UTC); };

/* Pole coordinate [rad]
 */
double IERS::x_pole(double Mjd_UTC) { return x_pole_; };

/* Pole coordinate [rad]
 */
double IERS::y_pole(double Mjd_UTC) { return y_pole_; };

/* Length of day [s/d]
 */
double IERS::lod(double Mjd_UTC) { return lod_; };

/* Elementary rotations
 */
Matrix3d R_x(
    double Angle) ///< Angle in radian		//todo aaron delete
{
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix3d U = Matrix3d::Zero();
    U(0, 0) = 1.0;
    U(0, 1) = 0.0;
    U(0, 2) = 0.0;
    U(1, 0) = 0.0;
    U(1, 1) = +C;
    U(1, 2) = +S;
    U(2, 0) = 0.0;
    U(2, 1) = -S;
    U(2, 2) = +C;
    return U;
}

Matrix3d R_y(
    double Angle) ///< Angle in radian
{
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix3d U = Matrix3d::Zero();
    U(0, 0) = +C;
    U(0, 1) = 0.0;
    U(0, 2) = -S;
    U(1, 0) = 0.0;
    U(1, 1) = 1.0;
    U(1, 2) = 0.0;
    U(2, 0) = +S;
    U(2, 1) = 0.0;
    U(2, 2) = +C;
    return U;
}

Matrix3d R_z(
    double Angle) ///< Angle in radian
{
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix3d U = Matrix3d::Zero();
    U(0, 0) = +C;
    U(0, 1) = +S;
    U(0, 2) = 0.0;
    U(1, 0) = -S;
    U(1, 1) = +C;
    U(1, 2) = 0.0;
    U(2, 0) = 0.0;
    U(2, 1) = 0.0;
    U(2, 2) = 1.0;
    return U;
}