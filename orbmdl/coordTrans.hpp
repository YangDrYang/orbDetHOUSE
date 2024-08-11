
#ifndef ___COORDTRANS_HPP__
#define ___COORDTRANS_HPP__

#include "satRefSys.hpp"
#include "constants.hpp"
#include "sofam.hpp"
#include "sofa.hpp"
#include <Eigen/Dense>
#include <iostream>
using namespace Eigen;
using namespace std;
/* coordinates transformation ------------------------------------------------*/

#define MAXLEAPS 18
const double leaps[MAXLEAPS + 1][7] =
	{
		/* leap seconds (y,m,d,h,m,s,utc-gpst) */
		{2017, 1, 1, 0, 0, 0, -18},
		{2015, 7, 1, 0, 0, 0, -17},
		{2012, 7, 1, 0, 0, 0, -16},
		{2009, 1, 1, 0, 0, 0, -15},
		{2006, 1, 1, 0, 0, 0, -14},
		{1999, 1, 1, 0, 0, 0, -13},
		{1997, 7, 1, 0, 0, 0, -12},
		{1996, 1, 1, 0, 0, 0, -11},
		{1994, 7, 1, 0, 0, 0, -10},
		{1993, 7, 1, 0, 0, 0, -9},
		{1992, 7, 1, 0, 0, 0, -8},
		{1991, 1, 1, 0, 0, 0, -7},
		{1990, 1, 1, 0, 0, 0, -6},
		{1988, 1, 1, 0, 0, 0, -5},
		{1985, 7, 1, 0, 0, 0, -4},
		{1983, 7, 1, 0, 0, 0, -3},
		{1982, 7, 1, 0, 0, 0, -2},
		{1981, 7, 1, 0, 0, 0, -1},
		{0}};

time_t convertMJD2Time_T(double mjd);
double getLeapSecond(time_t t);
void getIERS(double mjd, const erp_t erpt, IERS &iersInstance);

void eci2ecef_sofa(
	const double mjdUTC,
	IERS &iersIns,
	Matrix3d &U,
	Matrix3d &dU);

void eci2ecefVec_sofa(
	const double mjdUTC,
	IERS &iersIns,
	VectorXd &rvSat_eci,
	VectorXd &rvSat_ecef);

void ecef2eciVec_sofa(
	const double mjdUTC,
	IERS &iersIns,
	VectorXd &rvSat_ecef,
	VectorXd &rvSat_eci);

VectorXd mee2coe(const VectorXd &mee);
VectorXd coe2mee(const VectorXd &coe);
VectorXd coe2eci(const VectorXd &coe, double mu);
VectorXd eci2coe(const VectorXd &eci, double mu);
VectorXd eci2mee(const VectorXd &eci, double mu);
VectorXd mee2eci(const VectorXd &mee, double mu);
MatrixXd ric2eci(const MatrixXd &ric, const VectorXd &eci);

#endif
