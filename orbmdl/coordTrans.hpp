
#ifndef ___COORDTRANS_HPP__
#define ___COORDTRANS_HPP__

#include "satRefSys.hpp"
#include "constants.hpp"
#include "sofam.hpp"
#include "sofa.hpp"
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
/* coordinates transformation ------------------------------------------------*/

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

#endif
