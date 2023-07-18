#include "coordTrans.hpp"

void arrayMat2eigenMat(
	const double arrayMat[3][3], ///< c++ 2D array
	Matrix3d &eigenMat)			 ///< eigen 2d matrix
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			eigenMat(i, j) = arrayMat[i][j];
		}
	}
}

/** eci to ecef transformation matrix using the sofa library
 */
void eci2ecef_sofa(
	const double mjdUTC, ///< UTC, modified Julian date
	IERS &iersIns,		 ///< IERS instance
	Matrix3d &U,		 ///< eci to ecef transformation matrix (3 x 3)
	Matrix3d &dU)		 ///< eci to ecef transformation matrix (3 x 3)
{

	double mjdTT = mjdUTC + iersIns.TT_UTC(mjdUTC) / 86400;
	double mjdUT1 = mjdUTC + iersIns.UT1_UTC(mjdUTC) / 86400;

	/* Form bias-precession-nutation matrix */

	double arrNPB[3][3] = {};
	Matrix3d matNPB = Matrix3d::Zero();
	iauPnm06a(JD2MJD, mjdTT, arrNPB);
	arrayMat2eigenMat(arrNPB, matNPB);

	/* Form Earth rotation matrix */
	double arrTheta[3][3] = {};
	Matrix3d matTheta = Matrix3d::Zero();
	arrTheta[0][0] = 1;
	arrTheta[1][1] = 1;
	arrTheta[2][2] = 1;
	iauRz(iauGst06(JD2MJD, mjdUT1, JD2MJD, mjdTT, arrNPB), arrTheta);
	arrayMat2eigenMat(arrTheta, matTheta);

	/* Polar motion matrix (TIRS->ITRS, IERS 2003) */
	double xp = iersIns.x_pole(mjdUTC);
	double yp = iersIns.y_pole(mjdUTC);
	double arrPi[3][3] = {};
	Matrix3d matPi = Matrix3d::Zero();
	iauPom00(xp, yp, iauSp00(JD2MJD, mjdTT), arrPi);
	arrayMat2eigenMat(arrPi, matPi);

	/* ICRS to ITRS transformation matrix and derivative */
	Matrix3d matS = Matrix3d::Zero();
	matS(0, 1) = 1;
	matS(1, 0) = -1; // Derivative of Earth rotation
	double lod = iersIns.lod(mjdUTC);
	double Omega = OMGE - 0.843994809 * 1e-9 * lod; // IERS, [rad/s]
	// double Omega2  = 7292115.8553e-11 + 4.3e-15 * ( (mjdUTC - mjdJ2000) / 36525 ); // [rad/s]
	Matrix3d matdTheta = Omega * matS * matTheta; // matrix [1/s]

	U = matPi * matTheta * matNPB;
	dU = matPi * matdTheta * matNPB;
}

/** vector transformed from eci to ecef using the sofa library
 */
void eci2ecefVec_sofa(
	const double mjdUTC,  ///< UTC, modified Julian date
	IERS &iersIns,		  ///< IERS instance
	VectorXd &rvSat_eci,  ///< satellite state vector in eci
	VectorXd &rvSat_ecef) ///< satellite state vector in ecef
{

	Matrix3d matECI2ECEF = Matrix3d::Zero();
	Matrix3d matdECI2ECEF = Matrix3d::Zero();
	eci2ecef_sofa(mjdUTC, iersIns, matECI2ECEF, matdECI2ECEF);

	rvSat_ecef.head(3) = matECI2ECEF * rvSat_eci.head(3);
	rvSat_ecef.tail(3) = matECI2ECEF * rvSat_eci.tail(3) + matdECI2ECEF * rvSat_eci.head(3);
}

/** vector transformed from ecef to eci using the sofa library
 */
void ecef2eciVec_sofa(
	const double mjdUTC,  ///< UTC, modified Julian date
	IERS &iersIns,		  ///< IERS instance
	VectorXd &rvSat_ecef, ///< satellite state vector in ecef
	VectorXd &rvSat_eci)  ///< satellite state vector in eci
{

	Matrix3d matECI2ECEF = Matrix3d::Zero();
	Matrix3d matdECI2ECEF = Matrix3d::Zero();
	eci2ecef_sofa(mjdUTC, iersIns, matECI2ECEF, matdECI2ECEF);

	rvSat_eci.head(3) = matECI2ECEF.transpose() * rvSat_ecef.head(3);
	rvSat_eci.tail(3) = matECI2ECEF.transpose() * rvSat_ecef.tail(3) + matdECI2ECEF.transpose() * rvSat_ecef.head(3);
}

VectorXd mee2coe(const VectorXd &mee)
{
	double p = mee(0);
	double f = mee(1);
	double g = mee(2);
	double h = mee(3);
	double k = mee(4);
	double L = mee(5);

	double e = sqrt(f * f + g * g);
	double a = p / (1 - e * e);
	double w_ohm = atan2(g / e, f / e);
	double tani_2 = sqrt(h * h + k * k);
	double i = atan(tani_2) * 2;
	double ohm = atan2(k / tan(i / 2), h / tan(i / 2));
	double w = w_ohm - ohm;
	double nu = L - w_ohm;

	VectorXd coe(6);
	coe << a, e, i, w, ohm, nu;
	return coe;
}

VectorXd coe2eci(const VectorXd &coeEl, double mu)
{
	double a = coeEl(0);
	double e = coeEl(1);
	double i = coeEl(2);
	double w = coeEl(3);
	double ohm = coeEl(4);
	double nu = coeEl(5);

	if (e == 1.0)
	{
		throw runtime_error("e is 1");
	}

	double p = a * (1 - e * e);
	Vector3d rPQW(p * cos(nu) / (1 + e * cos(nu)),
				  p * sin(nu) / (1 + e * cos(nu)),
				  0.0);
	Vector3d vPQW(-sqrt(mu / p) * sin(nu),
				  sqrt(mu / p) * (e + cos(nu)),
				  0.0);
	Matrix3d tempMat;
	tempMat << cos(ohm) * cos(w) - sin(ohm) * sin(w) * cos(i),
		-cos(ohm) * sin(w) - sin(ohm) * cos(w) * cos(i),
		sin(ohm) * sin(i),
		sin(ohm) * cos(w) + cos(ohm) * sin(w) * cos(i),
		-sin(ohm) * sin(w) + cos(ohm) * cos(w) * cos(i),
		-cos(ohm) * sin(i),
		sin(w) * sin(i),
		cos(w) * sin(i),
		cos(i);
	Vector3d rIjk = tempMat * rPQW;
	Vector3d vIjk = tempMat * vPQW;

	VectorXd eciEl(6);
	eciEl << rIjk, vIjk;
	return eciEl;
}

VectorXd coe2mee(const VectorXd &coeEl)
{
	// coeEl = [a, e, i, w, Ohm, nu];
	double a = coeEl(0);
	double e = coeEl(1);
	double i = coeEl(2);
	double w = coeEl(3);
	double Ohm = coeEl(4);
	double nu = coeEl(5);

	if (abs(e - 1.0) < 1e-12)
	{
		throw runtime_error("The orbit is circular, the code will not work");
	}

	double p = a * (1 - e * e);
	double f = e * cos(w + Ohm);
	double g = e * sin(w + Ohm);
	double h = tan(i / 2) * cos(Ohm);
	double k = tan(i / 2) * sin(Ohm);
	double L = Ohm + w + nu;

	VectorXd meeEl(6);
	meeEl << p, f, g, h, k, L;

	return meeEl;
}

VectorXd eci2coe(const VectorXd &eciEl, double mu)
{
	// Format of ECI elements
	Vector3d rVec = eciEl.head(3);
	Vector3d vVec = eciEl.tail(3);

	// Step 1:
	Vector3d hVec = rVec.cross(vVec);
	double h = hVec.norm();

	// Step 2:
	Vector3d nVec = Vector3d::UnitZ().cross(hVec);
	Vector3d eVec = ((vVec.dot(vVec) - (mu / rVec.norm())) * rVec - rVec.dot(vVec) * vVec) / mu;
	double e = eVec.norm();

	// Step 3:
	double zeta = (vVec.dot(vVec) / 2) - (mu / rVec.norm());
	double a;
	if (abs(e - 1.0) > 1e-12)
	{
		a = -mu / (2 * zeta);
	}
	else
	{
		// cout << "THE ORBIT IS CIRCULAR" << endl;
		throw runtime_error("THE ORBIT IS parabolic");
		a = INFINITY;
	}

	// Step 4:
	double i = acos(hVec(2) / h);
	double Ohm = acos(nVec(0) / nVec.norm());
	if (nVec(1) < 0)
	{
		Ohm = 2 * M_PI - Ohm;
	}

	// Step 5:
	double w = acos(nVec.dot(eVec) / (nVec.norm() * eVec.norm()));
	if (eVec(2) < 0)
	{
		w = 2 * M_PI - w;
	}

	// Step 6:
	double nu = acos(eVec.dot(rVec) / (eVec.norm() * rVec.norm()));
	if (rVec.dot(vVec) < 0)
	{
		nu = 2 * M_PI - nu;
	}

	VectorXd coeEl(6);
	coeEl << a, e, i, w, Ohm, nu;

	return coeEl;
}

VectorXd eci2mee(const VectorXd &eciEl, double mu)
{
	return coe2mee(eci2coe(eciEl, mu));
}