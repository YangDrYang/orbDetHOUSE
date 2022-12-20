#include "common.hpp"
#include "sofam.hpp"
#include "sofa.hpp"

void arrayMat2eigenMat(
	const double	arrayMat[3][3],	///< c++ 2D array
	Matrix3d&		eigenMat)	///< eigen 2d matrix
{
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			eigenMat(i, j) = arrayMat[i][j];
		}
	}
}

/** eci to ecef transformation matrix using the sofa library
*/
void eci2ecef_sofa(
	const double	mjdUTC, ///< UTC, modified Julian date 
	IERS&			iersIns,///< IERS instance
	Matrix3d&		U,		///< eci to ecef transformation matrix (3 x 3)
	Matrix3d&		dU)		///< eci to ecef transformation matrix (3 x 3)
{

	double mjdTT	= mjdUTC + iersIns.TT_UTC	(mjdUTC) / 86400;
	double mjdUT1	= mjdUTC + iersIns.UT1_UTC	(mjdUTC) / 86400;

	/* Form bias-precession-nutation matrix */
	
	double arrNPB[3][3] = {0};
	Matrix3d matNPB = Matrix3d::Zero();
	iauPnm06a(JD2MJD, mjdTT, arrNPB);
	arrayMat2eigenMat(arrNPB, matNPB);

	/* Form Earth rotation matrix */
	double arrTheta[3][3] = {0};
	Matrix3d matTheta = Matrix3d::Zero();
	arrTheta[0][0] = 1; arrTheta[1][1] = 1; arrTheta[2][2] = 1;
	iauRz(iauGst06(JD2MJD, mjdUT1, JD2MJD, mjdTT, arrNPB), arrTheta);
	arrayMat2eigenMat(arrTheta, matTheta);
	
	/* Polar motion matrix (TIRS->ITRS, IERS 2003) */
	double xp = iersIns.x_pole(mjdUTC);
	double yp = iersIns.y_pole(mjdUTC);
	double arrPi[3][3] = {0};
	Matrix3d matPi = Matrix3d::Zero();
	iauPom00(xp, yp, iauSp00(JD2MJD, mjdTT), arrPi);
	arrayMat2eigenMat(arrPi, matPi);

	/* ICRS to ITRS transformation matrix and derivative */	
	Matrix3d matS = Matrix3d::Zero();
	matS(0, 1) = 1; matS(1, 0) = -1;       //Derivative of Earth rotation
	double lod = iersIns.lod(mjdUTC);
	double Omega = OMGE - 0.843994809*1e-9 * lod;  //IERS, [rad/s]
	double Omega2  = 7292115.8553e-11 + 4.3e-15 * ( (mjdUTC - mjdJ2000) / 36525 ); // [rad/s]
	Matrix3d matdTheta = Omega * matS * matTheta;   // matrix [1/s]

	U = matPi * matTheta * matNPB;
	dU = matPi * matdTheta * matNPB;
}

/** vector transformed from eci to ecef using the sofa library
*/
void eci2ecefVec_sofa(
	const double	mjdUTC,		///< UTC, modified Julian date
	IERS&			iersIns,	///< IERS instance	
	Vector6d&		rvSat_eci,	///< satellite state vector in eci
	Vector6d&		rvSat_ecef) ///< satellite state vector in ecef
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
	const double	mjdUTC,		///< UTC, modified Julian date
	IERS&			iersIns,	///< IERS instance	
	Vector6d&		rvSat_ecef,	///< satellite state vector in ecef
	Vector6d&		rvSat_eci) ///< satellite state vector in eci
{

	Matrix3d matECI2ECEF = Matrix3d::Zero();
	Matrix3d matdECI2ECEF = Matrix3d::Zero();
	eci2ecef_sofa(mjdUTC, iersIns, matECI2ECEF, matdECI2ECEF);

	rvSat_eci.head(3) = matECI2ECEF.transpose() * rvSat_ecef.head(3);
	rvSat_eci.tail(3) = matECI2ECEF.transpose() * rvSat_ecef.tail(3) + matdECI2ECEF.transpose() * rvSat_ecef.head(3);
}

