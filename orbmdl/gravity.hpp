#ifndef __GRAVITY_HPP__
#define __GRAVITY_HPP__

#include "configDefault.hpp"
#include "satRefSys.hpp"
#include "jplEph.hpp"
#include <string>
using std::string;

// ToDo: fix this, this is already defined in the config header file

/** Struct to save Earth gravity coefficients
 */
struct EarthGravityDeg
{
	int mMax = GRAVITY_DEG_M;
	int nMax = GRAVITY_DEG_N;
};

struct EGMCoef
{
	MatrixXd smn = MatrixXd::Zero(GRAVITY_DEG_N + 1, GRAVITY_DEG_N + 1);
	MatrixXd cmn = MatrixXd::Zero(GRAVITY_DEG_M + 1, GRAVITY_DEG_N + 1);
};

MatrixXd Legendre(
	int n,
	int m,
	double phi);

MatrixXd LegendreD(
	int m,		  ///< Maximum degree
	int n,		  ///< Maximum order
	MatrixXd pnm, ///< Normalised Legendre polinomial matrix
	double phi);  ///< Geocentric latitude in radian

struct EarthGravMdlOpt
{
	E_GravMdl earthGravMdl;
	EarthGravityDeg earthGravAccDeg;
	EarthGravityDeg earthGravSTMDeg;
};

struct EarthTidesOpt
{
	bool flagSolidEarthTides = 0;					///< Consider solid Earth tides or not
	bool flagOceanTides = 0;						///< Consider ocean tides or not
	E_TidesMdl earthModel = E_TidesMdl::ElasticMdl; ///< Elastic model or anelastic model
};

/** The gravity model of the central body
 */
struct GravityModel
{
public:
	// GravityModel();
	GravityModel(
		EarthGravMdlOpt gravMdlOpt,
		EarthTidesOpt tidesOpt,
		EGMCoef egmCoef,
		void *pJPLEph);
	// virtual ~GravityModel();

public:
	void setModelType(E_GravMdl modelType);

	E_GravMdl getModelType();

	// EarthGravityDeg mEarthGravAccDeg = {15, 15};
	EarthGravityDeg mEarthGravAccDeg;

	// EarthGravityDeg mEarthGravSTMDeg = {2, 2};
	EarthGravityDeg mEarthGravSTMDeg;

public:
	/* Computes the acceleration due to the central body gravity field of the central body
	 *
	 */
	Vector3d centralBodyGravityAcc(
		const double mjdUTC,
		const double *erpv,
		const Vector3d &rSat,	 ///< Satellite position vector in the inertial system
		const Matrix3d &mECI2BF, ///< Transformation matrix from ECI to central body fixed system
		bool bVarEq);			 ///< bVarEq = 1 if for the variational equation

	void solidEarthTidesCorrection(
		double mjdUTC,				 ///<
		EGMCoef &egmCoef,			 ///< Struct of Earth gravity coefficients
		IERS iersIns,				 ///< Instance of IERS class
		const Vector3d &vecRAESun,	 ///< Rho, azimuth and altitude information of Sun
		const Vector3d &vecRAEMoon); ///< Rho, azimuth and altitude information of Moon

	void oceanTidesCorrection(
		double mjdUTC,				 ///<
		EGMCoef &egmCoef,			 ///< Struct of Earth gravity coefficients
		const Vector3d &vecRAESun,	 ///< Rho, azimuth and altitude information of Sun
		const Vector3d &vecRAEMoon); ///< Rho, azimuth and altitude information of Moon

	Vector3d relativityEffectsAcc(
		const Vector3d &rSat,  ///< Inertial position of satellite (m)
		const Vector3d &vSat); ///< Inertial velocity of satellite (m/s)

	JPLEphemeris *mJPLEphemeris;

protected:
	E_GravMdl mGravModelType;
	EarthTidesOpt mEarthTidesOpt;
	double mGM;
	double mRadius;
	EGMCoef mEGMCoef;
};

struct ThirdBodyGravityOpt
{
	bool flagMercuryGrav = false;
	bool flagVenusGrav = false;
	bool flagEarthGrav = false;
	bool flagMarsGrav = false;
	bool flagJupiterGrav = false;
	bool flagSaturnGrav = false;
	bool flagUranusGrav = false;
	bool flagNeptuneGrav = false;
	bool flagPlutoGrav = false;
	bool flagMoonGrav = true;
	bool flagSunGrav = true;
};

extern EGMCoef egm;

#endif
