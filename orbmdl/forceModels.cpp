
#include "forceModels.hpp"
#include "satRefSys.hpp"
#include "constants.hpp"
#include "jplEph.hpp"
#include "enum.h"

#include <iomanip>

// check that this is correct
#define PI2 M_PI * 2
#define PSOL 1367 / CLIGHT /* solar radiation pressure at 1 AU, [N/m^2] (1367 W/m^2); IERS 96 */
#define CLIGHT 299792458.0 /* speed of light (m/s) */

/*
 * Class type: The class of the third body attraction
 */
ThirdBodyAttraction::ThirdBodyAttraction()
{
	mGM = 0;
	mThirdBodyStarType = eEarth;
	mCenterStarType = eEarth;
}

ThirdBodyAttraction::ThirdBodyAttraction(
	E_SolarSysStarType thirdBodyStarType,
	E_SolarSysStarType centerStarType,
	void *pJPLEph)
{
	mCenterStarType = centerStarType;
	mThirdBodyStarType = thirdBodyStarType;

	switch (thirdBodyStarType)
	{
	case eMercury:
		mGM = GM_Mercury;
		break;
	case eVenus:
		mGM = GM_Venus;
		break;
	case eEarth:
		mGM = GM_Earth;
		break;
	case eMars:
		mGM = GM_Mars;
		break;
	case eJupiter:
		mGM = GM_Jupiter;
		break;
	case eSaturn:
		mGM = GM_Saturn;
		break;
	case eUranus:
		mGM = GM_Uranus;
		break;
	case eNeptune:
		mGM = GM_Neptune;
		break;
	case ePluto:
		mGM = GM_Pluto;
		break;
	case eMoon:
		mGM = GM_Moon;
		break;
	case eSun:
		mGM = GM_Sun;
		break;
	default:
		std::cout << "ThirdBodyAttraction: E_SolarSysStarType Unsupport " << endl;
	}

	mJPLEphemeris = new JPLEphemeris(pJPLEph);
}

void ThirdBodyAttraction::setThirdBodyStar(
	E_SolarSysStarType thirdBodyStarType)
{
	mThirdBodyStarType = thirdBodyStarType;
	switch (thirdBodyStarType)
	{
	case eMercury:
		mGM = GM_Mercury;
		break;
	case eVenus:
		mGM = GM_Venus;
		break;
	case eEarth:
		mGM = GM_Earth;
		break;
	case eMars:
		mGM = GM_Mars;
		break;
	case eJupiter:
		mGM = GM_Jupiter;
		break;
	case eSaturn:
		mGM = GM_Saturn;
		break;
	case eUranus:
		mGM = GM_Uranus;
		break;
	case eNeptune:
		mGM = GM_Neptune;
		break;
	case ePluto:
		mGM = GM_Pluto;
		break;
	case eMoon:
		mGM = GM_Moon;
		break;
	case eSun:
		mGM = GM_Sun;
		break;
	default:
		std::cout << "ThirdBodyAttraction: E_SolarSysStarType Unsupport " << endl;
	}
}

void ThirdBodyAttraction::setCenterStar(
	E_SolarSysStarType centerStarType)
{
	mCenterStarType = centerStarType;
}

E_SolarSysStarType ThirdBodyAttraction::getThirdBodyStar()
{
	return mThirdBodyStarType;
}

E_SolarSysStarType ThirdBodyAttraction::getCenterStar()
{
	return mCenterStarType;
}

Vector3d ThirdBodyAttraction::accelPointMassGravity(
	const double mjdTT,
	const Vector3d &pos)
{
	Vector3d relativePos;
	Vector3d bodyPos;

	/*  Relative position vector of satellite w.r.t. point mass
	 */
	mJPLEphemeris->getJPLEphemeris(mjdTT + JD2MJD, mThirdBodyStarType, mCenterStarType, bodyPos);
	relativePos = pos - bodyPos;

	/* Acceleration
	 */
	return -mGM * (relativePos / pow(relativePos.norm(), 3) + bodyPos / pow(bodyPos.norm(), 3));
}

Vector3d relativityEffectsAcc(
	const Vector3d &rSat,
	const Vector3d &vSat)
{
	// Relativistic Effects

	double rSatNorm = rSat.norm();
	double vSatNorm = vSat.norm();

	Vector3d acc = EARTH_G_CONST / (pow(CLIGHT, 2) * pow(rSatNorm, 3)) * ((4 * EARTH_G_CONST / rSatNorm - pow(vSatNorm, 2)) * rSat + 4 * rSat.dot(vSat) * vSat);
	return acc;
}

/* Frac: Gives the fractional part of a number
 */
double Frac(double x)
{
	return x - floor(x);
}

/* Computes the Sun's geocentric position using a low precision analytical series
 */
Vector3d Sun(
	double mjdTT) ///< Terrestrial time: modified Julian date
{
	/* Constants
	 */
	const double eps = 23.43929111 * D2R;		   // Obliquity of J2000 ecliptic
	const double T = (mjdTT - mjdJ2000) / 36525.0; // Julian cent. since J2000

	/* Mean anomaly, ecliptic longitude and radius
	 */

	// PI2, this this meant to be PI * 2,
	double M = PI2 * Frac(0.9931267 + 99.9973583 * T);												 // [rad]
	double L = PI2 * Frac(0.7859444 + M / PI2 + (6892.0 * sin(M) + 72.0 * sin(2.0 * M)) / 1296.0e3); // [rad]
	double r = 149.619e9 - 2.499e9 * cos(M) - 0.021e9 * cos(2 * M);									 // [m]

	/* Solar position vector [m] with respect to the mean equator and equinox of J2000 (EME2000, ICRF)
	 */
	return R_x(-eps) * Vector3d(r * cos(L), r * sin(L), 0.0);
}

/* Computes the fractional illumination of a spacecraft in the vicinity of the Earth assuming a cylindrical shadow model
 *
 */
double Illumination(
	const Vector3d &rSat, ///< Spacecraft position vector [m]
	const Vector3d &rSun) ///< Sun position vector [m]
{
	Vector3d eSun = rSun / rSun.norm(); // Sun direction unit vector
	double s = rSat.dot(eSun);			// Projection of s/c position

	/* Illumination factor:
	 *  1: Spacecraft fully illuminated by the Sun
	 *  0: Spacecraft in Earth shadow
	 */
	return ((s > 0 || (rSat - s * eSun).norm() > RE_WGS84) ? 1 : 0);
}

/* Class type: The class of solar radiation pressure
 *
 */
SolarRadPressure::SolarRadPressure(
	SRPOpt optSRP,
	SRPPara paraSRP,
	void *pJPLEph)
{
	mSRPOpt = optSRP;
	mSRPPara = paraSRP;
	mJPLEphemeris = new JPLEphemeris(pJPLEph);
}

Vector3d SolarRadPressure::directSolarRadiationAcc(
	double mjdTT,		  ///< Terrestrial time (modified Julian date)
	const Vector3d &rSat) ///< Satellite position vector, unit: m, m/s
{

	/* Relative position vector of spacecraft w.r.t. Sun
	 */
	Vector3d rSun;
	mJPLEphemeris->getJPLEphemeris(mjdTT + JD2MJD, eSun, eEarth, rSun);

	// cout << "Calculated sun position: " << setw(14) << mjdTT << setw(14) << rSun.transpose() << endl;

	Vector3d rDis = rSat - rSun;
	switch (mSRPPara.srpMdlName)
	{
	case E_SRPModels::CANNONBALL:
		return Illumination(rSat, rSun) * mSRPPara.srpCoef * (mSRPPara.srpArea / mSRPPara.satMass) * PSOL * (AU * AU) * rDis / pow(rDis.norm(), 3);
	case E_SRPModels::BOXWING:
		return Vector3d::Zero();
	case E_SRPModels::ECOM:
		return Vector3d::Zero();
	case E_SRPModels::ECOM2:
		return Vector3d::Zero();
	default:
		return Vector3d::Zero();
	}
}

Vector3d SolarRadPressure::indirectSolarRadiationAcc(
	double mjdTT,		  ///< Terrestrial time (modified Julian date)
	const Vector3d &rSat) ///< Satellite position vector, unit: m, m/s
{
	return Vector3d::Zero();
}

Vector3d antennarThrustAcc()
{
	return Vector3d::Zero();
}

Vector3d empiricalAcc()
{
	return Vector3d::Zero();
}

Vector3d manoeuvreAcc()
{
	return Vector3d::Zero();
}

Propagator::Propagator()
{
	pmGravityModel = nullptr;
	pmThirdBodyGrva = nullptr;
	pmSolarRadPressure = nullptr;
}

/* Set options for force models in the propagator
 *
 */
void Propagator::setPropOption(
	ForceModels forceMdl)
{
	/* Options from yaml file */
	propOpt.centerStarType = eEarth;
	propOpt.flagEarthGravity = forceMdl.earth_gravity;
	propOpt.optEarthGravMdl.earthGravMdl = E_GravMdl::GGM03SModel; // to do: yaml input
	propOpt.optTides.earthModel = E_TidesMdl::ElasticMdl;		   // to do: yaml input
	propOpt.optTides.flagOceanTides = forceMdl.ocean_tide_loading;
	propOpt.optTides.flagSolidEarthTides = forceMdl.solid_earth_tide;
	propOpt.flagRelativityEffect = forceMdl.relativity_effect;

	propOpt.flagThirdBodyGravity = forceMdl.third_body_attraction;
	propOpt.optThirdBodyGravity.flagSunGrav = forceMdl.third_body_sun;
	propOpt.optThirdBodyGravity.flagMoonGrav = forceMdl.third_body_moon;
	propOpt.optThirdBodyGravity.flagJupiterGrav = forceMdl.third_body_planet;
	propOpt.optThirdBodyGravity.flagVenusGrav = forceMdl.third_body_planet;
	propOpt.optThirdBodyGravity.flagMarsGrav = forceMdl.third_body_planet;
	propOpt.optThirdBodyGravity.flagSaturnGrav = forceMdl.third_body_planet;
	propOpt.optThirdBodyGravity.flagUranusGrav = forceMdl.third_body_planet;
	propOpt.optThirdBodyGravity.flagNeptuneGrav = forceMdl.third_body_planet;
	propOpt.optThirdBodyGravity.flagPlutoGrav = forceMdl.third_body_planet;

	propOpt.flagSolarRadiationPressure = forceMdl.solar_radiation_pressure;
	propOpt.optSRP.flagDirectSRP = forceMdl.solar_radiation_pressure;
	propOpt.optSRP.flagThermalEmission = forceMdl.thermal_emission;
	propOpt.optSRP.flagEarthAlbedo = forceMdl.earth_albedo;
	propOpt.optSRP.flagInfraredRadiation = forceMdl.infrared_radiation;

	propOpt.flagAntennaThrust = forceMdl.antenna_thrust;
	propOpt.flagEmpiricalAcceleration = forceMdl.empirical_acceleration;
	propOpt.flagSatelliteManoeuvre = forceMdl.satellite_manoeuvre;

	propOpt.odeInteg = forceMdl.odeInteg;

	/* Parameters from yaml file */
	propOpt.optEarthGravMdl.earthGravAccDeg.mMax = forceMdl.egmAccDeg;
	propOpt.optEarthGravMdl.earthGravAccDeg.nMax = forceMdl.egmAccOrd;
	propOpt.optEarthGravMdl.earthGravSTMDeg.mMax = forceMdl.egmSTMDeg;
	propOpt.optEarthGravMdl.earthGravSTMDeg.nMax = forceMdl.egmSTMOrd;

	propOpt.paraSRP.srpMdlName = forceMdl.srpMdlName;
	propOpt.paraSRP.satMass = forceMdl.satMass;
	propOpt.paraSRP.srpArea = forceMdl.srpArea;
	propOpt.paraSRP.srpCoef = forceMdl.srpCoef;
}

/* initialise the propagator with time, state and necessary parameters
 *
 */
void Propagator::initPropagator(
	VectorXd rvSatECI,
	double mjdUTC,
	double leapSec,
	erp_t *erpSrc,
	EGMCoef egmCoe,
	void *pJPLEph)
{
	mRSat = rvSatECI.head(3);
	mVSat = rvSatECI.tail(3);
	mMJDUTC0 = mjdUTC;
	mMJDUTC = mjdUTC;
	mERP = erpSrc;
	mpJPLEph = pJPLEph;
	geterp_from_utc(mERP, leapSec, mMJDUTC, mERPv);

	pmGravityModel = new GravityModel(propOpt.optEarthGravMdl, propOpt.optTides, egmCoe, pJPLEph);

	if (propOpt.flagSolarRadiationPressure)
	{
		pmSolarRadPressure = new SolarRadPressure(propOpt.optSRP, propOpt.paraSRP, pJPLEph);
	}

	if (propOpt.flagThirdBodyGravity)
	{
		pmThirdBodyGrva = new ThirdBodyAttraction(propOpt.centerStarType, propOpt.centerStarType, pJPLEph);
	}
}

void Propagator::updPropagator(
	double mjdUTC)
{
	/* update time
	 */
	mMJDUTC = mjdUTC;
	/* update erp
	 */
	geterp_from_utc(mERP, mERPv[4], mMJDUTC, mERPv);
	double dUT1_UTC = mERPv[2];
	double dUTC_TAI = -(19 + mERPv[4]);
	double xp = mERPv[0];
	double yp = mERPv[1];
	double lod = mERPv[3];
	/* update iers
	 */
	mIERS->Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
}

void Propagator::updPropagator(
	double mjdUTC,
	double leapSec,
	erp_t *erpSrc)
{
	/* update time
	 */
	mMJDUTC = mjdUTC;
	/* update erp
	 */
	mERP = erpSrc;
	geterp_from_utc(mERP, leapSec, mMJDUTC, mERPv);
	double dUT1_UTC = mERPv[2];
	double dUTC_TAI = -(19 + leapSec);
	double xp = mERPv[0];
	double yp = mERPv[1];
	double lod = mERPv[3];
	/* update iers
	 */
	mIERS->Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
}

Vector3d Propagator::calculateAcceleration(
	const Vector3d &rSat,	   ///< Inertial position of satellite (m)
	const Vector3d &vSat,	   ///< Inertial velocity of satellite (m/s)
	const Matrix3d &mECI2ECEF) ///< Transformation matrix from ECI coordinate to ECEF
{
	Vector3d acc = Vector3d::Zero();

	/*calculate acceleration components
	 */
	Vector3d earthgravityAcc = Vector3d::Zero();
	bool bVarEq = false;
	if (propOpt.flagEarthGravity)
	{
		// Extremely expensive computation
		earthgravityAcc = pmGravityModel->centralBodyGravityAcc(mMJDUTC, mERPv, rSat, mECI2ECEF, bVarEq);
		// earthgravityAcc = -MU / (pow(rSat.norm(), 3)) * rSat;
		if (0)
		{
			cout << "Calculated accleration due to the Earth's central body gravity: " << setw(14) << mMJDUTC << setw(14) << earthgravityAcc.transpose() << endl;
		}
		acc += earthgravityAcc;
	}

	Vector3d thirdbodyAcc = Vector3d::Zero();
	if (propOpt.flagThirdBodyGravity)
	{
		double mjdTT = mIERS->TT_UTC(mMJDUTC) / 86400.0 + mMJDUTC;
		if (propOpt.optThirdBodyGravity.flagSunGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(eSun);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Sun's attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}

		if (propOpt.optThirdBodyGravity.flagMoonGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(eMoon);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Moon's attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}

		if (propOpt.optThirdBodyGravity.flagJupiterGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(eJupiter);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Jupiter's attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}

		if (propOpt.optThirdBodyGravity.flagVenusGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(eVenus);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Venus' attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}

		if (propOpt.optThirdBodyGravity.flagMarsGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(eMars);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Mars' attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}

		if (propOpt.optThirdBodyGravity.flagSaturnGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(eSaturn);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Saturn's attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}

		if (propOpt.optThirdBodyGravity.flagUranusGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(eUranus);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Uranus' attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}

		if (propOpt.optThirdBodyGravity.flagNeptuneGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(eNeptune);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Neptune' attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}

		if (propOpt.optThirdBodyGravity.flagPlutoGrav)
		{
			pmThirdBodyGrva->setThirdBodyStar(ePluto);
			thirdbodyAcc = pmThirdBodyGrva->accelPointMassGravity(mjdTT, rSat);
			acc += thirdbodyAcc;
			if (0)
			{
				cout << "Calculated acceleration due to Pluto' attraction: " << setw(14) << mMJDUTC << setw(14) << thirdbodyAcc.transpose() << endl;
			}
		}
	}

	Vector3d relativityAcc = Vector3d::Zero();
	if (propOpt.flagRelativityEffect)
	{

		relativityAcc = pmGravityModel->relativityEffectsAcc(rSat, vSat);
		acc += relativityAcc;

		if (0)
		{
			cout << "Calculated acceleration due to the relativity effect: " << setw(14) << mMJDUTC << setw(14) << relativityAcc.transpose() << endl;
		}
	}

	Vector3d directSRPAcc = Vector3d::Zero();
	if (propOpt.flagSolarRadiationPressure)
	{
		pmSolarRadPressure->mSRPPara = propOpt.paraSRP;
		directSRPAcc = pmSolarRadPressure->directSolarRadiationAcc(mMJDUTC, rSat);
		acc += directSRPAcc;

		if (1)
		{
			cout << "Calculated acceleration due to the direct solar radiation: " << setw(14) << mMJDUTC << setw(14) << directSRPAcc.transpose() << endl;
		}
	}
	return acc;
}

Vector3d Propagator::calculateAccelNGradient(
	const Vector3d &rSat,	   ///< Inertial position of satellite (m)
	const Vector3d &vSat,	   ///< Inertial velocity of satellite (m/s)
	Matrix3d &mGradient,	   ///< Gradient (G=da/dr) in the ICRF/J2000 system
	Vector3d &dadCr,		   ///< Partials of acceleration w.r.t. to the solar radiation coefficient
	const Matrix3d &mECI2ECEF) ///< Transformation matrix from ECI coordinate to ECEF
{

	bool bVarEq = true;
	Vector3d earthGravityAcc = pmGravityModel->centralBodyGravityAcc(mMJDUTC, mERPv, rSat, mECI2ECEF, bVarEq);
	Vector3d drSat;
	const double dinc = 1.0; // Position increment [m]
	Vector3d daSat;

	/* Gradient
	 */
	for (int i = 0; i < 3; i++)
	{
		drSat = Vector3d::Zero();
		drSat(i) = dinc;																								  // Set offset in i-th component of the position vector
		daSat = pmGravityModel->centralBodyGravityAcc(mMJDUTC, mERPv, rSat + drSat, mECI2ECEF, bVarEq) - earthGravityAcc; // Acceleration difference
		mGradient.col(i) = daSat / dinc;																				  // Derivative with respect to i-th component
	}

	/* Radiation pressure coefficient partials
	 */
	double mjdTT = mIERS->TT_UTC(mMJDUTC) / 86400.0 + mMJDUTC;
	pmSolarRadPressure->mSRPPara.srpCoef = 1.0;
	dadCr = pmSolarRadPressure->directSolarRadiationAcc(mjdTT, rSat);
}
