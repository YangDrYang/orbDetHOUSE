
#ifndef __FORCE_MODELS_HPP__
#define __FORCE_MODELS_HPP__

#include <map>

using std::map;

#include "atmosphereModel.hpp"
#include "configDefault.hpp"
#include "auxillaryData.hpp"
#include "coordTrans.hpp"
#include "satRefSys.hpp"
#include "gravity.hpp"
#include "jplEph.hpp"

Vector3d calculateAcceleration(
	double mjdTT,
	Vector3d &rSat,
	Vector3d &vSat,
	Matrix3d &mECI2ECEF,
	MatrixXd &cnm,
	MatrixXd &snm,
	double n_max,
	double m_max);

Vector3d earthCentralBodyAcc(
	Vector3d &rSat,		 ///< Inertial position of satellite (m)
	Vector3d &vSat,		 ///< Inertial velocity of satellite (m/s)
	Matrix3d &mECI2ECEF, ///< Transformation matrix from ECI coordinate to ECEF
	MatrixXd &cnm,		 ///< Earth model coefficients
	MatrixXd &snm,		 ///< Earth model coefficients
	double n_max,		 ///< Maximum Earth gravity degree
	double m_max);		 ///< Maximum Earth gravity order

Vector3d thirdBodyAttractionAcc();

struct ThirdBodyAttraction
{
	JPLEphemeris *mJPLEphemeris;

	ThirdBodyAttraction();

	ThirdBodyAttraction(
		E_SolarSysStarType thirdBodyStarType,
		E_SolarSysStarType centerStarType,
		void *pJPLEph);

	void setThirdBodyStar(
		E_SolarSysStarType thirdBodyStarType);

	void setCenterStar(
		E_SolarSysStarType centerStarType);

	E_SolarSysStarType getThirdBodyStar();
	E_SolarSysStarType getCenterStar();

	/*
	 *	Calculation of Three-body Gravitational Acceleration of Point Mass Model
	 */
	Vector3d accelPointMassGravity(
		const double mjdTT,	  ///< Terrestrial Time (Modified Julian Date)
		const Vector3d &pos); ///< Satellite Position Vector in the inertial system

protected:
	double mGM;
	E_SolarSysStarType mThirdBodyStarType;
	E_SolarSysStarType mCenterStarType;
};

Vector3d relativityEffectsAcc();

struct SRPOpt
{
	bool flagDirectSRP = 1;
	bool flagThermalEmission = 0;
	bool flagEarthAlbedo = 0;
	bool flagInfraredRadiation = 0;
};

struct SRPPara
{
	double satMass = 0; ///< Satellite mass, unit: kg
	double srpArea = 0; ///< SRP projected area, unit: m^2
	double srpCoef = 0; ///< SRP coefficient (cannonball model)
	E_SRPModels srpMdlName = E_SRPModels::CANNONBALL;
};

struct DragPara
{
	double satMass = 0;	 ///< Satellite mass, unit: kg
	double dragArea = 0; ///< SRP projected area, unit: m^2
	double dragCoef = 0; ///< SRP coefficient (cannonball model)
};

/*
 * Class type: The class of solar radiation pressure
 */
struct SolarRadPressure
{
	SolarRadPressure(
		SRPOpt optSRP,
		SRPPara paraSRP,
		void *pJPLEph);

	SRPOpt mSRPOpt;
	SRPPara mSRPPara;

	/*
	 * Calculation of direct solar radiation pressure acceleration
	 */
	Vector3d directSolarRadiationAcc(
		const double mjdTT,	   ///< Terrestrial time (modified Julian date)
		const Vector3d &rSat); ///< Satellite position vector, unit: m, m/s

	/*
	 * Calculation of indirect solar radiation pressure acceleration
	 */
	Vector3d indirectSolarRadiationAcc(
		const double mjdTT,	   ///< Terrestrial time (modified Julian date)
		const Vector3d &rSat); ///< Satellite position vector, unit: m, m/s

	JPLEphemeris *mJPLEphemeris;
};

Vector3d empiricalAcc();
Vector3d manoeuvreAcc();

struct PropOpt
{
	E_SolarSysStarType centerStarType;
	bool flagEarthGravity;
	EarthGravMdlOpt optEarthGravMdl;
	EarthTidesOpt optTides;
	bool flagRelativityEffect;
	bool flagThirdBodyGravity;
	ThirdBodyGravityOpt optThirdBodyGravity;
	bool flagSolarRadiationPressure;
	SRPOpt optSRP;
	SRPPara paraSRP;
	bool flagAtmosphericDrag;
	DragPara paraDrag;
	bool flagAntennaThrust;
	bool flagEmpiricalAcceleration;
	bool flagSatelliteManoeuvre;
	string odeInteg; ///< Integrator
};

struct Propagator
{
	Propagator();

	PropOpt propOpt;

	double mMJDUTC0;  ///< UTC modified Julian date at the starting epoch
	double mMJDUTC;	  ///< UTC modified Julian date
	IERS *mIERS;	  ///< Instance of IERS class
	Vector3d mRSat;	  ///< Inertial position of satellite (m)
	Vector3d mVSat;	  ///< Inertial velocity of satellite (m/s)
	erp_t *mERP;	  ///< Earth rotational parameters
	double mERPv[5];  ///< Earth rotational parameters (xp, yp, ut1_utc, lod, leapsec) for the current epoch
	EGMCoef mEGMCoef; ///< Earth model coefficients
	void *mpJPLEph;	  ///< Pointer to JPL ephemerides

	void setPropOption(
		ForceModels forceMdl);

	void printPropOption();

	void initPropagator(
		VectorXd rSatECI,
		double mjdUTC = 0,
		double leapSec = 0,
		erp_t *erpSrc = nullptr,
		EGMCoef egmCoe = {},
		void *pJPLEph = nullptr);

	void updPropagator(
		double mjdUTC); ///< mjdUTC at the current moment

	void updPropagator(
		double mjdUTC,	///< mjdUTC at the current moment
		double leapSec, ///< update the leap second if needed
		erp_t *erpSrc); ///< update the erp data source if needed

	Vector3d calculateAcceleration(
		const Vector3d &rSat,		///< Inertial position of satellite (m)
		const Vector3d &vSat,		///< Inertial velocity of satellite (m/s)
		const Matrix3d &mECI2ECEF); ///< Transformation matrix from ECI coordinate to ECEF

protected:
	GravityModel *pmGravityModel;
	ThirdBodyAttraction *pmThirdBodyGrva;
	SolarRadPressure *pmSolarRadPressure;
};

extern Propagator orbitProp;

#endif
