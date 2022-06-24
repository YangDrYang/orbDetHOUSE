#ifndef __GRAVITY_HPP__
#define __GRAVITY_HPP__

#include "eigenIncluder.hpp"
#include "streamTrace.hpp"
#include "satRefSys.hpp"
#include "jplEph.hpp"
#include "enums.h"
#include <string>
using std::string;

/** Struct to save Earth gravity coefficients
*/
struct EGMCoef
{
	MatrixXd	smn = MatrixXd::Zero(361,361);
	MatrixXd	cmn = MatrixXd::Zero(361,361);
};


MatrixXd Legendre(
	int 	n,		
	int 	m,		
	double 	phi);

MatrixXd LegendreD(
	int             m,		///< Maximum degree
	int             n,		///< Maximum order
    MatrixXd        pnm,    ///< Normalised Legendre polinomial matrix
	double          phi);   ///< Geocentric latitude in radian

struct EarthGravityDeg
{
	int		mMax = 15;
	int		nMax = 15;
};

struct EarthGravMdlOpt
{
	E_GravMdl			earthGravMdl;
	EarthGravityDeg		earthGravAccDeg;
	EarthGravityDeg		earthGravSTMDeg;		
};

struct EarthTidesOpt
{
	bool		flagSolidEarthTides	= 0;						///< Consider solid Earth tides or not   
	bool		flagOceanTides		= 0;						///< Consider ocean tides or not   
	E_TidesMdl	earthModel			= E_TidesMdl::ElasticMdl;	///< Elastic model or anelastic model
};

/** The gravity model of the central body
*/
struct GravityModel
{
public:

    // GravityModel();
    GravityModel(
        EarthGravMdlOpt         gravMdlOpt,
        EarthTidesOpt           tidesOpt,
        EGMCoef         		egmCoef,
		void*					pJPLEph);
    // virtual ~GravityModel();

public:
    void setModelType(E_GravMdl modelType);

    E_GravMdl getModelType();

    EarthGravityDeg mEarthGravAccDeg = {15, 15};
	EarthGravityDeg mEarthGravSTMDeg = { 2,  2};

public:
    /* Computes the acceleration due to the central body gravity field of the central body
    *
    */
    Vector3d centralBodyGravityAcc(
        Trace&                  trace,          ///< Trace to output to (similar to cout)
        const double            mjdUTC,
        const double*           erpv,
        const Vector3d&         rSat,           ///< Satellite position vector in the inertial system
        const Matrix3d&         mECI2BF,        ///< Transformation matrix from ECI to central body fixed system
		bool				 	bVarEq);		///< bVarEq = 1 if for the variational equation


    void solidEarthTidesCorrection(
        Trace&                  trace,
		double				    mjdUTC,			    ///< 
		EGMCoef&             	egmCoef,			///< Struct of Earth gravity coefficients
        IERS		 		    iersIns,			///< Instance of IERS class
		const Vector3d&		    vecRAESun,			///< Rho, azimuth and altitude information of Sun
		const Vector3d&		    vecRAEMoon);		///< Rho, azimuth and altitude information of Moon
    
    void oceanTidesCorrection(
        Trace&                  trace,
		double				    mjdUTC,			    ///< 
		EGMCoef&             	egmCoef,			///< Struct of Earth gravity coefficients
		const Vector3d&		    vecRAESun,			///< Rho, azimuth and altitude information of Sun
		const Vector3d&		    vecRAEMoon);		///< Rho, azimuth and altitude information of Moon
	

	Vector3d relativityEffectsAcc(
		Trace&    				trace,	 		///< Trace to output to (similar to cout)
		const Vector3d& 		rSat,			///< Inertial position of satellite (m)
		const Vector3d& 		vSat);			///< Inertial velocity of satellite (m/s)	    
	
	JPLEphemeris*		mJPLEphemeris;
	
protected:

    E_GravMdl     		mGravModelType;
    EarthTidesOpt       mEarthTidesOpt;
    double              mGM;
    double              mRadius;
    EGMCoef             mEGMCoef;
};

struct ThirdBodyGravityOpt
{
	bool flagMercuryGrav	= false;
	bool flagVenusGrav		= false;
	bool flagEarthGrav		= false;
	bool flagMarsGrav		= false;
	bool flagJupiterGrav	= false;
	bool flagSaturnGrav		= false;
	bool flagUranusGrav		= false;
	bool flagNeptuneGrav	= false;
	bool flagPlutoGrav		= false;
	bool flagMoonGrav		= true;
	bool flagSunGrav		= true;
};

extern	EGMCoef		egm;

#endif
