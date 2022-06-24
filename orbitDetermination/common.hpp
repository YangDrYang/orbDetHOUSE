
#ifndef ___COMMON_HPP__
#define ___COMMON_HPP__

#include <memory>

/* constants/macros ----------------------------------------------------------*/
#define SQR(x)      ((x)*(x))
#define POW2(x)     ((x)*(x))
#define POW4(x)     ((x)*(x)*(x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define ROUND(x)    (int)floor((x)+0.5)
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)
#define SGN(x)      ((x)<=0.0?-1.0:1.0)

#include "eigenIncluder.hpp"
#include "satRefSys.hpp"
#include "gTime.hpp"
#include "enums.h"



/* coordinates transformation ------------------------------------------------*/

void ecef2enu(const double *pos, const double *r, double *e);
void enu2ecef(const double *pos, const double *e, double *r);
void xyz2enu (const double *pos, double *E);


void eci2ecef(
	const GTime		tutc,	
	const double*	erpv,	
	Matrix3d&		U,		
	Matrix3d*		dU		= nullptr,		
	double*			gmst	= nullptr);

void eci2ecef_sofa(
	const double	mjdUTC,
	IERS&			iersIns,	
	Matrix3d&		U,		
	Matrix3d&		dU);

void eci2ecefVec_sofa(
	const double	mjdUTC,
	IERS&			iersIns,	
	Vector6d&		rvSat_eci,
	Vector6d&		rvSat_ecef);					

void ecef2eciVec_sofa(
	const double	mjdUTC,
	IERS&			iersIns,	
	Vector6d&		rvSat_ecef,
	Vector6d&		rvSat_eci);

void ecef2pos(const double *r, double *pos);
void ecef2pos(Vector3d& r, double *pos);
void pos2ecef(const double *pos, double *r);

double geodist(Vector3d& rs, Vector3d& rr, Vector3d& e);

double sagnac(
	Vector3d& rSource,
	Vector3d& rDest);

//forward declarations
struct prcopt_t;
struct erp_t;
struct StationOptions;
struct Obs;
struct SatSys;

/* satellites, systems, codes functions --------------------------------------*/

double satazel(const double *pos, const double *e, double *azel);
double geodist(const double *rs, const double *rr, double *e);

unsigned int getbitu	(const unsigned char *buff, int  pos, int len);
int          getbits	(const unsigned char *buff, int  pos, int len);
unsigned int getbituInc	(const unsigned char *buff, int& pos, int len);
int          getbitsInc	(const unsigned char *buff, int& pos, int len);
int setbitsInc(unsigned char *buff,int pos,int len,const int var);
int setbituInc(unsigned char *buff,int pos,int len,const unsigned int var);

unsigned int crc24q (const unsigned char *buff, int len);

/* positioning models --------------------------------------------------------*/
void dops(int ns, const double *azel, double elmin, double *dop);

int  readblq(string file, const char *sta, double *otlDisplacement);
int  readerp(string file, erp_t *erp);
int  geterp (const erp_t *erp, GTime time, double *val);
int  geterp_from_utc (const erp_t *erp, double leapSec, double mjdUTC, double *val);


int satexclude(SatSys& sat, E_Svh svh);

extern int		epoch;
extern GTime	tsync;

void replaceTimes(
	string&						str,		///< String to replace macros within
	boost::posix_time::ptime	time_time);	///< Time to use for replacements

#endif
