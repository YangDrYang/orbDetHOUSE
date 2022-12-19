
#ifndef ___COMMON_HPP__
#define ___COMMON_HPP__

#include <memory>

/* constants/macros ----------------------------------------------------------*/
#define SQR(x) ((x) * (x))
#define POW2(x) ((x) * (x))
#define POW4(x) ((x) * (x) * (x) * (x))
#define SQRT(x) ((x) <= 0.0 ? 0.0 : sqrt(x))
#define ROUND(x) (int)floor((x) + 0.5)
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define SWAP(x, y)   \
	do               \
	{                \
		double tmp_; \
		tmp_ = x;    \
		x = y;       \
		y = tmp_;    \
	} while (0)
#define SGN(x) ((x) <= 0.0 ? -1.0 : 1.0)

#include "satRefSys.hpp"

/* coordinates transformation ------------------------------------------------*/

void eci2ecef_sofa(
	const double mjdUTC,
	IERS &iersIns,
	Matrix3d &U,
	Matrix3d &dU);

void eci2ecefVec_sofa(
	const double mjdUTC,
	IERS &iersIns,
	Vector6d &rvSat_eci,
	Vector6d &rvSat_ecef);

void ecef2eciVec_sofa(
	const double mjdUTC,
	IERS &iersIns,
	Vector6d &rvSat_ecef,
	Vector6d &rvSat_eci);

void ecef2pos(const double *r, double *pos);
void ecef2pos(Vector3d &r, double *pos);
void pos2ecef(const double *pos, double *r);

double geodist(Vector3d &rs, Vector3d &rr, Vector3d &e);

// forward declarations
struct prcopt_t;
struct erp_t;
struct StationOptions;
struct SatSys;

int readerp(string file, erp_t *erp);
int geterp_from_utc(const erp_t *erp, double leapSec, double mjdUTC, double *val);

#endif
