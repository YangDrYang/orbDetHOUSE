
#ifndef __CONSTANTS___HPP_
#define __CONSTANTS___HPP_

#include <map>

using std::map;

#define CLIGHT 299792458.0 /* speed of light (m/s) */

#define PI 3.141592653589793238462643383279502884197169399375105820974 /* pi */
#define D2R (PI / 180.0)                                               /* deg to rad */
#define R2D (180.0 / PI)                                               /* rad to deg */
#define SC2RAD 3.1415926535898                                         /* semi-circle to radian (IS-GPS) */
#define AU 149597870691.0                                              /* 1 AU (m) */
#define AS2R (D2R / 3600.0)                                            /* arc sec to radian */

#define OMGE 7.2921151467E-5 /* earth angular velocity (IS-GPS) (rad/s) */

#define EARTH_MASS 5.9722E24
#define EARTH_G_CONST 3.9860050E14 /* gravitational constant of earth (G * M_e) */
#define G_CONST 6.67408E-11
#define RE_GLO 6378136.0      /* radius of earth (m)            ref [2] */
#define MU_GPS 3.9860050E14   /* gravitational constant         ref [1] */
#define MU_GLO 3.9860044E14   /* gravitational constant         ref [2] */
#define MU_GAL 3.986004418E14 /* earth gravitational constant   ref [7] */
#define MU_CMP 3.986004418E14 /* earth gravitational constant   ref [9] */

#define MU MU_GPS
#define MOMENTUM_SCALE 10000000.0

#define RE_WGS84 6378137.0             /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84 (1.0 / 298.257223563) /* earth flattening (WGS84) */

#define JD2MJD 2400000.5 /* JD to MJD */
#define DayToSec 8.64e04
#define AUPerDay (AU / DayToSec) /* AU/Day (IAU 2009)[m/s] */

#define GM_Earth 3.986004418e14 ///< Geocentric gravitation constant (WGS84) [m^3/s^2]
/*
// moon parameter
*/
#define GM_Moon 4.9027949e12 ///< moon gravitation constant [m^3/s^2]
#define MoonRadius 1738200.0 ///< Equatorial radius of the Moon [m]
#define MoonMinRadius 1738200.0
#define Moon_J2 0.2027e-3

/*
// Jupiter parameter
*/
#define GM_Jupiter 1.26712000000000e+017
#define JupiterRadius 7.14920000000000e+007
#define JupiterMinRadius 6.68540000000000e+007
#define Jupiter_J2 0.01475

/*
// Mars parameter
*/
#define GM_Mars 4.28282868534000e+013
#define MarsRadius 3.39700000000000e+006
#define MarsMinRadius 3.37500000000000e+006
#define Mars_J2 1.96045e-3

/*
// Mercury parameter
*/
#define GM_Mercury 2.20320800000000e+013
#define MercuryRadius 2.43970000000000e+006
#define MercuryMinRadius 2.43970000000000e+006
#define Mercury_J2 50.3e-6

/*
// Neptune parameter
*/
#define GM_Neptune 6.87130000000000e+015
#define NeptuneRadius 2.52690000000000e+007
#define NeptuneMinRadius 2.48000000000000e+007
#define Neptune_J2 3.411e-3
/*
// Pluto parameter
*/
#define GM_Pluto 1.00907600000000e+012
#define PlutoRadius 1.16200000000000e+006
#define PlutoMinRadius 1.16200000000000e+006
#define Pluto_J2 0.0
/*
// Saturn parameter
*/
#define GM_Saturn 3.79340000000000e+016
#define SaturnRadius 6.02680000000000e+007
#define SaturnMinRadius 5.43640000000000e+007
#define Saturn_J2 0.016298

/*
// Sun parameter
*/
#define GM_Sun 1.327122E20    ///< Heliocentric gravitation constant [m^3/s^2]
#define SunRadius 695990000.0 ///< Equatorial radius of the Sun [m], Seidelmann 1992
#define SunMinRadius 695990000.0
#define SolarRadPreAtAU 4.560E-6 ///< Solar Radiation pressure at 1 AU [N/m^2] (~1367 W/m^2); IERS 96
#define Sun_J2 0.0
/*
// Uranus parameter
*/
#define GM_Uranus 5.80320000000000e+015
#define UranusRadius 2.55590000000000e+007
#define UranusMinRadius 2.49730000000000e+007
#define Uranus_J2 3.34343e-3

/*
 * Venus parameter
 */
#define GM_Venus 3.24858800000000e+014
#define VenusRadius 6.05190000000000e+006
#define VenusMinRadius 6.05190000000000e+006
#define Venus_J2 4.458e-6

#define mjdJ2000 51544.5         /* modified Julian date of J2000.0 */
#define Arcs 3600.0 * 180.0 / PI /* Arcseconds per radian */
#endif
