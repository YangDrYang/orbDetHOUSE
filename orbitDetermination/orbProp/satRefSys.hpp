/*
* Modified version of satRefSys.cpp by O. Montenbruck, E. Gill
* Transformations betweeen celestial and terrestrial reference systems
*/

#ifndef INC_SATREFSYS_HPP
#define INC_SATREFSYS_HPP

#include "eigenIncluder.hpp"
#include "constants.hpp"


/** Management of IERS time and polar motion data  
*/

struct IERS 
{

	static void Set (double UT1_UTC,          // Set Earth rotation parameters
					double UTC_TAI,          // (UT1-UTC [s],UTC-TAI [s],
					double x_pole,           //  x [radian], y [radian])
					double y_pole,
					double lod);             //length of day [second/day]

	static const double  TT_TAI;              //  TT-TAI time difference 32.184s
	static const double GPS_TAI;              // GPS-TAI time difference -19s
	static double UTC_TAI(double Mjd_UTC);    // UTC_TAI time difference [s]
	static double UT1_TAI(double Mjd_UTC);    // UT1-UTC time difference [s]

	static const double  TT_GPS;              //  TT-GPS time difference 51.184s
	static const double TAI_GPS;              // TAI-GPS time difference 19s
	static double UTC_GPS(double Mjd_UTC);    // UTC_GPS time difference [s]
	static double UT1_GPS(double Mjd_UTC);    // UT1-GPS time difference [s]

	static double  TT_UTC(double Mjd_UTC);    //  TT-UTC time difference [s]
	static double TAI_UTC(double Mjd_UTC);    // TAI-UTC time difference [s]
	static double GPS_UTC(double Mjd_UTC);    // GPS-UTC time difference [s]
	static double UT1_UTC(double Mjd_UTC);    // UT1-UTC time difference [s]

	static double x_pole(double Mjd_UTC);     // Pole coordinate [rad]
	static double y_pole(double Mjd_UTC);     // Pole coordinate [rad]
	static double lod(double Mjd_UTC);        // Length of day [s/d]

private:
	static double UT1_TAI_;                   // UT1-TAI time difference [s]
	static double UTC_TAI_;                   // UTC-TAI time difference [s]
	static double x_pole_;                    // Pole coordinate [rad]
	static double y_pole_;                    // Pole coordinate [rad]
	static double lod_;                        // length of day [s/d]
};


Matrix3d R_x(
	double    Angle);

Matrix3d R_y(
	double    Angle);

Matrix3d R_z(
	double    Angle);


#endif 
