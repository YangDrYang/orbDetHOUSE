#include <iostream>

#include "constants.hpp"
#include "jplEph.hpp"
#include "jpl_eph.hpp"



/*
* Class type: The class of Read JPL Ephemeris
*/
JPLEphemeris::JPLEphemeris(void* pJPLEph)
{
	char nams[2000][6];
	double vals[2000];
	mpJPLEph = pJPLEph;
	
	int result = jpl_init_error_code();
	switch (result)
	{
		case 0:
			mEphStartJD	= jpl_get_double(mpJPLEph, JPL_EPHEM_START_JD);
			mEphEndJD	= jpl_get_double(mpJPLEph, JPL_EPHEM_END_JD);
			mEphStep	= jpl_get_double(mpJPLEph, JPL_EPHEM_STEP);
			break;
		case -1:			std::cout << "jplEphemeris Init: JPL_INIT_FILE_NOT_FOUND"		<< std::endl;		
		case -2:			std::cout << "jplEphemeris Init: JPL_INIT_FSEEK_FAILED"			<< std::endl;				//todo aaron, no breaks
		case -3:			std::cout << "jplEphemeris Init: JPL_INIT_FREAD_FAILED"			<< std::endl;		
		case -4:			std::cout << "jplEphemeris Init: JPL_INIT_FREAD2_FAILED"		<< std::endl;		
		case -5:			std::cout << "jplEphemeris Init: JPL_INIT_FILE_CORRUPT"			<< std::endl;		
		case -6:			std::cout << "jplEphemeris Init: JPL_INIT_MEMORY_FAILURE"		<< std::endl;		
		case -7:			std::cout << "jplEphemeris Init: JPL_INIT_FREAD3_FAILED"		<< std::endl;		
		case -8:			std::cout << "jplEphemeris Init: JPL_INIT_FREAD4_FAILED"		<< std::endl;		
		case -9:			std::cout << "jplEphemeris Init: JPL_INIT_NOT_CALLED"			<< std::endl;	
		case -10:			std::cout << "jplEphemeris Init: JPL_INIT_FREAD5_FAILED"		<< std::endl;	
		default:			std::cout << "jplEphemeris: Result is out of Known Situation"	<< std::endl;
	}
}

double JPLEphemeris::getJPLEphemerisStartJD()
{
	return mEphStartJD;
}

double JPLEphemeris::getJPLEphemerisEndJD()
{
	return mEphEndJD;
}

double JPLEphemeris::getJPLEphemerisStep()
{
	return mEphStep;
}

double JPLEphemeris::getJPLEphData(
	int type)
{
	return jpl_get_double(mpJPLEph, type);
}

void JPLEphemeris::getJPLEphemeris(
	double jd,
	E_SolarSysStarType target, 
	E_SolarSysStarType center, 
	Vector3d &pos)
{
	double r_p[6];
	int result = jpl_pleph(mpJPLEph, jd, target, center, r_p, 0);
	
	switch (result)
	{
		case 0:
			pos(0) = r_p[0] * AU;
			pos(1) = r_p[1] * AU;
			pos(2) = r_p[2] * AU;
			break;
		case -1:			std::cout << "getJPLEphemeris:JPL_EPH_OUTSIDE_RANGE"				<< std::endl;	
		case -2:			std::cout << "getJPLEphemeris:JPL_EPH_READ_ERROR"					<< std::endl;		
		case -3:			std::cout << "getJPLEphemeris:JPL_EPH_QUANTITY_NOT_IN_EPHEMERIS"	<< std::endl;	
		case -5:			std::cout << "getJPLEphemeris:JPL_EPH_INVALID_INDEX"				<< std::endl;		
		case -6:			std::cout << "getJPLEphemeris:JPL_EPH_FSEEK_ERROR"					<< std::endl;		
		default:			std::cout << "jplEphemeris: Result is out of Known Situation"		<< std::endl;
	}
}

void JPLEphemeris::getJPLEphemeris(
	double jd, 
	E_SolarSysStarType target, 
	E_SolarSysStarType center, 
	Vector3d &pos, 
	Vector3d &vel)
{
	double r_p[6];
	int result = jpl_pleph(mpJPLEph, jd, target, center, r_p, 1);
	
	switch (result)
	{
		case 0:
			pos(0) = r_p[0] * AU;
			pos(1) = r_p[1] * AU;
			pos(2) = r_p[2] * AU;

			vel(0) = r_p[3] * AUPerDay;
			vel(1) = r_p[4] * AUPerDay;
			vel(2) = r_p[5] * AUPerDay;
			break;
		case -1:		std::cout << "getJPLEphemeris: JPL_EPH_OUTSIDE_RANGE"				<< std::endl;	
		case -2:		std::cout << "getJPLEphemeris:JPL_EPH_READ_ERROR"					<< std::endl;	
		case -3:		std::cout << "getJPLEphemeris:JPL_EPH_QUANTITY_NOT_IN_EPHEMERIS"	<< std::endl;	
		case -5:		std::cout << "getJPLEphemeris:JPL_EPH_INVALID_INDEX"				<< std::endl;
		case -6:		std::cout << "getJPLEphemeris:JPL_EPH_FSEEK_ERROR"					<< std::endl;	
		default:		std::cout << "jplEphemeris: Result is out of Known Situation"		<< std::endl;
	}
}

