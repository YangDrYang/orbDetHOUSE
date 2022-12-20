
#include "auxillaryData.hpp"
#include "forceModels.hpp"
#include "jplEph.hpp"
#include "config.hpp"

#include <iostream>
#include <fstream>
using namespace std;

// /** Options associated with orbital force models
//  */
// struct ForceModels
// {
// 	bool earth_gravity = true;
// 	E_GravMdl gravity_model = E_GravMdl::GGM03SModel;
// 	bool solid_earth_tide = false;
// 	bool ocean_tide_loading = false;
// 	bool third_body_attraction = false;
// 	bool third_body_sun = false;
// 	bool third_body_moon = false;
// 	bool third_body_planet = false;
// 	bool relativity_effect = false;
// 	bool solar_radiation_pressure = false;
// 	bool thermal_emission = false;
// 	bool earth_albedo = false;
// 	bool infrared_radiation = false;
// 	bool antenna_thrust = false;
// 	bool empirical_acceleration = false;
// 	bool satellite_manoeuvre = false;
// 	double satMass = 100;
// 	E_SRPModels srpMdlName = E_SRPModels::CANNONBALL;
// 	double srpArea = 5;
// 	double srpCoef = 1;
// 	int egmAccDeg = 12;
// 	int egmAccOrd = 12;
// 	int egmSTMDeg = 4;
// 	int egmSTMOrd = 4;
// 	string odeInteg = "rkf78";
// };

int main()
{
	erp_t *erpt;
	/* initial condition for differential equation*/

	// // Initial value in ECEF, to be transformed to ECI
	// Vector6d rvECEF;
	// rvECEF = {-12053.996853e3,  21412.163227e3,  -9849.315798e3, -1032.8190896,   707.2236495,  2889.2728872};

	double leapSec = 32;
	double mjdUTC = 53300; // random number
	// double erpv[0] = {};
	geterp_from_utc(&erpt, leapSec, mjdUTC, erpv);

	// double dUT1_UTC = erpv[2];
	// double dUTC_TAI = -(19 + leapSec);
	// double xp = erpv[0];
	// double yp = erpv[1];
	// double lod = erpv[3];
	// IERS iersInstance;
	// iersInstance.Set(dUT1_UTC, dUTC_TAI, xp, yp, lod);
	// double mjdTT = mjdUTC + iersInstance.TT_UTC(mjdUTC) / 86400;

	// Matrix3d mECI2ECEF = Matrix3d::Identity();
	// Matrix3d mdECI2ECEF = Matrix3d::Identity();
	// Vector6d rvECI;
	// ecef2eciVec_sofa(mjdUTC, iersInstance, rvECEF, rvECI);

	VectorXd rvECI(5);
	rvECI(0) = 14226165.4977822;
	rvECI(1) = 20021922.7850642;
	rvECI(2) = -9875597.15080248;
	rvECI(3) = -1254.27669098652;
	rvECI(4) = 2274.30031195604;
	rvECI(5) = 2891.66233001166;

	/* initial condition for variational equation*/

	int nState = 6;
	int nVar = 7; // 6 state variables and one orbital parameter Cr
	VectorXd rvPhiS = VectorXd::Zero(nState + nState * nVar);
	// Create combined vector from epoch state, epoch transition matrix (=1) and epoch sensitivity matrix (=0)
	for (int i = 0; i < nState; i++)
	{
		rvPhiS(i) = rvECI(i);
		for (int j = 0; j < nVar; j++)
		{
			rvPhiS(nState * (j + 1) + i) = (i == j ? 1 : 0);
		}
	}

	ForceModels forceModelsOpt = {};

	/* setting options for the propagator */
	orbitProp.setPropOption(forceModelsOpt);

	/* initialising the propagator */
	orbitProp.initPropagator(rvECI, rvPhiS, mjdUTC, leapSec, erpt, egm, pJPLEph);
	orbitProp.updPropagator(mjdUTC);

	ofstream fileOPResults;
	fileOPResults.open("OPResults.txt");
	if (fileOPResults)
	{
		fileOPResults << "MJD(UTC)			x(ECI)			y(ECI)			z(ECI)			vx(ECI)			vy(ECI)			vz(ECI)\n";
		fileOPResults.close();
	}
	else
		cout << "Error openinng results file!\n";

	return 0;
}