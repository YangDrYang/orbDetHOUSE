#ifndef CONFIG_HPP
#define CONFIG_HPP

#define BETTER_ENUMS_DEFAULT_CONSTRUCTOR(Enum) \
public:                                        \
    Enum() = default;

#include "enum.h" //BETTER_ENUM
#include <string>
using std::string;

BETTER_ENUM(E_SRPModels, int,
            CANNONBALL,
            BOXWING,
            ECOM,
            ECOM2)

BETTER_ENUM(E_GravMdl, short int,
            EGM08Model = 0,
            GGM03SModel = 1,
            GGM05SModel = 2)

BETTER_ENUM(E_TidesMdl, short int,
            ElasticMdl = 0,
            AnelasticMdl = 1)

/** Options associated with orbital force models
 */
struct ForceModels
{
    bool earth_gravity = true;
    E_GravMdl gravity_model = E_GravMdl::GGM03SModel;
    bool solid_earth_tide = false;
    bool ocean_tide_loading = false;
    bool third_body_attraction = false;
    bool third_body_sun = false;
    bool third_body_moon = false;
    bool third_body_planet = false;
    bool relativity_effect = false;
    bool solar_radiation_pressure = false;
    bool thermal_emission = false;
    bool earth_albedo = false;
    bool infrared_radiation = false;
    bool antenna_thrust = false;
    bool empirical_acceleration = false;
    bool satellite_manoeuvre = false;
    double satMass = 100;
    E_SRPModels srpMdlName = E_SRPModels::CANNONBALL;
    double srpArea = 5;
    double srpCoef = 1;
    int egmAccDeg = 12;
    int egmAccOrd = 12;
    int egmSTMDeg = 4;
    int egmSTMOrd = 4;
    string odeInteg = "rkf78";
};
#endif