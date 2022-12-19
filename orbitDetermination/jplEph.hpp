/* JPL planetary and lunar ephemerides */

#ifndef JPLEPH_HPP
#define JPLEPH_HPP

// #include "eigenIncluder.hpp"
#include "Eigen/Dense"

/*
*                          
*   5 = jupiter          12 = solar-system barycenter             
*   6 = saturn           13 = earth-moon barycenter               
*   7 = uranus           14 = nutations (longitude and obliq)     
*                        15 = librations, if on eph. file         
*                        16 = lunar mantle omega_x,omega_y,omega_z
*                        17 = TT-TDB, if on eph. file             
*   If nutations are wanted, set ntarg = 14.                      
*   For librations, set ntarg = 15. set ncent= 0.                 
*   For TT-TDB,  set ntarg = 17.                                  
*   I've not actually seen an ntarg = 16 case yet.)               
*/
enum E_SolarSysStarType
{
    eMercury                =   1,   ///< Mercury
    eVenus                  =   2,   ///< Venus
    eEarth                  =   3,   ///< Earth
    eMars                   =   4,   ///< Mars
    eJupiter                =   5,   ///< Jupiter
    eSaturn                 =   6,   ///< Saturn
    eUranus                 =   7,   ///< Uranus
    eNeptune                =   8,   ///< Neptune
    ePluto                  =   9,   ///< Pluto
    eMoon                   =   10,  ///< Moon
    eSun                    =   11,  ///< Sun

    eSolarSysBarycenter     =   12,   ///< Solar System barycenter
    eEarthMoonBaryCenter    =   13,   ///< Earth Moon barycenter
    eJplNutation            =   14,   ///< Nutations
    eJPLLnrLibration        =   15,   ///< Lunar Librations
    eJPLLunarMantle         =   16,   ///< Lunar Mantle omega_x,omega_y,omega_z
    eJPLTT_TDB              =   17,   ///< TT-TDB, if on eph. file

    eSelfDefinedStarType    =   99,
};

/** Read JPL Ephemeris
*/
struct JPLEphemeris
{
    JPLEphemeris(void* pJPLEph);
    void*       mpJPLEph;          ///< a Pointer to The jpl_eph_data Structure

    double      getJPLEphemerisStartJD();
    double      getJPLEphemerisEndJD();
    double      getJPLEphemerisStep();

    /*
    *                 JPL_EPHEM_START_JD               0
    *                 JPL_EPHEM_END_JD                 8
    *                 JPL_EPHEM_STEP                  16
    *                 JPL_EPHEM_N_CONSTANTS           24
    *                 JPL_EPHEM_AU_IN_KM              28
    *                 JPL_EPHEM_EARTH_MOON_RATIO      36
    *                 JPL_EPHEM_IPT_ARRAY             44
    *                 JPL_EPHEM_EPHEMERIS_VERSION    224
    *                 JPL_EPHEM_KERNEL_SIZE          228
    *                 JPL_EPHEM_KERNEL_RECORD_SIZE   232
    *                 JPL_EPHEM_KERNEL_NCOEFF        236
    *                 JPL_EPHEM_KERNEL_SWAP_BYTES    240
    */
    double getJPLEphData(
        int                     type);

    /*
    * Computes the position of the target object with respect to reference object at the specified julian ephemeris date.
    *  The target and reference codes should be chosen from the JPL planetary ephemeris data request codes.
    */
    void getJPLEphemeris (
        double                      jd,         ///< Julian_TT
        E_SolarSysStarType          target,     ///< Star Needs to Calculate the Velocity and Position
        E_SolarSysStarType          center,     ///< Reference Body
        Vector3d&                   pos);       ///< Unit: m

    void getJPLEphemeris (
        double                      jd,
        E_SolarSysStarType          target,
        E_SolarSysStarType          center,
        Vector3d&                   pos,
        Vector3d&                   vel);       ///< vel (m/s)

private:
    double      mEphStartJD;       ///< The start dates are given in Julian Ephemeris Date (JD format of TDT)
    double      mEphEndJD;         ///< The stop dates are given in Julian Ephemeris Date (JD format of TDT)
    double      mEphStep;

};

extern  void*   pJPLEph;           ///< a pointer to JPL ephemerides

#endif //JPLEPH_HPP
