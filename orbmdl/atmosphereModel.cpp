#include "atmosphereModel.hpp"

#define R_EARTH 6378.137e3 // from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

static void MJDtoDOY(double mjd, int &year, int &doy, double &sec);
static double calculateDensity(double t, const Eigen::VectorXd &X, IERS iersInstance);
static long ymd_to_mjd(int year, int month, int day);
static void mjd_to_ymd(long mjd, int *year, int *month, int *day);
static vector<double> cartesianGeodetic(VectorXd Xeci, double t, IERS iersInstance);

void eci2ecefVec_sofa(
    const double mjdUTC,
    IERS &iersIns,
    VectorXd &rvSat_eci,
    VectorXd &rvSat_ecef);

Vector3d calculateDragForce(const Vector3d &rSat,
                            const Vector3d &vSat,
                            IERS iers,
                            double Area,
                            double Cd,
                            double t)
{

    VectorXd X(6);
    X.head(3) = rSat;
    X.tail(3) = vSat;

    // calculate density from model
    double rho = calculateDensity(t, X, iers);
    Vector3d V = X.tail(3);

    // calculate drag force from density
    Vector3d drag = -0.5 * rho * Area * Cd * pow(V.norm(), 1) * V;
    // ToDo: currently the orbit is decaying excessively, (below 0 altitude),
    // check units/figure out why this is occuring
    // drag << 0, 0, 0;
    return drag;
}

// from http://www.leapsecond.com/tools/gdate.c
void MJDtoDOY(double mjd, int &year, int &doy, double &sec)
{
    int month, day;

    long lnMJD = (long)mjd;
    mjd_to_ymd(lnMJD, &year, &month, &day);
    doy = lnMJD - ymd_to_mjd(year, 1, 1) + 1;
    sec = fmod(mjd, 1) * 86400;
    return;
}

double calculateDensity(double t, const VectorXd &X, IERS iersInstance)
{
    // get density from atmospheric model
    // https://github.com/magnific0/nrlmsise-00/blob/master/nrlmsise-00.h

    // input & output structs
    struct nrlmsise_input input;
    struct nrlmsise_output output;
    struct nrlmsise_flags flags;
    // struct ap_array aph;

    // set aph array
    // for (int i = 0; i < 7; i++) aph.a[i] = 100;

    // set flags
    for (int i = 0; i < 24; i++)
        flags.switches[i] = 1;

    // convert MJD
    int year, doy;
    double sec;
    MJDtoDOY(t, year, doy, sec);

    // <lat, long, alt>
    vector<double> position = cartesianGeodetic(X, t, iersInstance);
    double g_lat = position[0];
    double g_long = position[1];
    double alt = position[2];

    // cout << "lat = " << g_lat << ", long = " << g_long <<  ", alt = " << alt << endl;

    double lst = sec / 3600 + g_long / 15;

    input = {
        .year = year,
        .doy = doy,
        .sec = sec,
        .alt = alt / 1000,
        .g_lat = g_lat,
        .g_long = g_long,
        .lst = lst,
        .f107A = 150.0,
        .f107 = 150.0,
        .ap = 4.0};
    gtd7(&input, &flags, &output);
    return output.d[5];
}

vector<double> cartesianGeodetic(VectorXd Xeci, double t, IERS iersInstance)
{
    VectorXd Xecef(6);

    // convert from eci to ecef
    eci2ecefVec_sofa(t, iersInstance, Xeci, Xecef);

    // calculate geodetic coords from cartesian
    double latitude = atan2(Xecef(2), Xecef.head(2).norm());
    double longitude = atan(Xecef(1) / Xecef(0));
    double altitude = Xecef.norm() - R_EARTH;

    // convert from rad to deg
    latitude *= 180 / M_PI;
    longitude *= 180 / M_PI;

    // check that value are reasonable
    assert(latitude <= 90 && latitude >= -90);
    assert(longitude <= 180 && longitude >= -180);
    // assert(altitude >= 0);

    // structure return value
    vector<double> returnVal{latitude, longitude, altitude};
    return returnVal;
}

// helper functions
long ymd_to_mjd(int year, int month, int day)
{
    long Y = year, M = month, D = day;
    long mjd =
        367 * Y - 7 * (Y + (M + 9) / 12) / 4 - 3 * ((Y + (M - 9) / 7) / 100 + 1) / 4 + 275 * M / 9 + D + 1721029 - 2400001;
    return mjd;
}
void mjd_to_ymd(long mjd, int *year, int *month, int *day)
{
    long J, C, Y, M;

    J = mjd + 2400001 + 68569;
    C = 4 * J / 146097;
    J = J - (146097 * C + 3) / 4;
    Y = 4000 * (J + 1) / 1461001;
    J = J - 1461 * Y / 4 + 31;
    M = 80 * J / 2447;
    *day = (int)(J - 2447 * M / 80);
    J = M / 11;
    *month = (int)(M + 2 - (12 * J));
    *year = (int)(100 * (C - 49) + Y + J);
    return;
}