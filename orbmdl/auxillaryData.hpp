
#ifndef AUXILLARYDATA_HPP
#define AUXILLARYDATA_HPP

#include "constants.hpp"
#include <iostream>
#include <string>

using namespace std;

struct erpd_t
{
    /// earth rotation parameter data type
    double mjd;      ///< mjd (days)
    double xp, yp;   ///< pole offset (rad)
    double xpr, ypr; ///< pole offset rate (rad/day)
    double ut1_utc;  ///< ut1-utc (s)
    double lod;      ///< length of day (s/day)
};
struct erp_t
{
    /// earth rotation parameter type
    int n, nmax;  ///< number and max number of data
    erpd_t *data; ///< earth rotation parameter data
};

int readerp(string file, erp_t *erp);
int geterp_from_utc(
    const erp_t *erp,
    double leapSec,
    double mjd,
    double *erpv);
#endif