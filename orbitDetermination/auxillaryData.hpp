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