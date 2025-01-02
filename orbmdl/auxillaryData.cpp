#include "auxillaryData.hpp"

/* read earth rotation parameters ----------------------------------------------
 * read earth rotation parameters
 * args   : char   *file       I   IGS ERP file (IGS ERP ver.2)
 *          erp_t  *erp        O   earth rotation parameters
 * return : status (1:ok,0:file open error)
 *-----------------------------------------------------------------------------*/
int readerp(string file, erp_t *erp)
{
    FILE *fp;
    erpd_t *erp_data;
    double v[14] = {0};
    char buff[256];

    if (!(fp = fopen(file.c_str(), "r")))
    {
        return 0;
    }

    while (fgets(buff, sizeof(buff), fp))
    {
        if (sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7, v + 8, v + 9, v + 10, v + 11, v + 12, v + 13) < 5)
        {
            continue;
        }
        if (erp->n >= erp->nmax)
        {
            erp->nmax = erp->nmax <= 0 ? 128 : erp->nmax * 2;

            erp_data = (erpd_t *)realloc(erp->data, sizeof(erpd_t) * erp->nmax);
            if (!erp_data)
            {
                free(erp->data);
                erp->data = nullptr;
                erp->n = 0;
                erp->nmax = 0;
                fclose(fp);
                return 0;
            }
            erp->data = erp_data;
        }
        erp->data[erp->n].mjd = v[0];
        erp->data[erp->n].xp = v[1] * 1E-6 * AS2R;
        erp->data[erp->n].yp = v[2] * 1E-6 * AS2R;
        erp->data[erp->n].ut1_utc = v[3] * 1E-7;
        erp->data[erp->n].lod = v[4] * 1E-7;
        erp->data[erp->n].xpr = v[12] * 1E-6 * AS2R;
        erp->data[erp->n].ypr = v[13] * 1E-6 * AS2R;
        erp->n++;
    }
    fclose(fp);
    return 1;
}
/* get earth rotation parameter values -----------------------------------------
 * get earth rotation parameter values
 * args   : erp_t  *erp        I   earth rotation parameters
 *          double time        I   time (modified julian date)
 *          double *erpv       O   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
 * return : status (1:ok,0:error)
 *-----------------------------------------------------------------------------*/
int geterp_from_utc(
    const erp_t *erp,
    double leapSec,
    double mjd,
    double *erpv)
{
    double day;

    if (erp->n <= 0)
        return 0;
    // cout << "erp data[0] \t" << erp->data[0].mjd << "\t"
    //      << "erp data[n-1]" << erp->data[erp->n - 1].mjd << endl;
    if (mjd <= erp->data[0].mjd)
    {
        day = mjd - erp->data[0].mjd;

        erpv[0] = erp->data[0].xp + erp->data[0].xpr * day;
        erpv[1] = erp->data[0].yp + erp->data[0].ypr * day;
        erpv[2] = erp->data[0].ut1_utc - erp->data[0].lod * day;
        erpv[3] = erp->data[0].lod;
        erpv[4] = leapSec;

        return 1;
    }

    if (mjd >= erp->data[erp->n - 1].mjd)
    {
        day = mjd - erp->data[erp->n - 1].mjd;

        erpv[0] = erp->data[erp->n - 1].xp + erp->data[erp->n - 1].xpr * day;
        erpv[1] = erp->data[erp->n - 1].yp + erp->data[erp->n - 1].ypr * day;
        erpv[2] = erp->data[erp->n - 1].ut1_utc - erp->data[erp->n - 1].lod * day;
        erpv[3] = erp->data[erp->n - 1].lod;
        erpv[4] = leapSec;

        return 1;
    }

    int i;
    int j = 0;
    int k = erp->n - 1;
    while (j < k - 1)
    {
        i = (j + k) / 2;

        if (mjd < erp->data[i].mjd)
            k = i;
        else
            j = i;
    }

    double a;
    if (erp->data[j].mjd == erp->data[j + 1].mjd)
    {
        a = 0.5;
    }
    else
    {
        a = (mjd - erp->data[j].mjd) / (erp->data[j + 1].mjd - erp->data[j].mjd);
    }

    erpv[0] = (1 - a) * erp->data[j].xp + a * erp->data[j + 1].xp;
    erpv[1] = (1 - a) * erp->data[j].yp + a * erp->data[j + 1].yp;
    erpv[2] = (1 - a) * erp->data[j].ut1_utc + a * erp->data[j + 1].ut1_utc;
    erpv[3] = (1 - a) * erp->data[j].lod + a * erp->data[j + 1].lod;
    erpv[4] = leapSec;

    return 1;
}
