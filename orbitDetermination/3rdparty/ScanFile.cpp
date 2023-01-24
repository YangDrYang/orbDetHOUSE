#ifndef _INCLUDE_POSIX_SOURCE
#define _INCLUDE_POSIX_SOURCE
#endif

/*
 *   include files
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include "atest.hpp"
#include "gTime.hpp"
#include "constants.hpp"

/** YMDHMS to Julian Day
 * return   :       julian day
 */
double ymdhms2jd(
    const double time[6]) ///< civil date time [YMDHMS]
{
    double i;
    double j;

    double yr = time[0];
    double mon = time[1];
    //  double day  = time[2];
    double hr = time[3];
    double min = time[4];
    double sec = time[5];

    if (yr <= 0 || yr >= 2099)
        return 0;

    if (mon > 2)
    {
        i = yr;
        j = mon;
    }
    else
    {
        i = yr - 1;
        j = mon + 12;
    }

    double day = time[2] + hr / 24 + min / 24 / 60 + sec / 24 / 60 / 60;

    double jd = (floor(365.25 * i) + floor(30.6001 * (j + 1)) + day + 1720981.5);

    //  BOOST_LOG_TRIVIAL(debug) << "YMDH to JD: year=" << yr << " mon=" << mon << " day=" << day << " hour=" << hr << ", jd=" << jd;
    return jd;
}

int ScanFile(
    EpochDatum **data, int max_data, char *file_name, int prn)
{

    char buf[LINE_MAX + 1];
    char tmp[LINE_MAX + 1];

    double epoch;
    double second;
    double sig;

    int day;
    int hour;
    int minute;
    int month;
    int year;

    int count_data = 0;
    int i;
    int id;
    int id_count = 0;
    int ios = 0;
    int store_id_count = 0;
    int sig_count = 0;

    FILE *filptr = NULL;

    /*
     *  1.0  Open file
     */

    if ((filptr = fopen(file_name, "r")) == NULL)
    {
        ios = -1;
        return (ios);
    }

    /*
     *  2.0  Read all lines
     */

    while (fgets(buf, LINE_MAX, filptr) != NULL)
    {

        /*
         *  3.0  ID location in PRN list
         */

        if (buf[0] == '+' && buf[1] == ' ' && strlen(buf) >= 60)
        {
            for (i = 9; i < strlen(buf) - 3; i = i + 3)
            {
                strncpy(tmp, buf + i, 3);
                tmp[3] = '\0';
                if (tmp[0] == ' ' || tmp[0] == 'G' || tmp[0] == 'g')
                {
                    id = atoi(tmp + 1);
                    if (id == prn)
                    {
                        store_id_count = id_count;
                    }
                }
                ++id_count;
            }

            /*
             *  4.0  Trap the sigma
             */
        }
        else if (buf[0] == '+' && buf[1] == '+' && strlen(buf) >= 60)
        {
            for (i = 9; i < strlen(buf) - 3; i = i + 3)
            {
                if (sig_count == store_id_count)
                {
                    strncpy(tmp, buf + i, 3);
                    tmp[3] = '\0';
                    id = atoi(tmp);
                    if (id > 0)
                    {
                        sig = pow(2, id) / 1000.0;
                    }
                    else
                    {
                        sig = 1.0;
                    }
                }
                ++sig_count;
            }

            /*
             *  5.0  Trap epoch lines
             */
        }
        else if (buf[0] == '*')
        {
            sscanf(buf + 1, "%d %d %d %d %d %lf",
                   &year, &month, &day, &hour, &minute, &second);
            double ep[6] = {1.0 * year, 1.0 * month, 1.0 * day, 1.0 * hour, 1.0 * minute, second};
            epoch = ymdhms2jd(ep) - JD2MJD;
            // epoch = ymdmjd( year, month, day ) + hmsday( hour, minute, second );

            /*
             *  6.0  Trap data lines
             */
        }
        else if (
            buf[0] == 'P' &&
            atof(buf + 2) == prn &&
            count_data < max_data &&
            (data[count_data] = (EpochDatum *)malloc(sizeof(EpochDatum))) != NULL)
        {
            data[count_data]->epoch = epoch;
            sscanf(
                buf + 4, "%lf %lf %lf",
                &(data[count_data]->xyz[0]),
                &(data[count_data]->xyz[1]),
                &(data[count_data]->xyz[2]));
            data[count_data]->sig[0] = sig;
            data[count_data]->sig[1] = sig;
            data[count_data]->sig[2] = sig;

            ++count_data;
        }
    }

    fclose(filptr);

    return (count_data);
}