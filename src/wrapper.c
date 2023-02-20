#include <stddef.h>
#include "wbgt.h"
/*
 * All parameters below are the same as those defined in wbgt.c, with the
 * exception of num_obs, which counts the number of observations and is used
 * for properly vectorizing the R code
 */
void calc_wbgt_vector(size_t n, int *year, int *month, int *day, int *hour, int *minute, int *gmt, int *avg,
    double *lat, double *lon, double *solar, double *pres, double *Tair, double *relhum, double *speed, double *zspeed, 
    double *dT, int *urban, double *est_speed,
    double *Tg, double *Tnwb, double *Tpsy, double *Twbg,
    int *status)
{
  // size_t n = *num_obs;
  // double est_speed = 0.;
  for (size_t i = 0; i < n; ++i)
  {
    status[i] = calc_wbgt(year[i], month[i], day[i], hour[i], minute[i], gmt[i], avg[i],
        lat[i], lon[i], solar[i], pres[i], Tair[i], relhum[i], speed[i], zspeed[i],
        dT[i], urban[i],
        &est_speed[i], &Tg[i], &Tnwb[i], &Tpsy[i], &Twbg[i]);
  }
}