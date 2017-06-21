/*
 * The MIT License (MIT)
 * 
 * Copyright (c) 2016 Max Lieblich
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *   
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#include "wbgt.h"
#include <R.h>

/*
 * All parameters below are the same as those defined in wbgt.c, with the
 * exception of num_obs, which counts the number of observations and is used
 * for properly vectorizing the R code
 */
void wbgt(int *num_obs, int *year, int *month, int *day, int *hour, int *minute, int *gmt, int *avg, 
    float *lat, float *lon, float *solar, float *pres, float *Tair, float *relhum, float *speed, float *zspeed, 
    float *dT, int *urban, 
    float *est_speed, float *Tg, float *Tnwb, float *Tpsy, float *Twbg,
    int *status)
{
  int n = *num_obs;
  for (int i = 0; i < n; ++i)
  {
    status[i] = calc_wbgt(year[i], month[i], day[i], hour[i], minute[i], gmt[i], avg[i],
        lat[i], lon[i], solar[i], pres[i], Tair[i], relhum[i], speed[i], zspeed[i],
        dT[i], urban[i],
        est_speed + i, Tg + i, Tnwb + i, Tpsy + i, Twbg + i);
  }
}
