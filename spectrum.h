#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define RESE 600
#define REST 1000

/* normalized LL energy spectrum - 15.4 average energy, 3.8: betha, 4802: normalization factor*/
#define LL_energy_spectrum(E) pow(E, 3.8)*exp(-(1.0+3.8)*E/15.4)/4802.516160
/* normalized LL arrival time spectrum - approximated by log-normal-dist */
#define MY -1.0324
#define SIGMA 0.9134
#define LL_time_spectrum(t) exp( - (log(t)-MY)*(log(t)-MY)/(2*SIGMA*SIGMA) ) / (t*SIGMA*sqrt(2*M_PI))
/*gaussian (around 0): error for the energy smearing: sqrt(E) -> sigma2=E*/
#define GAUSS(i, sigma2) 1/sqrt(2.0*M_PI*sigma2)*exp(-(i*i)/(2.0*sigma2))

#ifndef SPECTRUM_H_
#define SPECTRUM_H_

double time_shift(double t, double E, double mass, double dist);
double LLSpectrumTotal (double t, double E, double mass, double dist);
void ProbFirstHitDist (double mass, double dist, double events, double *result);
void correlation(double mass, double dist, double events, double *newSpec);
void generateDist(double mass, double dist, double events, double *distribution, double *triggerEfficiency, bool energyRes);

#endif 
