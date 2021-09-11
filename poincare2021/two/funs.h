#ifndef FUNS_H
#define FUNS_H

#include "poincare.h"

void find_peak(const double t, double peak_v[], double peak_now_v[], double peak_past_v[], double peak_old_v[], double peak_d[], double peak_now_d[], double peak_past_d[], double peak_old_d[], const double y[], int ind_v[], int ind_d[]);
int kuramoto(double t, const double y[], double f[], void *params);
void write_curve(const double y[], const double t, int light);
void write_peaks(const double peak_v[], const double peak_d[], double period_v[], double period_d[], double phase_gap[], const int peak_num);
void write_peaks_p(const double k_vd, const double peak_v[], const double peak_d[], double period_v[], double period_d[], double phase_gap[], const int peak_num);

#endif