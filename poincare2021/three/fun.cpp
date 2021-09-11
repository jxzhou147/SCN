#include <fstream>
#include <iostream>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "poincare.h"
#include "funs.h"

using namespace std;

int kuramoto(double t, const double y[], double f[], void *params)
{
    double TL = *(double *)params;
    //double k_vd = k_vd0 * (1 + y[0] / 2.0);
    double k_vd = k_vd0;
    /*
    if (fmod(t, TL) < (18.0))
        k_vd = 0.01;
    else
        k_vd = 0.05;
        */
    int light = 0;

    if ((t < 1000) && (fmod(t, TL) < (TL / 2.0)))
        light = 1;
    else
        light = 0;

    //f[0] = omega_v + K_f * light * phase(theta[0]) + k1 * sin(y[1] - y[0]);
    //f[0] = omega_v + K_f * light * phase_light(theta[0], daylength) + k1 * sin(y[1] - y[0]);
    //f[0] = omega_v + k1 * sin(y[1] - y[0]) + light * sin(Omega * t - y[0]);
    //f[0] = omega_v + k1 * sin(y[1] - y[0]) + light * 0.05;
    //f[0] = omega_v + light * 0.05;
    //f[0] = gam * y[0] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) - y[2] * omega_v + K * y[1] + K_f * light * phase_light(theta[0], daylength);
    //f[1] = gam * y[1] * (a - sqrt(pow(y[1], 2) + pow(y[3], 2))) - y[3] * omega_d + K * y[0];

    f[0] = gam * y[0] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) - y[2] * omega_v + k_dv * y[1] + K_f * light;
    f[1] = gam * y[1] * (a - sqrt(pow(y[1], 2) + pow(y[3], 2))) - y[3] * omega_d + k_vd * y[0];

    f[2] = gam * y[2] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) + y[0] * omega_v;
    f[3] = gam * y[3] * (a - sqrt(pow(y[1], 2) + pow(y[3], 2))) + y[1] * omega_d;

    return GSL_SUCCESS;
}

void find_peak(const double t, double peak_v[], double peak_now_v[], double peak_past_v[], double peak_old_v[], double peak_d[], double peak_now_d[], double peak_past_d[], double peak_old_d[], const double y[], int ind_v[], int ind_d[])
{
    peak_old_v[0] = peak_past_v[0]; // [0] is x, [1] is t
    peak_old_v[1] = peak_past_v[1];
    peak_past_v[0] = peak_now_v[0];
    peak_past_v[1] = peak_now_v[1];
    peak_now_v[0] = y[0];
    peak_now_v[1] = t;
    peak_old_d[0] = peak_past_d[0]; // [0] is x, [1] is t
    peak_old_d[1] = peak_past_d[1];
    peak_past_d[0] = peak_now_d[0];
    peak_past_d[1] = peak_now_d[1];
    peak_now_d[0] = y[1];
    peak_now_d[1] = t;

    if ((peak_past_v[0] > peak_old_v[0]) & (peak_past_v[0] > peak_now_v[0]))
    {
        peak_v[ind_v[0]] = peak_past_v[1];
        ind_v[0]++;
    }
    if ((peak_past_d[0] > peak_old_d[0]) & (peak_past_d[0] > peak_now_d[0]))
    {
        peak_d[ind_d[0]] = peak_past_d[1];
        ind_d[0]++;
    }
}

void write_curve(const double y[], const double t, int light)
{
    ofstream ofile;
    ofile.open("tcycle.csv", ios::app);
    ofile << t << '\t' << y[0] << '\t' << y[1] << '\t' << light << endl;
    ofile.close();
}

void write_peaks(const double peak_v[], const double peak_d[], double period_v[], double period_d[], double phase_gap[], const int peak_num)
{
    ofstream ofile;
    ofile.open("peaks.csv", ios::app);

    for (size_t i = 0; i < (peak_num - 1); i++)
    {
        period_v[i] = peak_v[i + 1] - peak_v[i];
        period_d[i] = peak_d[i + 1] - peak_d[i];
        phase_gap[i] = peak_v[i] - peak_d[i];
        ofile << peak_v[i] << '\t' << peak_d[i] << '\t' << period_v[i] << '\t' << period_d[i] << '\t' << phase_gap[i] << endl;
    }
    ofile.close();
}