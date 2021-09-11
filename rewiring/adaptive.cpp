#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "kuramoto.h"

using namespace std;
using namespace std::chrono;

double theta[2] = {0};
double daylength = 18.0;
int light;
//double k1 = 0.01;
//double k2 = 0.01;
//double K_f = 0.05;
double Kdv = 0.05;
double Kvd = -0.01;
double p = 0;
double p_c = 0.9;
double epsilon = 1.0;

double peak1[2] = {0};
double peak2[2] = {0};
double P[2] = {0, 0};
double peak[2][5];
int ind[2];

double delta_theta(const double (*peak)[5])
{
    double tv[5] = {0, 0, 0, 0, 0};
    double ts[5] = {0, 0, 0, 0, 0};
    double tvs[5] = {0, 0, 0, 0, 0};
    double tvsi;
    double delta_the = 0;

    for (size_t i = 0; i < 5; i++)
    {
        tvs[i] = abs((peak[0][i] - peak[1][i]));
        if (tvs[i] > 12)
            tvs[i] = abs(24 - tvs[i]);

        delta_the += tvs[i];
    }

    delta_the = delta_the / 5;

    return delta_the;
}

void period(double P[], const double (*peak)[5])
{
    P[0] = ((peak[0][4] - peak[0][0]) / 4);

    P[1] = (peak[1][4] - peak[1][0]) / 4;
}

void find_peak(const double t, double (*peak)[5], double peak1[], double peak2[], const double theta[], int ind[])
{
    int ifpeak = 1;
    for (size_t i = 0; i < 2; i++)
    {
        if (theta[i] > 0)
            ifpeak = 0;
    }

    //if (ifpeak == 1)
    //{
    int ind_min = 5;
    for (size_t i_min = 0; i_min < 2; i_min++)
    {
        if (ind_min > ind[i_min])
            ind_min = ind[i_min];
    }

    if (ind_min < 5)
    {
        for (size_t i = 0; i < 2; i++)
            peak1[i] = theta[i];

        for (size_t i_p = 0; i_p < 2; i_p++)
        {
            if (((peak2[i_p] - peak1[i_p]) > 6) & (ind[i_p] < 5))
            {
                peak[i_p][ind[i_p]] = t;
                //peak(i_p, ind[i_p]) = peak2[i_p];
                ind[i_p] = ind[i_p] + 1;
            }
        }

        for (size_t i = 0; i < 2; i++)
            peak2[i] = theta[i];
    }
    //}
}

double phase(double pha_pi)
{
    //double pha = 2.0 * PI / 24.0;
    double pha = pha_pi * 24.0 / (2.0 * PI);
    if ((pha > 3.0) & (pha < 9.0))
        return 0.0;
    else
    {
        if (pha > 9.0)
            return -1.0 * sin((pha - 9) * 2 * PI / 18.0);
        else
            return -1.0 * sin((pha + 24 - 9) * 2 * PI / 18.0);
    }
}

double phase_light(double pha_pi, double light)
{
    double t0 = 2;
    double pha = pha_pi * 24.0 / (2.0 * PI);
    if ((pha > t0) & (pha < (light - t0))) //// keep the morning 3 hours and dusk 3 hours having effect of entrainment
        return 0.0;
    else
    {
        if (pha > (light - t0))
            return -1 * sin((pha - (light - t0)) * 2 * PI / (24 - (light - 2 * t0)));
        else
            return -1 * sin((pha + 24 - (light - t0)) * 2 * PI / (24 - (light - 2 * t0)));
    }
}

int kuramoto(double t, const double y[], double f[], void *params)
{
    if (fmod(t, 24.0) < daylength)
        light = 1;
    else
        light = 0;

    theta[0] = fmod(y[0], 2 * PI);
    theta[1] = fmod(y[1], 2 * PI);
    p = sqrt((1 + cos(theta[0] - theta[1])) / 2);

    f[0] = omega_v + K_f * light * phase(theta[0]) + Kdv * y[2] * sin(y[1] - y[0]);
    //f[0] = omega_v + K_f * light * phase(theta[0]) + Kdv * sin(y[1] - y[0]);
    //f[0] = omega_v + K_f * light * phase_light(theta[0], daylength) + k1 * sin(y[1] - y[0]);
    //f[0] = omega_v + K_f * light * phase(theta[0]);
    //f[0] = omega_v + k1 * sin(y[1] - y[0]) + light * 0.05;
    //f[0] = omega_v + light * 0.5;
    //f[0] = omega_v + K_f * light * phase_light(theta[0], daylength);
    //f[1] = omega_d + Kvd * sin(y[0] - y[1]);
    f[1] = omega_d + Kvd * y[3] * sin(y[0] - y[1]);
    f[2] = epsilon * (p - p_c) * y[2] * (1 - y[2]);
    f[3] = epsilon * (p - p_c) * y[3] * (1 - y[3]);

    return GSL_SUCCESS;
}

void write_kuramoto(const double y[], const double t)
{
    double phase_sun = fmod(Omega * t, 2 * PI);

    ofstream ofile;
    ofile.open("demo.csv", ios::app);
    ofile << t << '\t' << y[2] << '\t' << y[3] << '\t' << y[0] << '\t' << y[1] << '\t' << light << '\t' << fmod(y[0] - y[1], 2 * PI) << '\t' << fmod(y[0] - Omega * t, 2 * PI) << endl;
    ofile.close();
}

void write_results()
{
    ofstream ofile;
    ofile.open("adaptive.csv", ios::app);
    ofile << daylength << '\t' << Kdv << '\t' << Kvd << '\t' << p_c << '\t' << epsilon << '\t' << delta_theta(peak) << '\t' << P[0] << '\t' << P[1] << endl;
    ofile.close();
}

int main()
{
    auto start = high_resolution_clock::now();

    ofstream ofile;

    ofile.open("demo.csv");
    ofile.close();

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 4);

    gsl_odeiv_system sys = {kuramoto, NULL, 4, NULL};

    double t = 0.0, t0 = 800.0, t1 = 100.0; // t0 is the time of light off
    double h = 0.1;
    double y_err[4];
    double dydt_in[4], dydt_out[4];
    double y[4] = {1, 2, 0.5, 0.5};

    theta[0] = fmod(y[0], 2 * PI);
    theta[1] = fmod(y[1], 2 * PI);

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    while (t < t0)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        if (t > (t0 - 200))
            find_peak(t, peak, peak1, peak2, theta, ind);

        write_kuramoto(y, t);

        dydt_in[0] = dydt_out[0];
        dydt_in[1] = dydt_out[1];
        t += h;
    }

    period(P, peak);
    write_results();

    daylength = 12;

    while (t < t1)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        write_kuramoto(y, t);

        dydt_in[0] = dydt_out[0];
        dydt_in[1] = dydt_out[1];
        dydt_in[2] = dydt_out[2];
        dydt_in[3] = dydt_out[3];
        t += h;
    }

    gsl_odeiv_step_free(s);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Time taken by program: " << duration.count() << " seconds" << endl;

    return 0;
}