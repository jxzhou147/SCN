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

#include "poincare.h"

using namespace std;
using namespace std::chrono;

double theta[2] = {0};
double daylength = 12.0;
int light;
//double K = 0.1;
double Kvd = -0.01;
double Kdv = 0.1;
double K_f = 0.2;
double t0 = 1000.0, t1 = 2000.0;

double phase_light(double pha_pi, double light)
{
    double t0 = 2;
    double pha;
    if (pha_pi > 0)
        pha = pha_pi * 24.0 / (2.0 * PI);
    else
        pha = (pha_pi + 2 * PI) * 24.0 / (2.0 * PI);

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
    if (t > t0)
        daylength = 0;

    if (fmod(t, 24.0) < daylength)
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

    f[0] = gam * y[0] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) - y[2] * omega_v + Kdv * y[1] + K_f * light;
    f[1] = gam * y[1] * (a - sqrt(pow(y[1], 2) + pow(y[3], 2))) - y[3] * omega_d + Kvd * y[0];

    f[2] = gam * y[2] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) + y[0] * omega_v;
    f[3] = gam * y[3] * (a - sqrt(pow(y[1], 2) + pow(y[3], 2))) + y[1] * omega_d;

    theta[0] = atan2(y[2], y[0]);
    theta[1] = atan2(y[3], y[1]);

    return GSL_SUCCESS;
}

void write_kuramoto(const double y[], const double t)
{
    double phase_sun = fmod(Omega * t, 2 * PI);

    ofstream ofile;
    ofile.open("demo.csv", ios::app);
    ofile << t << '\t' << phase_sun << '\t' << theta[0] << '\t' << theta[1] << '\t' << theta[0] - theta[1] << '\t' << y[0] << '\t' << y[1] << '\t' << light << endl;
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

    double t = 0.0; // t0 is the time of light off
    double h = 0.1;
    double y_err[4];
    double dydt_in[4], dydt_out[4];
    double y[4] = {1, 0.1, 0, 1};

    theta[0] = atan2(y[2], y[0]);
    theta[1] = atan2(y[3], y[1]);

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    while (t < t0)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        write_kuramoto(y, t);

        dydt_in[0] = dydt_out[0];
        dydt_in[1] = dydt_out[1];
        t += h;
    }

    daylength = 0;

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