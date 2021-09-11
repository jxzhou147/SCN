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

double period = 24;
//double omega = 2 * PI / 3;
double omega = 0;
double theta = 0;
double daylength = 12.0;
int light;
double K_f = 0;
double t0 = 500.0;

double peak[5] = {0};
double peak1, peak2;
int ind = 0;

double Period()
{
    double P = (peak[4] - peak[0]) / 4;
    return P;
}

void find_peak(const double t)
{
    if (ind < 5)
    {
        peak1 = theta;

        if (((peak2 - peak1) > 6) & (ind < 5))
        {
            peak[ind] = t;
            ind++;
        }

        peak2 = theta;
    }
}

int kuramoto(double t, const double y[], double f[], void *params)
{
    if (fmod(t, period) < daylength)
        light = 1;
    else
        light = 0;

    f[0] = gam * y[0] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) - y[2] * omega + K_f * light;

    f[2] = gam * y[2] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) + y[0] * omega;

    f[1] = 0;

    theta = atan2(y[2], y[0]);

    return GSL_SUCCESS;
}

void write_kuramoto(const double y[], const double t)
{
    double phase_sun = fmod(2 * PI * t / period, 2 * PI);

    ofstream ofile;
    ofile.open("demo.csv", ios::app);
    //ofile << t << '\t' << phase_sun << '\t' << y[0] << '\t' << y[2] << '\t' << theta << '\t' << light << endl;
    ofile << t << '\t' << theta << endl;
    ofile.close();
}

void write_results(const double t)
{
    ofstream ofile;
    ofile.open("simple.csv", ios::app);

    ofile << t << '\t' << K_f << '\t' << omega << '\t' << period << '\t' << Period() << endl;

    ofile.close();
}

int main()
{
    auto start = high_resolution_clock::now();

    ofstream ofile;

    ofile.open("demo.csv");
    ofile.close();

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 3);

    gsl_odeiv_system sys = {kuramoto, NULL, 3, NULL};

    double t = 0.0; // t0 is the time of light off
    double h = 0.1;
    double y_err[3];
    double dydt_in[3], dydt_out[3];
    double y[3] = {1, 0.1, 0};

    theta = atan2(y[2], y[0]);

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    while (t < t0)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        write_kuramoto(y, t);

        if (t > (t0 - 200))
        {
            find_peak(t);
        }

        dydt_in[0] = dydt_out[0];
        dydt_in[1] = dydt_out[1];
        dydt_in[2] = dydt_out[2];
        t += h;
    }

    write_results(t);

    gsl_odeiv_step_free(s);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Time taken by program: " << duration.count() << " seconds" << endl;

    return 0;
}