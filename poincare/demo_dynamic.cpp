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
#include <pthread.h>

#include "poincare.h"

using namespace std;
using namespace std::chrono;
pthread_rwlock_t rwlock;

double theta[2] = {0};
double daylength = 12.0;
int light;
double K = 0.005;
double K_f = 0;
double omega_1 = 2 * PI / 24.0;
double omega_2 = 2 * PI / 24.1;
//double omega = 2 * PI / 4.0;
double t0 = 2000.0;

void synchronization(double, double y[], double &R, int &Ri)
{
    Ri++;
    double ave = (y[0] + y[1]) / 2.0;
    R += (abs(y[0] - ave) + abs(y[1] - ave)) / 2.0;
}

int kuramoto(double t, const double y[], double f[], void *params)
{
    double omega = *(double *)params;
    if (fmod(t, 24.0) < daylength)
        light = 1;
    else
        light = 0;

    f[0] = gam * y[0] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) - y[2] * omega_1 + K * (1 + sin(omega * t)) * (y[0] + y[1]) / 2.0 + K_f * light;
    f[1] = gam * y[1] * (a - sqrt(pow(y[1], 2) + pow(y[3], 2))) - y[3] * omega_2 + K * (1 + sin(omega * t)) * (y[0] + y[1]) / 2.0;

    f[2] = gam * y[2] * (a - sqrt(pow(y[0], 2) + pow(y[2], 2))) + y[0] * omega_1;
    f[3] = gam * y[3] * (a - sqrt(pow(y[1], 2) + pow(y[3], 2))) + y[1] * omega_2;

    theta[0] = atan2(y[2], y[0]);
    theta[1] = atan2(y[3], y[1]);

    return GSL_SUCCESS;
}

void write_kuramoto(const double y[], const double t)
{
    ofstream ofile;
    ofile.open("demo.csv", ios::app);
    ofile << t << '\t' << theta[0] << '\t' << theta[1] << '\t' << theta[0] - theta[1] << '\t' << y[0] << '\t' << y[1] << '\t' << light << endl;
    ofile.close();
}

void write_results(const double y[], const double t, double omega, double R)
{
    ofstream ofile;
    ofile.open("demo_re.csv", ios::app);
    ofile << omega << '\t' << R << endl;
    ofile.close();
}

void *threadfunc(void *arg)
{
    double P = *(double *)arg;
    double omega = 2 * PI / P;
    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 4);

    gsl_odeiv_system sys = {kuramoto, NULL, 4, &omega};

    double t = 0.0; // t0 is the time of light off
    double h = 0.1;
    double y_err[4];
    double dydt_in[4], dydt_out[4];
    double y[4] = {1, 0.1, 0, 1};

    theta[0] = atan2(y[2], y[0]);
    theta[1] = atan2(y[3], y[1]);

    double R = 0;
    int Ri = 0;

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    while (t < t0)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        if (t > 1880)
            synchronization(t, y, R, Ri);

        //write_kuramoto(y, t);

        dydt_in[0] = dydt_out[0];
        dydt_in[1] = dydt_out[1];
        dydt_in[2] = dydt_out[2];
        dydt_in[3] = dydt_out[3];
        t += h;
    }

    R = R / double(Ri);

    pthread_rwlock_wrlock(&rwlock);
    write_results(y, t, P, R);
    pthread_rwlock_unlock(&rwlock);

    gsl_odeiv_step_free(s);

    return 0;
}

int main()
{
    auto start = high_resolution_clock::now();

    ofstream ofile;

    ofile.open("demo.csv");
    ofile.close();

    int nthread = 400;
    int same_num = 1;
    double values[nthread];
    pthread_t threads[nthread];

    pthread_rwlock_init(&rwlock, NULL);

    for (size_t i = 0; i < nthread; i++)
    {
        values[i] = floor((double)i / (double)same_num) * 0.02 + 8;
        //values[i] = pow(100, i) * 0.000001;
        //values[0] = 0;
        pthread_create(&(threads[i]), NULL, threadfunc, values + i);
    }

    for (size_t i = 0; i < nthread; i++)
        pthread_join(threads[i], NULL);

    pthread_rwlock_destroy(&rwlock);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Time taken by program: " << duration.count() << " seconds" << endl;
}