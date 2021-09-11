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

using namespace std;
using namespace std::chrono;
pthread_rwlock_t rwlock;

double eps0 = 1.6;
double Omega = 1.0;

void synchronization(double, double y[], double &R, int &Ri)
{
    Ri++;
    double ave = (y[0] + y[3] + y[6]) / 3.0;
    R += (abs(y[0] - ave) + abs(y[3] - ave) + abs(y[6] - ave)) / 3.0;
}

int kuramoto(double t, const double y[], double f[], void *params)
{
    double omega = *(double *)params;

    f[0] = -1 * Omega * y[1] - y[2] + eps0 * (1 + sin(omega * t)) * (y[3] + y[6] - 2 * y[0]);
    f[1] = Omega * y[0] + 0.2 * y[1];
    f[2] = 0.2 + y[0] * y[2] - 9 * y[2];

    f[3] = -1 * (Omega + 0.001) * y[4] - y[5] + eps0 * (1 + sin(omega * t)) * (y[0] + y[6] - 2 * y[3]);
    f[4] = Omega * y[3] + 0.2 * y[4];
    f[5] = 0.2 + y[3] * y[5] - 9 * y[5];

    f[6] = -1 * (Omega - 0.001) * y[7] - y[8] + eps0 * (1 + sin(omega * t)) * (y[0] + y[3] - 2 * y[6]);
    f[7] = Omega * y[6] + 0.2 * y[7];
    f[8] = 0.2 + y[6] * y[8] - 9 * y[8];

    return GSL_SUCCESS;
}

void write_curve(const double y[], const double t)
{
    ofstream ofile;
    ofile.open("demo_rossler.csv", ios::app);
    ofile << t << '\t' << y[0] << '\t' << y[3] << '\t' << y[6] << endl;
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
    double omega = *(double *)arg;
    double t0 = 1000.0;

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 9);

    gsl_odeiv_system sys = {kuramoto, NULL, 9, &omega};

    double t = 0.0; // t0 is the time of light off
    double h = 0.1;
    double y_err[9];
    double dydt_in[9], dydt_out[9];
    double y[9] = {1, 0.1, 0, 1, 0, 1, 4, 5, 6};

    double R = 0;
    int Ri = 0;

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    while (t < t0)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        if (t > 900)
            synchronization(t, y, R, Ri);

        //write_curve(y, t);

        for (int i = 0; i < 9; i++)
        {
            dydt_in[i] = dydt_out[i];
        }
        t += h;
    }

    R = R / double(Ri);

    pthread_rwlock_wrlock(&rwlock);
    write_results(y, t, omega, R);
    pthread_rwlock_unlock(&rwlock);

    gsl_odeiv_step_free(s);

    return 0;
}

int main()
{
    auto start = high_resolution_clock::now();

    ofstream ofile;

    ofile.open("demo_rossler.csv");
    ofile.close();

    int nthread = 2000;
    int same_num = 1;
    double values[nthread];
    pthread_t threads[nthread];

    pthread_rwlock_init(&rwlock, NULL);

    for (size_t i = 0; i < nthread; i++)
    {
        values[i] = floor((double)i / (double)same_num) * 0.001;
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