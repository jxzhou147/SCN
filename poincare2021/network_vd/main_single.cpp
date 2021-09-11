#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <random>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <pthread.h>

#include "parameter.h"
#include "funs.h"

using namespace std;
using namespace std::chrono;

pthread_rwlock_t rwlock;
int tnum = 0;

void *threadfunc(void *arg)
{
    double day = *(double *)arg;
    double beta_light = 0.8;
    int tid = tnum++;
    double daylength = 24;

    double y[Ns];

    double peak[N + 3][5]; // peak[N] is ventral average, peak[N+1] is dorsal average, peak[N+2] is scn average
    double peak_now[N + 3][2];
    double peak_past[N + 3][2];
    double peak_old[N + 3][2];
    int ind[N + 3];
    double R[3] = {0, 0, 0};
    double P[3] = {0, 0, 0};
    double P_var[3] = {0, 0, 0};

    int light = 0;

    for (int i = 0; i < (N + 3); i++)
    {
        for (int j = 0; j < 5; j++)
            peak[i][j] = 0;
        for (int j = 0; j < 2; j++)
        {
            peak_now[i][j] = 0;
            peak_past[i][j] = 0;
            peak_old[i][j] = 0;
        }
        ind[i] = 0;
    }

    for (int i = 0; i < N; i++)
    {
        y[0] = 2.80;
        y[1] = 2.00;
        y[2] = 7.94;
        y[3] = 0.40;
        y[4] = 12.0;
        y[5] = 0.13;
        y[6] = 9.00;
        y[7] = 1.26;
        y[8] = 0.16;
        y[9] = 0.20;
        y[10] = 0.091;
        y[11] = 2.41;
        y[12] = 0.48;
        y[13] = 1.94;
        y[14] = 0.32;
        y[15] = 0.05;
        y[16] = 0.10;
        y[17] = 0.10;
        y[18] = 0.00;
        y[19] = 0.12;
        y[20] = 0.00;
    }

    double t = 0, t0 = 600;

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, Ns);
    gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-6, 1e-3);
    gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(Ns);

    gsl_odeiv_system sys = {scn_firing, NULL, Ns, &(daylength)};

    double h = 0.1;

    while (t < t0)
    {
        int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t0, &h, y);
        if (status != GSL_SUCCESS)
            break;

        if (fmod(t, 24.0) >= (24.0 - day))
            light = 1;
        else
            light = 0;

        write_n1(y, t, light);
        /*
        if (t > 450)
        {
            find_peak(t, peak, peak_now, peak_past, peak_old, y, ind);
            //write_MP(y, t, light);
        }
        */
    }

    //period(P, P_var, peak);
    //syn_degree(R, P, peak);

    //write_results(P, P_var, R, p_vd, connect_num, beta_light, p_cut, day);

    //write_period_N(peak);

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);

    return (void *)0;
}

int main(int argc, char *argv[])
{
    auto start = high_resolution_clock::now();

    ofstream ofile;
    ofile.open("n1.csv");
    ofile.close();

    //double param_in = stod(argv[1]);

    int nthread = 1;
    int same_num = 1;
    int iv;
    double values[nthread];
    pthread_t threads[nthread];

    pthread_rwlock_init(&rwlock, NULL);

    for (size_t i = 0; i < nthread; i++)
    {
        /*
        iv = i / same_num;
        switch (iv)
        {
        case 0:
            values[i] = 0.8;
            break;
        case 1:
            values[i] = 0.6;
            break;
        case 2:
            values[i] = 0.7;
            break;
        case 3:
            values[i] = 0.8;
            break;
        case 4:
            values[i] = 0.9;
            break;
        case 5:
            values[i] = 1;
            break;
        }
        */
        //if (i <= 20)
        values[i] = floor((double)i / (double)same_num) * 0.005 + 18;
        //else
        //   values[i] = (floor((double)i / (double)same_num) - 20) * 0.1;
        //values[i] = pow(10, i) * 0.1;
        //values[i] = 0.8;
        pthread_create(&(threads[i]), NULL, threadfunc, values + i);
    }

    for (size_t i = 0; i < nthread; i++)
        pthread_join(threads[i], NULL);

    pthread_rwlock_destroy(&rwlock);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << "Time taken by program: " << duration.count() << " seconds" << endl;

    return 0;
}