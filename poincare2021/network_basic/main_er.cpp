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
//#include "funs.h"
#include "new_funs.h"

using namespace std;
using namespace std::chrono;

pthread_rwlock_t rwlock;
int tnum = 1;
int numEntrain = 0;
int ifEntrain = 0;

void *threadfunc(void *arg)
{
    double tcycle = *(double *)arg;
    int tid = tnum++;

    double daylength = tcycle / 2.0;
    double K_f = 0.15;

    double K = 0.1;
    double sigma = 0.05;

    double p_random = 0.01;
    //double p_random = 0.015;
    double d = 1;
    double p_link = 0.01;
    double p_vd = 0.01;
    double p_dv = 0;
    double p_cut = 0.1;
    double t_re = 1;

    double y[2 * N];
    double ed[N];
    double eta[N];
    double **a;
    a = new double *[N];
    for (int i = 0; i < N; i++)
    {
        a[i] = new double[N];
    }

    double **a_static;
    a_static = new double *[N];
    for (int i = 0; i < N; i++)
    {
        a_static[i] = new double[N];
    }

    /*
    int peak_ave_num = 100;
    double peak_ave_v[peak_ave_num] = {0};
    double peak_now_v[2] = {0};
    double peak_past_v[2] = {0};
    double peak_old_v[2] = {0};
    int ind_ave_v[1] = {0};
    double peak_ave_d[peak_ave_num] = {0};
    double peak_now_d[2] = {0};
    double peak_past_d[2] = {0};
    double peak_old_d[2] = {0};
    int ind_ave_d[1] = {0};
    */

    double peak[N + 1][5] = {0};
    double peak_now[N + 1][2] = {0};
    double peak_past[N + 1][2] = {0};
    double peak_old[N + 1][2] = {0};
    int ind[N + 1] = {0};

    int light = 0;
    int light_on = 0, light_off = 0;
    int t1, t2;
    double compute_t1 = 0, compute_t2 = 0;

    //srand((unsigned)time(NULL));

    construct_connection_random(a, p_random, tid);
    //construct_connection_vd(a, p_link, d, p_vd, p_dv, tid);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            a_static[i][j] = a[i][j];
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        ed[i] = 0;
        for (size_t j = 0; j < N; j++)
            ed[i] += fabs(a[j][i]);
    }

    unsigned seed = system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    normal_distribution<double> norm{1, sigma};
    for (size_t i = 0; i < N; i++)
    {
        eta[i] = norm(gen);
    }

    for (size_t i = 0; i < N; i++)
    {
        y[i] = (rand() % (N_rand + 1) / (float)(N_rand + 1)) * 2.0 - 1.0;
        y[i + N] = (rand() % (N_rand + 1) / (float)(N_rand + 1)) * 2.0 - 1.0;
    }

    Params param =
        {
            a,
            ed,
            eta,
            daylength,
            tcycle,
            K,
            K_f};

    double t = 0, t0 = 1000;

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 2 * N);

    gsl_odeiv_system sys = {poincare, NULL, 2 * N, &(param)};

    double h = 0.1;
    double y_err[2 * N];
    double dydt_in[2 * N], dydt_out[2 * N];

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);

    while (t < t0)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        if (fmod(t, tcycle) < daylength)
            light = 1;
        else
            light = 0;

        //if (t > 600)
        //  light = 1;

        //////////// dynamic network
        /*
        t1 = floor(t / t_re);
        if (t1 != t2)
        {
            dynamicWireAndCut(a, p_cut, tid);

            for (size_t i = 0; i < N; i++)
            {
                ed[i] = 0;
                for (size_t j = 0; j < N; j++)
                    ed[i] += fabs(a[j][i]);
            }
        }
        t2 = floor(t / t_re);
*/
        light_on = light;

        if (light_on > light_off)
        {
            for (int i = 0; i < N_v; i++)
                for (int j = 0; j < N_v; j++)
                {
                    a[i][j] = a_static[i][j];
                }

            for (size_t i = 0; i < N; i++)
            {
                ed[i] = 0;
                for (size_t j = 0; j < N; j++)
                    ed[i] += fabs(a[j][i]);
            }
        }
        else if (light_on < light_off)
        {
            dynamicNet(a, p_cut);
            for (size_t i = 0; i < N; i++)
            {
                ed[i] = 0;
                for (size_t j = 0; j < N; j++)
                    ed[i] += fabs(a[j][i]);
            }
        }

        light_off = light;
        /*
        compute_t1 = floor(t / 150.0);

        if (compute_t1 != compute_t2)
        {
            double period_ = period(peak);

            pthread_rwlock_wrlock(&rwlock);
            write_rt(t, period_, syn_degree(period_, peak));
            pthread_rwlock_unlock(&rwlock);

            for (int i = 0; i < (N + 1); i++)
            {
                for (int j = 0; j < 5; j++)
                    peak[i][j] = 0;
                for (int j = 0; j < 2; j++)
                {
                    peak_now[i][j] = 10;
                    peak_past[i][j] = 10;
                    peak_old[i][j] = 10;
                }
                ind[i] = 0;
            }
        }

        compute_t2 = floor(t / 150.0);
        */

        if (t > (t0 - 300))
        {
            find_peak(t, peak, peak_now, peak_past, peak_old, y, ind);
            write_curve(y, t, light);
        }

        for (size_t i = 0; i < 2 * N; i++)
            dydt_in[i] = dydt_out[i];
        t += h;
    }

    double period_ = period(peak);

    pthread_rwlock_wrlock(&rwlock);
    numEntrain++;
    if (fabs(period_ - tcycle) < 0.1)
        ifEntrain++;
    //ifEntrain += ifEntrainFun(peak, tcycle);
    pthread_rwlock_unlock(&rwlock);

    //if (numEntrain == 30)
    write_er(tcycle, K_f, period_, ifEntrain);

    /*
    pthread_rwlock_wrlock(&rwlock);
    write_peaks_vd_p(peak_ave_v, peak_ave_d, peak_ave_num, p_vd, p_dv, daylength);
    pthread_rwlock_unlock(&rwlock);
    //write_peaks(peak_ave, peak_ave_num);
*/
    for (int i = 0; i < N; i++)
    {
        delete[] a[i];
        delete[] a_static[i];
    }
    delete[] a;
    delete[] a_static;

    gsl_odeiv_step_free(s);

    return (void *)0;
}

int main(int argc, char *argv[])
{
    auto start = high_resolution_clock::now();

    ofstream ofile;
    ofile.open("curve.csv");
    ofile.close();

    //double param_in = stod(argv[1]);

    int nthread = 1;
    int same_num = 30;
    int iv;
    double values[nthread];
    pthread_t threads[nthread];

    pthread_rwlock_init(&rwlock, NULL);

    for (size_t i = 0; i < nthread; i++)
    {
        iv = i / same_num;
        //if (i <= 20)
        values[i] = floor((double)i / (double)same_num) * 0.1 + 21.9;
        //else
        //   values[i] = (floor((double)i / (double)same_num) - 20) * 0.1;
        //values[i] = pow(10, i) * 0.1;
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