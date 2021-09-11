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
    double beta_light = *(double *)arg;
    double day = 18;
    double p_cut = 0;
    double p_vd = 0.05;
    int tid = tnum++;
    double daylength = 24;
    double pulse = 0;

    double y[Ns * N];
    //static double a_vip[N][N];
    //static double a_gaba[N][N];
    double **a_vip;
    double **a_gaba;
    a_vip = new double *[N];
    a_gaba = new double *[N];
    for (int i = 0; i < N; i++)
    {
        a_vip[i] = new double[N];
        a_gaba[i] = new double[N];
    }

    double **a_vip_static;
    double **a_gaba_static;
    a_vip_static = new double *[N];
    a_gaba_static = new double *[N];
    for (int i = 0; i < N; i++)
    {
        a_vip_static[i] = new double[N];
        a_gaba_static[i] = new double[N];
    }

    double peak[N + 3][5]; // peak[N] is ventral average, peak[N+1] is dorsal average, peak[N+2] is scn average
    double peak_now[N + 3][2];
    double peak_past[N + 3][2];
    double peak_old[N + 3][2];
    int ind[N + 3];
    double R[3] = {0, 0, 0};
    double P[3] = {0, 0, 0};
    double P_var[3] = {0, 0, 0};

    int light = 0;
    int light_on = 0, light_off = 0;

    double compute_t1 = 0, compute_t2 = 0;

    int connect_num[2] = {0, 0}; // connect_num[0]: vip; connect_num[1]:gaba

    for (int i = 0; i < (N + 3); i++)
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

    //srand(tid * (unsigned int)(time(NULL)));

    //construct_connection(a_vip, a_gaba, p_link, d);
    //construct_connection_vip(a_vip, a_gaba, p_link, d, p_vd);
    construct_connection_vd(a_vip, a_gaba, p_link, d, p_vd, tid);
    connection_number(a_vip, a_gaba, connect_num);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            a_vip_static[i][j] = a_vip[i][j];
            a_gaba_static[i][j] = a_gaba[i][j];
        }
    }
    //if (tid == 0)
    //  cout << p_vd << '/' << connect_num[0] << '/' << connect_num[1] << endl;

    //write_connect(a_vip);

    double v_sP[N],
        v_sB[N], v_mB[N];

    unsigned seed = system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    normal_distribution<double> norm(0, 1);
    for (int i = 0; i < N; i++)
    {
        v_sP[i] = v_sP0 + v_sP0 * 0.06 * norm(gen);
        v_sB[i] = v_sB0 + v_sB0 * 0.01 * norm(gen);
        v_mB[i] = v_mB0 + v_mB0 * 0.01 * norm(gen);
    }

    Params param =
        {
            a_vip,
            a_gaba,
            v_sP,
            v_sB,
            v_mB,
            beta_light,
            day,
            daylength,
            pulse};

    for (int i = 0; i < N; i++)
    {
        y[i * Ns] = 2.80;
        y[i * Ns + 1] = 2.00;
        y[i * Ns + 2] = 7.94;
        y[i * Ns + 3] = 0.40;
        y[i * Ns + 4] = 12.0;
        y[i * Ns + 5] = 0.13;
        y[i * Ns + 6] = 9.00;
        y[i * Ns + 7] = 1.26;
        y[i * Ns + 8] = 0.16;
        y[i * Ns + 9] = 0.20;
        y[i * Ns + 10] = 0.091;
        y[i * Ns + 11] = 2.41;
        y[i * Ns + 12] = 0.48;
        y[i * Ns + 13] = 1.94;
        y[i * Ns + 14] = 0.32;
        y[i * Ns + 15] = 0.05;
        y[i * Ns + 16] = 0.10;
        y[i * Ns + 17] = 0.10;
        y[i * Ns + 18] = 0.00;
        y[i * Ns + 19] = 0.12;
        y[i * Ns + 20] = 0.00;
    }

    double t = 0, t0 = 2000;

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, Ns * N);
    gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-6, 1e-3);
    gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(Ns * N);

    gsl_odeiv_system sys = {scn_firing, NULL, Ns * N, &(param)};

    double h = 0.1;

    while (t < t0)
    {
        int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t0, &h, y);
        if (status != GSL_SUCCESS)
            break;

        double delay = 0;
        //if (t > 400)
        //  delay = 10;

        if ((t > 150) && (t < 1000) && (fmod(t - delay, daylength) >= (daylength - day)))
            light = 1;
        else
            light = 0;

        /*
        //////////// dynamic network
        light_on = light;

        if (light_on > light_off)
        {
            for (int i = 0; i < N_v; i++)
                for (int j = 0; j < N_v; j++)
                {
                    a_vip[i][j] = a_vip_static[i][j];
                }
        }
        else if (light_on < light_off)
            dynamicNet(a_vip, a_gaba, a_vip_static, a_gaba_static, t, delay, p_cut, day);

        light_off = light;
*/
        compute_t1 = floor(t / 90.0);

        if (compute_t1 != compute_t2)
        {
            period_N1(P, P_var, peak);
            syn_degree_N1(R, P, peak);

            pthread_rwlock_wrlock(&rwlock);
            write_rt(t, P, P_var, R);
            pthread_rwlock_unlock(&rwlock);

            for (size_t i = 0; i < 3; i++)
            {
                R[i] = 0;
                P[i] = 0; //periods
                P_var[i] = 0;
            }

            for (int i = 0; i < (N + 3); i++)
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

        compute_t2 = floor(t / 90.0);

        find_peak(t, peak, peak_now, peak_past, peak_old, y, ind);

        if ((t > 800) && (t < 1300))
            write_MP(y, t, light);
    }

    for (int i = 0; i < N; i++)
    {
        delete[] a_vip[i];
        delete[] a_gaba[i];
        delete[] a_vip_static[i];
        delete[] a_gaba_static[i];
    }
    delete[] a_vip;
    delete[] a_gaba;
    delete[] a_vip_static;
    delete[] a_gaba_static;

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);

    return (void *)0;
}

int main(int argc, char *argv[])
{
    auto start = high_resolution_clock::now();

    ofstream ofile;
    ofile.open("mp.csv");
    ofile.close();

    //double param_in = stod(argv[1]);

    int nthread = 1;
    int same_num = 5;
    int iv;
    double values[nthread];
    pthread_t threads[nthread];

    pthread_rwlock_init(&rwlock, NULL);

    for (size_t i = 0; i < nthread; i++)
    {
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
        //if (i <= 20)
        //values[i] = floor((double)i / (double)same_num) * 0.01;
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