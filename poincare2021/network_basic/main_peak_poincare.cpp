#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <random>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "parameter.h"
#include "funs.h"

using namespace std;
using namespace std::chrono;

int main()
{
    auto start = high_resolution_clock::now();
    ofstream ofile;
    ofile.open("curve.csv");
    ofile.close();

    double daylength = 18;

    double K = 0.15;
    double K_f = 0.1;
    double sigma = 0.05;

    double d = 1;
    double p_link = 0.1;
    double p_vd = 0.01;
    double p_dv = 0.0;

    double y[2 * N];
    double ed[N];
    double eta[N];
    double **a;
    a = new double *[N];
    for (int i = 0; i < N; i++)
    {
        a[i] = new double[N];
    }

    int peak_ave_num = 350;
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

    double peak_ave[peak_ave_num] = {0};
    double peak_now[2] = {0};
    double peak_past[2] = {0};
    double peak_old[2] = {0};
    int ind_ave[1] = {0};

    int light = 0;
    int tid = 0;

    int connect_num[2] = {0, 0}; // connect_num[0]: vip; connect_num[1]:gaba

    //srand((unsigned)time(NULL));

    construct_connection_vd(a, p_link, d, p_vd, p_dv, tid);

    for (size_t i = 0; i < N; i++)
    {
        ed[i] = 0;
        for (size_t j = 0; j < N; j++)
            ed[i] += abs(a[j][i]);
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
            K,
            K_f};

    double t = 0, t0 = 8000;

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

        if (fmod(t, 24.0) >= (24 - daylength))
            light = 1;
        else
            light = 0;

        //if (t > 600)
        //  light = 1;

        //////////// dynamic network
        //light_on = light;
        /*
        if (light_on > light_off)
        {
            for (int i = 0; i < N_v; i++)
                for (int j = 0; j < N; j++)
                {
                    a_vip[i][j] = a_vip_static[i][j];
                }
        }
        else if (light_on < light_off)
            dynamicNet(a_vip, a_gaba, a_vip_static, a_gaba_static, t, delay, p_cut, day);
            */
        /*
        {
            for (int i = 0; i < N_v; i++)
            {
                for (int j = 0; j < N_v; j++)
                {
                    if ((a_vip[i][j] == 1))
                        a_vip[i][j] = 0;
                }
            }
        }*/
        /*
        //// vtod 1/0
        if (light_on > light_off)
        {
            dynamicNet(a_vip, a_gaba, a_vip_static, a_gaba_static, t, delay, p_cut, day);
        }
        else if (light_on < light_off)
        {

            for (int i = 0; i < N_v; i++)
                for (int j = 0; j < N; j++)
                {
                    a_vip[i][j] = a_vip_static[i][j];
                }
        }
*/
        //light_off = light;

        /*
        if ((t > 600) & (t < 610))
        {
            for (int i = 0; i < N_v; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    a_vip[i][j] = 0;
                }
            }
        }

        if ((t > 800) & (t < 810))
        {
            for (int i = 0; i < N_v; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    a_vip[i][j] = a_vip_static[i][j];
                }
            }
        }
*/
        if (t > 300)
        {
            find_peak_ave_vd(t, peak_ave_v, peak_now_v, peak_past_v, peak_old_v, peak_ave_d, peak_now_d, peak_past_d, peak_old_d, y, ind_ave_v, ind_ave_d);
            //            write_curve(y, t, light);
        }

        for (size_t i = 0; i < 2 * N; i++)
            dydt_in[i] = dydt_out[i];
        t += h;
    }

    write_peaks_vd(peak_ave_v, peak_ave_d, peak_ave_num);
    //write_peaks(peak_ave, peak_ave_num);

    for (int i = 0; i < N; i++)
    {
        delete[] a[i];
    }
    delete[] a;

    gsl_odeiv_step_free(s);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << "Time taken by program: " << duration.count() << " seconds" << endl;

    return 0;
}