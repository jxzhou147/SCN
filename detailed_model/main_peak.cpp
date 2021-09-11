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
    ofile.open("mp.csv");
    ofile.close();
    ofile.open("firingRate.csv");
    ofile.close();

    double day = 18;
    double p_cut = 0.3;
    double p_vd = 0.05;
    //double p_vd = 0.5;
    double beta_light = 0.8;
    double daylength = 24;

    double y[Ns * N];
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

    int peak_ave_num = 40;
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

    double pulse = 0;
    int light = 0;
    int light_on = 0, light_off = 0;

    int connect_num[2] = {0, 0}; // connect_num[0]: vip; connect_num[1]:gaba

    //srand((unsigned)time(NULL));

    construct_connection_vd(a_vip, a_gaba, p_link, d, p_vd, 1);
    /*
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
        {
            a_vip[i][j] = 0;
            a_gaba[i][j] = 0;
        }
        */
    /*
    ////// static 0.5/0.5
    for (int i = 0; i < N_v; i++)
    {
        for (int j = 0; j < N_v; j++)
        {
            if ((a_vip[i][j] == 1) && ((rand() / double(RAND_MAX)) < p_cut))
                a_vip[i][j] = 0;
        }
    }
*/
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            a_vip_static[i][j] = a_vip[i][j];
            a_gaba_static[i][j] = a_gaba[i][j];
        }
    }

    //write_connect(a_vip);

    double v_sP[N],
        v_sB[N], v_mB[N];

    unsigned seed = system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    normal_distribution<double> norm(0, 1);

    for (int i = 0; i < N_v; i++)
    {
        v_sP[i] = v_sP0 + v_sP0 * 0.01 * norm(gen);
        v_sB[i] = v_sB0 + v_sB0 * 0.01 * norm(gen);
        v_mB[i] = v_mB0 + v_mB0 * 0.01 * norm(gen);
    }

    for (int i = N_v; i < N; i++)
    {
        v_sP[i] = v_sP0_d + v_sP0_d * 0.01 * norm(gen);
        v_sB[i] = v_sB0 + v_sB0 * 0.01 * norm(gen);
        v_mB[i] = v_mB0 + v_mB0 * 0.01 * norm(gen);
    }
    /*
    for (int i = 0; i < N_v; i++)
    {
        v_sP[i] = v_sP0;
        v_sB[i] = v_sB0;
        v_mB[i] = v_mB0;
    }

    for (int i = N_v; i < N; i++)
    {
        v_sP[i] = v_sP0_d;
        v_sB[i] = v_sB0;
        v_mB[i] = v_mB0;
    }
    */
    /*
    for (int i = 0; i < N_v; i++)
    {
        v_sP[i] = v_sP0 + v_sP0 * 0.06 * norm(gen);
        v_sB[i] = v_sB0 + v_sB0 * 0.01 * norm(gen);
        v_mB[i] = v_mB0 + v_mB0 * 0.01 * norm(gen);
    }

    for (int i = N_v; i < N; i++)
    {
        v_sP[i] = v_sP0_d + v_sP0_d * 0.06 * norm(gen);
        v_sB[i] = v_sB0 + v_sB0 * 0.01 * norm(gen);
        v_mB[i] = v_mB0 + v_mB0 * 0.01 * norm(gen);
    }
*/
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

    double t = 0, t0 = 1000;

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
        /*
        for (int i = 0; i < N; i++)
        {
            if (y[i * Ns + 18] < 1e-7)
                y[i * Ns + 18] = 1e-7;
            if (y[i * Ns + 20] < 1e-7)
                y[i * Ns + 20] = 1e-7;
        }
*/
        //double pulse = 500;
        double delay = 0;

        //if (t > 396)
        //  delay = 6;

        //if ((t > 150) & (fmod(t - delay, 24.0) >= (12.0)))
        //if ((t > 480) & (fmod(t - delay, 24.0) <= (12.0)))

        if ((t > 150) && (fmod(t - delay, 24.0) >= (24 - day)))
            light = 1;
        else
            light = 0;

        //if (t > 600)
        //  light = 1;

        //////////// dynamic network
        light_on = light;
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
        light_off = light;

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
        //if (t > 300)
        //find_peak_ave(t, peak_ave, peak_now, peak_past, peak_old, y, ind_ave);

        if (t > 700)
        {
            find_peak_ave_vd(t, peak_ave_v, peak_now_v, peak_past_v, peak_old_v, peak_ave_d, peak_now_d, peak_past_d, peak_old_d, y, ind_ave_v, ind_ave_d);
            write_MP(y, t, light);
        }
    }

    write_peaks_vd(peak_ave_v, peak_ave_d, peak_ave_num);
    //write_peaks(peak_ave, peak_ave_num);

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

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << "Time taken by program: " << duration.count() << " seconds" << endl;

    return 0;
}