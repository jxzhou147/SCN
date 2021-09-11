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

    double y[Ns * N];
    double a_vip[N][N];
    double a_gaba[N][N];

    construct_connection(a_vip, a_gaba, p_link, d);
    /*for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            a_vip[i][j] =0;
            a_gaba[i][j] = 0;
        }*/
    Params param =
    {
        a_vip,
        a_gaba
    };

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
        y[i * Ns + 9]= 0.20;
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

    double t = 0, t0 = 600;

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, Ns * N);

    gsl_odeiv_system sys ={ scn_firing, NULL, Ns * N, &(param) };

    double h = 0.1;
    double y_err[Ns * N];
    double dydt_in[Ns * N], dydt_out[Ns * N];

    //for (int i = 0; i < Ns * N; i++)
      //  dydt_in[i] = 0;
    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    while (t < t0)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        write_MP(y, t);

        for (size_t i = 0; i < Ns * N; i++)
            dydt_in[i] = dydt_out[i];
        t += h;
    }

    gsl_odeiv_step_free(s);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Time taken by program: " << duration.count() << " seconds" << endl;

    return 0;
}