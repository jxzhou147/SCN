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
#include "funs.h"

using namespace std;
using namespace std::chrono;

int main()
{
    auto start = high_resolution_clock::now();

    ofstream ofile;

    ofile.open("tcycle.csv");
    ofile.close();

    double TL = 24.0;
    int light;
    double t0 = 1000.0, t1 = 2300.0;

    int peak_num = 100;
    double peak_v[peak_num] = {0};
    double peak_d[peak_num] = {0};
    double period_v[peak_num] = {0};
    double period_d[peak_num] = {0};
    double phase_gap[peak_num] = {0};
    double peak_now_v[2] = {0};
    double peak_past_v[2] = {0};
    double peak_old_v[2] = {0};
    double peak_now_d[2] = {0};
    double peak_past_d[2] = {0};
    double peak_old_d[2] = {0};
    int ind_v[1] = {0};
    int ind_d[1] = {0};

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 4);

    gsl_odeiv_system sys = {kuramoto, NULL, 4, &TL};

    double t = 0.0; // t0 is the time of light off
    double h = 0.1;
    double y_err[4];
    double dydt_in[4], dydt_out[4];
    double y[4] = {1, 0.1, 0, 1};

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);

    while (t < t1)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        if ((t < 1000) && (fmod(t, TL) < (TL / 2.0)))
            light = 1;
        else
            light = 0;

        //write_curve(y, t, light);

        if (t > (t0 - 500))
            find_peak(t, peak_v, peak_now_v, peak_past_v, peak_old_v, peak_d, peak_now_d, peak_past_d, peak_old_d, y, ind_v, ind_d);

        dydt_in[0] = dydt_out[0];
        dydt_in[1] = dydt_out[1];
        dydt_in[2] = dydt_out[2];
        dydt_in[3] = dydt_out[3];
        t += h;
    }

    write_peaks(peak_v, peak_d, period_v, period_d, phase_gap, peak_num);

    gsl_odeiv_step_free(s);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Time taken by program: " << duration.count() << " seconds" << endl;

    return 0;
}