#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "kuramoto.h"
#include "fun.h"

using namespace std;

int main(void)
{
    ofstream ofile;

    ofile.open("kuramoto.csv");
    ofile.close();
    /*
	ofile.open("peaks.csv");
	ofile.close();

	ofile.open("connection.csv");
	ofile.close();
    */

    Kuramoto kuramoto;

    kuramoto.construct_connection();
    kuramoto.initialization();

    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, N);

    gsl_odeiv_system sys = {func, NULL, N, &kuramoto};

    double t = 0.0, t0 = 500.0;
    double h = 0.1;
    double y_err[N];
    double dydt_in[N], dydt_out[N];
    double y[N];

    for (size_t i = 0; i < N; i++)
    {
        y[i] = 2 * PI * (rand() % (N_rand + 1) / (float)(N_rand + 1));
    }

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    while (t < t0)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        if (t > 350)
        {
            if (kuramoto.times < 100)
            {
                kuramoto.times++;
                kuramoto.syn_degree();
            }
        }

        ///////////rewiring
        kuramoto.t1 = floor(t / t_re);
        if (kuramoto.t1 != kuramoto.t2)
            kuramoto.rewiring();
        kuramoto.t2 = floor(t / t_re);

        if (t > 300)
            kuramoto.find_peak(t);

        kuramoto.write_kuramoto(y, t);

        dydt_in[0] = dydt_out[0];
        dydt_in[1] = dydt_out[1];
        t += h;
    }

    kuramoto.R[0] = kuramoto.R_[0] / (N_v * kuramoto.times);
    kuramoto.R[1] = kuramoto.R_[1] / ((N_v - N_v) * kuramoto.times);
    kuramoto.R[2] = kuramoto.R_[2] / (N * kuramoto.times);

    kuramoto.period();

    //kuramoto.write_peaks();
    //kuramoto.write_connect();
    kuramoto.write_results();

    gsl_odeiv_step_free(s);

    return 0;
}