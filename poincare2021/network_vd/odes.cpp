#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>
#include <iostream>
#include <random>

#include "parameter.h"
#include "funs.h"

//////// 2N-dim, 0:N-1 -> x, N:2N-1 -> y
int poincare(double t, const double y[], double f[], void *params)
{
    Params *param;
    param = (struct Params *)params;

    double interaction[N];
    int light = 0;

    if (fmod(t, 24) < (param->daylength))
        light = 1;
    else
        light = 0;

    for (size_t i = 0; i < N; i++)
    {
        interaction[i] = 0;
        if (param->ed[i] != 0)
        {
            for (size_t j = 0; j < N; j++)
            {
                //double flu = (rand() % (N_rand + 1) / (float)(N_rand + 1)) + 0.5; //random strength of links
                double flu = 1; //constant strength of links
                //double flu = sin(Omega * t) + 1;
                /*
                if (i < N_v)
                {
                    if (light == 1)
                        flu = 1;
                    else
                    {
                        //if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < 0.2)
                        flu = 0;
                    }
                }
*/
                //if (param->light == 1)
                //  flu = 1.5;
                interaction[i] += param->a[j][i] * flu * y[j] / param->ed[i];
                //param->interaction[i] += param->a[j][i] * flu * y[j];
            }
        }
    }

    for (size_t i = 0; i < N_v; i++)
    {
        //f[i] = omega_v * param->eta[i] + K * param->interaction[i] + K_f * param->light * phase(param->theta[i]);
        //f[i] = omega_v * param->eta[i] + K * param->interaction[i] + K_f * sin(Omega * t - y[i]);
        //f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_v * param->eta[i] + param->K * param->interaction[i] + param->K_f * param->light * phase_light(param->theta[i], param->daylength);
        f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_v * param->eta[i] + param->K * interaction[i] + param->K_f * light;
        //f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_v * param->eta[i];
        f[i + N] = gam * y[i + N] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) + y[i] * omega_v * param->eta[i];
        //f[i + N] = gam * y[i + N] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) + y[i] * omega_v * param->eta[i];
    }

    for (size_t i = N_v; i < N; i++)
    {
        f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_d * param->eta[i] + param->K * interaction[i];

        //f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_d * param->eta[i];
        f[i + N] = gam * y[i + N] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) + y[i] * omega_d * param->eta[i];
    }

    return GSL_SUCCESS;
}