#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>
#include <iostream>
#include <random>

#include "parameter.h"
#include "new_funs.h"

//////// 2N-dim, 0:N-1 -> x, N:2N-1 -> y
int poincare(double t, const double y[], double f[], void *params)
{
    Params *param;
    param = (struct Params *)params;

    double interaction[N];
    int light = 0;

    if (fmod(t, param->tcycle) < (param->daylength))
        light = 1;
    else
        light = 0;

    // no average
    for (size_t i = 0; i < N; i++)
    {
        interaction[i] = 0;
        for (size_t j = 0; j < N; j++)
        {
            interaction[i] += param->a[j][i] * y[j];
        }
    }
    /*
    // 0/1
    for (size_t i = 0; i < N; i++)
    {
        interaction[i] = 0;
        if (param->ed[i] != 0)
        {
            for (size_t j = 0; j < N; j++)
            {
                interaction[i] += param->a[j][i] * y[j] / param->ed[i];
            }
        }
    }
*/
    /*
    // 0/0.8/0.2
    for (size_t i = 0; i < N; i++)
    {
        interaction[i] = 0;
        double interaction_8 = 0;
        double interaction_2 = 0;
        int e8 = 0;
        int e2 = 0;

        for (size_t j = 0; j < N; j++)
        {
            if (abs(param->a[j][i] - 0.8) < 1e-6)
            {
                interaction_8 += y[j];
                e8++;
            }
            else if (abs(param->a[j][i] - 0.2) < 1e-6)
            {
                interaction_2 += y[j];
                e2++;
            }
        }

        if (e8 != 0)
        {
            interaction[i] += 0.8 * interaction_8 / e8;
        }
        if (e2 != 0)
        {
            interaction[i] += 0.2 * interaction_2 / e2;
        }
    }
*/

    // uniform
    /*
    for (size_t i = 0; i < N; i++)
    {
        interaction[i] = 0;
        double interaction_each[10] = {0};
        int e_each[10] = {0};
        int _which_;

        for (size_t j = 0; j < N; j++)
        {
            if (param->a[j][i] > 1e-3)
            {
                _which_ = ((int)(param->a[j][i] * 10)) % 10;
                interaction_each[_which_] += param->a[j][i] * y[j] / 5.0; // 0.05 + 0.15 + ... + 0.95 = 5.0
                e_each[_which_]++;
            }
        }

        for (int j = 0; j < 10; j++)
        {
            if (e_each[j] != 0)
            {
                interaction[i] += interaction_each[j] / e_each[j];
            }
        }
    }
*/
    /*
// normal
    for (size_t i = 0; i < N; i++)
    {
        interaction[i] = 0;
        double interaction_sum = 0;
        int e = 0;

        for (size_t j = 0; j < N; j++)
        {
            if (fabs(param->a[j][i] - 0.0) > 1e-3)
            {
                interaction_sum += param->a[j][i] * y[j];
                e++;
            }
        }

        if (e != 0)
        {
            interaction[i] += interaction_sum / e;
        }
    }
*/
    for (size_t i = 0; i < N_v; i++)
    {
        f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_v * param->eta[i] + param->K * interaction[i] + param->K_f * light;
        f[i + N] = gam * y[i + N] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) + y[i] * omega_v * param->eta[i];
    }

    for (size_t i = N_v; i < N; i++)
    {
        f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_d * param->eta[i] + param->K * interaction[i];
        f[i + N] = gam * y[i + N] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) + y[i] * omega_d * param->eta[i];
    }

    return GSL_SUCCESS;
}