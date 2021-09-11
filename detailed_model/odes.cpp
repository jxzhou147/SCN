#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>
#include <iostream>
#include <random>

#include "parameter.h"
#include "funs.h"

/*
for ith neuron (i = 0, 1, ..., N-1)
y[i * Ns] = M_P;
y[i * Ns + 1] = M_c;
y[i * Ns + 2] = M_B;
y[i * Ns + 3] = P_C;
y[i * Ns + 4] = C_C;
y[i * Ns + 5] = P_CP;
y[i * Ns + 6] = C_CP;
y[i * Ns + 7] = PC_C;
y[i * Ns + 8] = PC_N;
y[i * Ns + 9] = PC_CP;
y[i * Ns + 10] = PC_NP;
y[i * Ns + 11] = B_C;
y[i * Ns + 12] = B_CP;
y[i * Ns + 13] = B_N;
y[i * Ns + 14] = B_NP;
y[i * Ns + 15] = I_N;
y[i * Ns + 16] = Ca;
y[i * Ns + 17] = Ca_store;
y[i * Ns + 18] = VIP;
y[i * Ns + 19] = CB_a;
y[i * Ns + 20] = GABA_f
*/

int scn_firing(double t, const double y[], double f[], void *params)
{
    Params *param;
    param = (struct Params *)params;

    double fr_w[N];
    int light = 0;

    //double p_input = 0.5 * sin(2 * PI * (t - 12.0) / 24.0) + 0.5;

    for (int i = 0; i < N; i++)
    {
        //// Glutamate pathways
        double Glu = v_Glu * y[i * Ns] / (K_Glu + y[i * Ns]);
        double b_GluR = Glu / (K_GluR + Glu);
        /*
        //// Intercellular pathways
        //// For VIP
        int k_connect = 0;
        double VIP_connect = 0;
        double gamma = 0;
        for (size_t j = 0; j < N; j++)
        {
            if (param->a_vip[j][i] == 1)
            {
                k_connect += 1;
                VIP_connect += y[j * Ns + 18];
            }
        }
        if (k_connect != 0)
        {
            gamma = VIP_connect / k_connect;
        }
        double beta = gamma / (K_D + gamma);

        /// For GABA
        int mu_connect = 0;
        double GABA_connect = 0;
        double delta = 0;
        for (size_t j = 0; j < N; j++)
        {
            if (param->a_gaba[j][i] == 1)
            {
                mu_connect += 1;
                GABA_connect += GABA_0 + y[j * Ns + 20];
            }
        }
        if (mu_connect != 0)
            delta = GABA_connect / mu_connect;
*/

        // 0/0.8/0.2
        //// Intercellular pathways
        //// For VIP
        int k_connect_8 = 0;
        double VIP_connect_8 = 0;
        int k_connect_2 = 0;
        double VIP_connect_2 = 0;
        double gamma = 0;
        for (size_t j = 0; j < N; j++)
        {
            if (fabs(param->a_vip[j][i] - 0.8) < 1e-6)
            {
                k_connect_8 += 1;
                VIP_connect_8 += y[j * Ns + 18];
            }
            else if (fabs(param->a_vip[j][i] - 0.2) < 1e-6)
            {
                k_connect_2 += 1;
                VIP_connect_2 += y[j * Ns + 18];
            }
        }
        if (k_connect_8 != 0)
        {
            gamma += 0.8 * VIP_connect_8 / k_connect_8;
        }
        if (k_connect_2 != 0)
        {
            gamma += 0.2 * VIP_connect_2 / k_connect_2;
        }
        double beta = gamma / (K_D + gamma);

        /// For GABA
        int mu_connect_8 = 0;
        double GABA_connect_8 = 0;
        int mu_connect_2 = 0;
        double GABA_connect_2 = 0;
        double delta = 0;
        for (size_t j = 0; j < N; j++)
        {
            if (fabs(param->a_gaba[j][i] - 0.8) < 1e-6)
            {
                mu_connect_8 += 1;
                GABA_connect_8 += GABA_0 + y[j * Ns + 20];
            }
            else if (fabs(param->a_gaba[j][i] - 0.2) < 1e-6)
            {
                mu_connect_2 += 1;
                GABA_connect_2 += GABA_0 + y[j * Ns + 20];
            }
        }
        if (mu_connect_8 != 0)
            delta += 0.8 * GABA_connect_8 / mu_connect_8;
        if (mu_connect_2 != 0)
            delta += 0.2 * GABA_connect_2 / mu_connect_2;
        /*
        // uniform
        double gamma = 0;
        double delta = 0;
        double interaction_each_vip[10] = {0};
        double interaction_each_gaba[10] = {0};
        int e_each[10] = {0};
        int _which_;

        for (size_t j = 0; j < N; j++)
        {
            if (param->a_vip[j][i] > 1e-3)
            {
                _which_ = ((int)(param->a_vip[j][i] * 10)) % 10;
                interaction_each_vip[_which_] += param->a_vip[j][i] * y[j * Ns + 18] / 5.0;             // 0.05 + 0.15 + ... + 0.95 = 5.0
                interaction_each_gaba[_which_] += param->a_vip[j][i] * (GABA_0 + y[j * Ns + 20]) / 5.0; // 0.05 + 0.15 + ... + 0.95 = 5.0
                e_each[_which_]++;
            }
        }

        for (int j = 0; j < 10; j++)
        {
            if (e_each[j] != 0)
            {
                gamma += interaction_each_vip[j] / e_each[j];
                delta += interaction_each_gaba[j] / e_each[j];
            }
        }

        double beta = gamma / (K_D + gamma);
*/
        ///// The effect of light entrainment
        //double pulse = 500;
        double delay = 0;
        //if (t > 396)
        //  delay = 6;

        //if ((t > 150) & (fmod(t + delay, 24.0) >= (12.0)) & (i < N_v))
        //if ((t > 150) & (fmod(t - delay, 24.0) >= (24 - param->day)) & (i < N_v))
        //if ((t > 480) & (fmod(t - delay, 24.0) <= (24 - param->day)) & (i < N_v))

        //if ((t >= param->pulse) && (t < (param->pulse + delay)))

        //if ((t > 150) && (fmod(t, param->daylength) >= (param->daylength / 2.0)) && (i < (N_v * (0.5 * sin(2 * PI * (t - 10 - param->daylength / 2.0) / param->daylength) + 0.5))))
        //if ((t > 150) & (fmod(t, param->daylength) >= (param->daylength - param->day)) & (i < N_v * p_input))
        //if ((t > 150) && (fmod(t, 24.0) >= (24 - param->day)) && (i < N_v) && ((rand() / double(RAND_MAX)) < (0.5 * sin(2 * PI * (t - 12.0) / 24.0) + 0.5)))
        /*
        if ((t > 150) && (fmod(t - delay, param->daylength) >= (param->daylength - param->day)) && (i < (N_v)))
        {
            light = 1;
            b_GluR = 1;
            beta = param->beta_light;
        }
*/
        ///// v to d weight cutting 50% in days
        //if ((t > 150) && (fmod(t, param->daylength) >= (param->daylength - param->day)) && (i >= (N_v)))
        //  beta = beta * 0.95;
        /*
        if (t > 600)
        {
            light = 1;
            b_GluR = 1;
            beta = param->beta_light;
        }*/
        /*
        ///// cut all the links in ventral SCN at night
        if ((i < N_v) && (fmod(t - delay, 24.0) < (12.0)) && ((rand() / double(RAND_MAX)) < 0.2))
        {
            beta = 0;
            delta = GABA_0;
        }*/

        if (t < 150)
        {
            //gamma = 0;
            beta = 0;
            delta = 0;
        }

        double f_r = firingRate(y[i * Ns], y[i * Ns + 4], y[i * Ns + 11], y[i * Ns + 16], delta, t);

        fr_w[i] = f_r;

        double v_k = V_MK * y[i * Ns + 16] / (y[i * Ns + 16] + K_MK) + V_beta * beta / (beta + K_beta);
        if (t < 150)
            v_k = 0;

        //if ((t > 334) & (i == 400))
        //  std::cout << t << '\t' << delta << '\t' << beta << '\t' << f_r << '\t' << v_k << '\t' << y[i * Ns] << '\t' << y[i * Ns + 4] << '\t' << y[i * Ns + 11] << '\t' << y[i * Ns + 16] << std::endl;

        f[i * Ns] = (param->v_sP[i] + C_T * CB_T * y[i * Ns + 19] / (K_C + CB_T * y[i * Ns + 19])) * pow(y[i * Ns + 13], n) / (pow(K_AP, n) + pow(y[i * Ns + 13], n)) - v_mP * y[i * Ns] / (K_mP + y[i * Ns]) - k_dmp * y[i * Ns];
        f[i * Ns + 1] = v_sC * pow(y[i * Ns + 13], n) / (pow(K_AC, n) + pow(y[i * Ns + 13], n)) - v_mC * y[i * Ns + 1] / (K_mC + y[i * Ns + 1]) - k_dmc * y[i * Ns + 1];
        f[i * Ns + 2] = param->v_sB[i] * pow(K_IB, m) / (pow(K_IB, m) + pow(y[i * Ns + 13], m)) - param->v_mB[i] * y[i * Ns + 2] / (K_mB + y[i * Ns + 2]) - k_dmb * y[i * Ns + 2];
        f[i * Ns + 3] = k_sP * y[i * Ns] - V_1P * y[i * Ns + 3] / (K_p + y[i * Ns + 3]) + V_2P * y[i * Ns + 5] / (K_dp + y[i * Ns + 5]) + k_4 * y[i * Ns + 7] - k_3 * y[i * Ns + 3] * y[i * Ns + 4] - k_dn * y[i * Ns + 3];
        f[i * Ns + 4] = k_sC * y[i * Ns + 1] - V_1C * y[i * Ns + 4] / (K_p + y[i * Ns + 4]) + V_2C * y[i * Ns + 6] / (K_dp + y[i * Ns + 6]) + k_4 * y[i * Ns + 7] - k_3 * y[i * Ns + 3] * y[i * Ns + 4] - k_dnc * y[i * Ns + 4];
        f[i * Ns + 5] = V_1P * y[i * Ns + 3] / (K_p + y[i * Ns + 3]) - V_2P * y[i * Ns + 5] / (K_dp + y[i * Ns + 5]) - v_dPC * y[i * Ns + 5] / (K_d + y[i * Ns + 5]) - k_dn * y[i * Ns + 5];
        f[i * Ns + 6] = V_1C * y[i * Ns + 4] / (K_p + y[i * Ns + 4]) - V_2C * y[i * Ns + 6] / (K_dp + y[i * Ns + 6]) - v_dCC * y[i * Ns + 6] / (K_d + y[i * Ns + 6]) - k_dn * y[i * Ns + 6];
        f[i * Ns + 7] = -V_1PC * y[i * Ns + 7] / (K_p + y[i * Ns + 7]) + V_2PC * y[i * Ns + 9] / (K_dp + y[i * Ns + 9]) - k_4 * y[i * Ns + 7] + k_3 * y[i * Ns + 3] * y[i * Ns + 4] + k_2 * y[i * Ns + 8] - k_1 * y[i * Ns + 7] - k_dn * y[i * Ns + 7];
        f[i * Ns + 8] = -V_3PC * y[i * Ns + 8] / (K_p + y[i * Ns + 8]) + V_4PC * y[i * Ns + 10] / (K_dp + y[i * Ns + 10]) - k_2 * y[i * Ns + 8] + k_1 * y[i * Ns + 7] - k_7 * y[i * Ns + 13] * y[i * Ns + 8] + k_8 * y[i * Ns + 15] - k_dn * y[i * Ns + 8];
        f[i * Ns + 9] = V_1PC * y[i * Ns + 7] / (K_p + y[i * Ns + 7]) - V_2PC * y[i * Ns + 9] / (K_dp + y[i * Ns + 9]) - v_dPCC * y[i * Ns + 9] / (K_d + y[i * Ns + 9]) - k_dn * y[i * Ns + 9];
        f[i * Ns + 10] = V_3PC * y[i * Ns + 8] / (K_p + y[i * Ns + 8]) - V_4PC * y[i * Ns + 10] / (K_dp + y[i * Ns + 10]) - v_dPCN * y[i * Ns + 10] / (K_d + y[i * Ns + 10]) - k_dn * y[i * Ns + 10];
        f[i * Ns + 11] = k_sB * y[i * Ns + 2] - V_1B * y[i * Ns + 11] / (K_p + y[i * Ns + 11]) + V_2B * y[i * Ns + 12] / (K_dp + y[i * Ns + 12]) - k_5 * y[i * Ns + 11] + k_6 * y[i * Ns + 13] - k_dn * y[i * Ns + 11];
        f[i * Ns + 12] = V_1B * y[i * Ns + 11] / (K_p + y[i * Ns + 11]) - V_2B * y[i * Ns + 12] / (K_dp + y[i * Ns + 12]) - v_dBC * y[i * Ns + 12] / (K_d + y[i * Ns + 12]) - k_dn * y[i * Ns + 12];
        f[i * Ns + 13] = -V_3B * y[i * Ns + 13] / (K_p + y[i * Ns + 13]) + V_4B * y[i * Ns + 14] / (K_dp + y[i * Ns + 14]) + k_5 * y[i * Ns + 11] - k_6 * y[i * Ns + 13] - k_7 * y[i * Ns + 13] * y[i * Ns + 8] + k_8 * y[i * Ns + 15] - k_dn * y[i * Ns + 13];
        f[i * Ns + 14] = V_3B * y[i * Ns + 13] / (K_p + y[i * Ns + 13]) - V_4B * y[i * Ns + 14] / (K_dp + y[i * Ns + 14]) - v_dBN * y[i * Ns + 14] / (K_d + y[i * Ns + 14]) - k_dn * y[i * Ns + 14];
        f[i * Ns + 15] = -k_8 * y[i * Ns + 15] + k_7 * y[i * Ns + 13] * y[i * Ns + 8] - v_dIN * y[i * Ns + 15] / (K_d + y[i * Ns + 15]) - k_dn * y[i * Ns + 15];

        f[i * Ns + 16] = v_v0 * real_root_no(y[i * Ns + 11], nv0) / (K_v0 + real_root_no(y[i * Ns + 11], nv0)) + v_GluR * b_GluR + v_1 * beta_IP3 - (v_kk * real_root_no(y[i * Ns + 4], nkk) / (K_kk + real_root_no(y[i * Ns + 4], nkk))) * pow(y[i * Ns + 16], v) - V_M2 * real_root_transfer(y[i * Ns + 16], n1) / (real_root_transfer(K_2, n1) + real_root_transfer(y[i * Ns + 16], n1)) + V_M3 * (pow(y[i * Ns + 17], m1) / (pow(K_R, m1) + pow(y[i * Ns + 17], m1))) * real_root_transfer(y[i * Ns + 16], p1) / (real_root_transfer(K_A, p1) + real_root_transfer(y[i * Ns + 16], p1)) + k_f * y[i * Ns + 17];
        f[i * Ns + 17] = V_M2 * real_root_transfer(y[i * Ns + 16], n1) / (real_root_transfer(K_2, n1) + real_root_transfer(y[i * Ns + 16], n1)) - V_M3 * (pow(y[i * Ns + 17], m1) / (pow(K_R, m1) + pow(y[i * Ns + 17], m1))) * real_root_transfer(y[i * Ns + 16], p1) / (real_root_transfer(K_A, p1) + real_root_transfer(y[i * Ns + 16], p1)) - k_f * y[i * Ns + 17];
        f[i * Ns + 18] = v_VIP * pow(f_r, n_VIP) / (K_VIP + pow(f_r, n_VIP)) - k_dVIP * real_root_transfer(y[i * Ns + 18], n_dVIP);
        f[i * Ns + 19] = (v_P / CB_T) * ((v_k / v_P) * (1 - y[i * Ns + 19]) / (K_1 + 1 - y[i * Ns + 19]) - y[i * Ns + 19] / (K__2 + y[i * Ns + 19]));
        f[i * Ns + 20] = v_GABA * pow(f_r, n_GABA) / (K_GABA + pow(f_r, n_GABA)) - k_dGABA * real_root_transfer(y[i * Ns + 20], n_dGABA);
        //f[i * Ns + 20] = v_GABA * pow(f_r, n_GABA) / (K_GABA + pow(f_r, n_GABA)) - k_dGABA * pow(GABA_0 + y[i * Ns + 20], n_dGABA);
        /*
        if (y[i * Ns + 18] > 0)
            f[i * Ns + 18] = v_VIP * pow(f_r, n_VIP) / (K_VIP + pow(f_r, n_VIP)) - k_dVIP * pow(y[i * Ns + 18], n_dVIP);
        else
        {
            f[i * Ns + 18] = v_VIP * pow(f_r, n_VIP) / (K_VIP + pow(f_r, n_VIP)) + k_dVIP * pow(-1.0 * y[i * Ns + 18], n_dVIP);
        }

        if (y[i * Ns + 20] > 0)
        {
            f[i * Ns + 20] = v_GABA * pow(f_r, n_GABA) / (K_GABA + pow(f_r, n_GABA)) - k_dGABA * pow(y[i * Ns + 20], n_dGABA);
        }
        else
        {
            f[i * Ns + 20] = v_GABA * pow(f_r, n_GABA) / (K_GABA + pow(f_r, n_GABA)) + k_dGABA * pow(-1.0 * y[i * Ns + 20], n_dGABA);
            //if (t > 335.1)
            std::cout << y[i * Ns + 20] << '\t' << f[i * Ns + 20] << '\t';
        }
*/
        //if (t > 335.1)
        //  std::cout << f[i * Ns + 20] << '\t';
    }

    //if (t > 300)
    //  write_fr(fr_w, t, light);
    //write_input(p_input, t, light);

    return GSL_SUCCESS;
}