#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>
#include <iostream>
#include <random>

#include "parameter.h"
#include "funs.h"

/*
y[0] = M_P;
y[1] = M_c;
y[2] = M_B;
y[3] = P_C;
y[4] = C_C;
y[5] = P_CP;
y[6] = C_CP;
y[7] = PC_C;
y[8] = PC_N;
y[9] = PC_CP;
y[10] = PC_NP;
y[11] = B_C;
y[12] = B_CP;
y[13] = B_N;
y[14] = B_NP;
y[15] = I_N;
y[16] = Ca;
y[17] = Ca_store;
y[18] = VIP;
y[19] = CB_a;
y[20] = GABA_f
*/

int scn_firing(double t, const double y[], double f[], void *params)
{
    int light = 0;
    double daylength = *(double *)params;
    double day = 18;
    double beta_light = 0.8;

    //// Glutamate pathways
    double Glu = v_Glu * y[0] / (K_Glu + y[0]);
    double b_GluR = Glu / (K_GluR + Glu);

    double gamma = y[18];
    double beta = gamma / (K_D + gamma);

    double delta = GABA_0 + y[20];

    ///// The effect of light entrainment
    if (fmod(t, daylength) >= (daylength - day))
    {
        light = 1;
        b_GluR = 1;
        beta = beta_light;
    }

    double f_r = firingRate(y[0], y[4], y[11], y[16], delta, t);

    double v_k = V_MK * y[16] / (y[16] + K_MK) + V_beta * beta / (beta + K_beta);

    f[0] = (v_sP0 + C_T * CB_T * y[19] / (K_C + CB_T * y[19])) * pow(y[13], n) / (pow(K_AP, n) + pow(y[13], n)) - v_mP * y[0] / (K_mP + y[0]) - k_dmp * y[0];
    f[1] = v_sC * pow(y[13], n) / (pow(K_AC, n) + pow(y[13], n)) - v_mC * y[1] / (K_mC + y[1]) - k_dmc * y[1];
    f[2] = v_sB0 * pow(K_IB, m) / (pow(K_IB, m) + pow(y[13], m)) - v_mB0 * y[2] / (K_mB + y[2]) - k_dmb * y[2];
    f[3] = k_sP * y[0] - V_1P * y[3] / (K_p + y[3]) + V_2P * y[5] / (K_dp + y[5]) + k_4 * y[7] - k_3 * y[3] * y[4] - k_dn * y[3];
    f[4] = k_sC * y[1] - V_1C * y[4] / (K_p + y[4]) + V_2C * y[6] / (K_dp + y[6]) + k_4 * y[7] - k_3 * y[3] * y[4] - k_dnc * y[4];
    f[5] = V_1P * y[3] / (K_p + y[3]) - V_2P * y[5] / (K_dp + y[5]) - v_dPC * y[5] / (K_d + y[5]) - k_dn * y[5];
    f[6] = V_1C * y[4] / (K_p + y[4]) - V_2C * y[6] / (K_dp + y[6]) - v_dCC * y[6] / (K_d + y[6]) - k_dn * y[6];
    f[7] = -V_1PC * y[7] / (K_p + y[7]) + V_2PC * y[9] / (K_dp + y[9]) - k_4 * y[7] + k_3 * y[3] * y[4] + k_2 * y[8] - k_1 * y[7] - k_dn * y[7];
    f[8] = -V_3PC * y[8] / (K_p + y[8]) + V_4PC * y[10] / (K_dp + y[10]) - k_2 * y[8] + k_1 * y[7] - k_7 * y[13] * y[8] + k_8 * y[15] - k_dn * y[8];
    f[9] = V_1PC * y[7] / (K_p + y[7]) - V_2PC * y[9] / (K_dp + y[9]) - v_dPCC * y[9] / (K_d + y[9]) - k_dn * y[9];
    f[10] = V_3PC * y[8] / (K_p + y[8]) - V_4PC * y[10] / (K_dp + y[10]) - v_dPCN * y[10] / (K_d + y[10]) - k_dn * y[10];
    f[11] = k_sB * y[2] - V_1B * y[11] / (K_p + y[11]) + V_2B * y[12] / (K_dp + y[12]) - k_5 * y[11] + k_6 * y[13] - k_dn * y[11];
    f[12] = V_1B * y[11] / (K_p + y[11]) - V_2B * y[12] / (K_dp + y[12]) - v_dBC * y[12] / (K_d + y[12]) - k_dn * y[12];
    f[13] = -V_3B * y[13] / (K_p + y[13]) + V_4B * y[14] / (K_dp + y[14]) + k_5 * y[11] - k_6 * y[13] - k_7 * y[13] * y[8] + k_8 * y[15] - k_dn * y[13];
    f[14] = V_3B * y[13] / (K_p + y[13]) - V_4B * y[14] / (K_dp + y[14]) - v_dBN * y[14] / (K_d + y[14]) - k_dn * y[14];
    f[15] = -k_8 * y[15] + k_7 * y[13] * y[8] - v_dIN * y[15] / (K_d + y[15]) - k_dn * y[15];

    f[16] = v_v0 * real_root_no(y[11], nv0) / (K_v0 + real_root_no(y[11], nv0)) + v_GluR * b_GluR + v_1 * beta_IP3 - (v_kk * real_root_no(y[4], nkk) / (K_kk + real_root_no(y[4], nkk))) * pow(y[16], v) - V_M2 * real_root_transfer(y[16], n1) / (real_root_transfer(K_2, n1) + real_root_transfer(y[16], n1)) + V_M3 * (pow(y[17], m1) / (pow(K_R, m1) + pow(y[17], m1))) * real_root_transfer(y[16], p1) / (real_root_transfer(K_A, p1) + real_root_transfer(y[16], p1)) + k_f * y[17];
    f[17] = V_M2 * real_root_transfer(y[16], n1) / (real_root_transfer(K_2, n1) + real_root_transfer(y[16], n1)) - V_M3 * (pow(y[17], m1) / (pow(K_R, m1) + pow(y[17], m1))) * real_root_transfer(y[16], p1) / (real_root_transfer(K_A, p1) + real_root_transfer(y[16], p1)) - k_f * y[17];
    f[18] = v_VIP * pow(f_r, n_VIP) / (K_VIP + pow(f_r, n_VIP)) - k_dVIP * real_root_transfer(y[18], n_dVIP);
    f[19] = (v_P / CB_T) * ((v_k / v_P) * (1 - y[19]) / (K_1 + 1 - y[19]) - y[19] / (K__2 + y[19]));
    f[20] = v_GABA * pow(f_r, n_GABA) / (K_GABA + pow(f_r, n_GABA)) - k_dGABA * real_root_transfer(y[20], n_dGABA);

    //if (t > 300)
    //  write_fr(fr_w, t, light);
    //write_input(p_input, t, light);

    return GSL_SUCCESS;
}