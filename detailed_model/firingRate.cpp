#include <iostream>
#include <math.h>
#include "parameter.h"
#include "funs.h"

double firingRate(double M_P, double C_C, double B_C, double Ca, double delta, double t)
{
    if (M_P < 0)
        M_P = 0;
    double g_K = g_K0 + v_gk * M_P / (K_gk + M_P);
    double g_Ca = v_Ca * real_root_transfer(M_P, n_Ca) / (K_Ca + real_root_transfer(M_P, n_Ca));
    double g_KCa = v_KCa * real_root_zero(C_C, n_KCa) / (K_KCa + real_root_zero(C_C, n_KCa));

    if ((delta <= 0) | (isnan(delta)))
        delta = GABA_0;

    double Cl_in = Cl_0 + v_Cl1 * M_P / (K_Cl1 + M_P) + v_Cl2 * real_root_transfer(delta, n_Cl) / (K_Cl2 + real_root_transfer(delta, n_Cl));

    double E_Na = E_Na0 * T / T_room;
    double E_K = E_K0 * T / T_room;
    double E_L = E_L0 * T / T_room;
    double E_Ca = k_B * T * log(Ca_ex / Ca) / 2;
    double E_GABA = -1 * k_B * T * log(Cl_ex / Cl_in);
    //double E_GABA = k_B * T * log(Cl_ex / Cl_in);   my old code

    double K_in = K_ex / exp(E_K / (k_B * T));
    double Na_in = Na_ex / exp(E_Na / (k_B * T));
    double P_K = v_PK * real_root_zero(B_C, npk) / (K_PK + real_root_zero(B_C, npk));
    //double a_v = 4 * P_Ca * Ca + P_K * K_in + P_Na * Na_in + P_Cl * Cl_ex;
    double a_v = 4 * P_Ca * Ca * pow(10, -3) + P_K * K_in + P_Na * Na_in + P_Cl * Cl_ex; // umol to mmol
    double b_v = P_K * K_in - P_K * K_ex + P_Na * Na_in - P_Na * Na_ex + P_Cl * Cl_ex - P_Cl * Cl_in;
    //double c_v = -1 * (P_K * K_ex + 4 * P_Ca * Ca_ex + P_Na * Na_ex + P_Cl * Cl_in);
    double c_v = -1 * (P_K * K_ex + 4 * P_Ca * Ca_ex * pow(10, -3) + P_Na * Na_ex + P_Cl * Cl_in);
    double inlog = (-b_v + real_root_no((pow(b_v, 2) - (4 * a_v) * c_v), 0.5)) / (2 * a_v);
    double V_rest = k_B1 * T * 2.303 * log10(inlog); // ln or lg ? KT?
    //double V_rest = k_B1 * T * 2.303 * log10((-b_v + pow((pow(b_v, 2) - (4 * a_v) * c_v), 0.5)) / (2 * a_v)); // ln or lg ? KT?
    //double g_ex = v_ex1 * pow(abs(g_Na * (V_rest - E_Na) * 10e-6), n_ex1) / (K_ex1 + pow(abs(g_Na * (V_rest - E_Na) * 10e-6), n_ex1)) + v_ex2 * pow(Ca, n_ex2) / (K_ex2 + pow(Ca, n_ex2));
    double g_ex = v_ex1 * real_root_no(abs(g_Na * (V_rest - E_Na)), n_ex1) / (K_ex1 + real_root_no(abs(g_Na * (V_rest - E_Na)), n_ex1)) + v_ex2 * real_root_zero(Ca, n_ex2) / (K_ex2 + real_root_zero(Ca, n_ex2));
    double R_m = V_R * V_rest / (K__R + V_rest);
    double g_L = 1 / R_m;
    double V_reset = 4 + V_rest;
    double I_s = g_Na * E_Na + g_Ca * E_Ca + g_K * E_K + g_KCa * E_K + g_L * E_L - g_ex * E_ex - g_GABA * E_GABA;
    double R_s = 1 / (g_Na + g_Ca + g_K + g_KCa + g_L - g_ex - g_GABA);
    double tau_m = C_m * R_s;
    double theta = 20 + V_rest;
    double phi = (theta - I_s * R_s) / (V_reset - I_s * R_s);

    double f_r;
    if (phi > 0)
        f_r = -1 * real_root_zero((tau_m * log(phi)), -1);
    else
        f_r = 3;

    //double f_r = -1 * pow((tau_m * log((theta - I_s * R_s) / (V_reset - I_s * R_s))), -1);
    if (isnan(f_r))
    {
        std::cout << Ca << '/' << M_P << '/' << delta << '/' << Cl_in << '/' << a_v << '/' << b_v << '/' << '/' << c_v << '/' << inlog << '/' << tau_m << '/' << V_rest << '/' << '/' << M_P << '/' << g_Ca << '/' << g_K << '/' << g_KCa << '/' << g_L << '/' << g_ex << '/' << I_s << '/' << R_s << '/' << phi << '/' << f_r << '\t';
        //f_r = 3;
    }
    return f_r;
}