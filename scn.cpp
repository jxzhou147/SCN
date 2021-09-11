#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <math.h>

#include <boost/numeric/odeint.hpp>

#include "parameter.h"

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double, 20> state_type;

/*
x[0] = M_P;
x[1] = M_C;
x[2] = M_B;
x[3] = P_C;
x[4] = C_C;
x[5] = P_CP;
x[6] = C_CP;
x[7] = PC_C;
x[8] = PC_N;
x[9] = PC_CP;
x[10] = PC_NP;
x[11] = B_C;
x[12] = B_CP;
x[13] = B_N;
x[14] = B_NP;
x[15] = I_N;
x[16] = Ca;
x[17] = Ca_store;
x[18] = VIP;
x[19] = CB_a;
*/

double g_K;
double g_L;
double g_Ca;
double E_Ca;
double g_KCa;
double g_ex;
double E_GABA;
double I_s;
double R_s;
double tau_m;

double a_v;
double b_v;
double c_v;
double P_K;
double V_reset;
double V_rest;
double R_m;
double theta;
double f_r;

double GABA;
double Cl_in;
double K_in;
double Na_in;

double beta;
double v_k;


void scn_firing(const state_type& x, state_type& dxdt, double t)
{
	g_K = g_K0 + v_gk * x[0] / (K_gk + x[0]);
	g_Ca = v_Ca * pow(x[0], n_Ca) / (K_Ca + pow(x[0], n_Ca));
	g_KCa = v_KCa * pow(x[4], n_KCa) / (K_KCa + pow(x[4], n_KCa));
	
	GABA = GABA_0 + v_GABA * x[18] / (K_GABA + x[18]);
	Cl_in = Cl_0 + v_Cl1 * x[0] / (K_Cl1 + x[0]) + v_Cl2 * pow(GABA, n_Cl) / (K_Cl2 + pow(GABA, n_Cl));

	E_Ca = k_B * T * log(Ca_ex / x[16]) / 2;
	E_GABA = k_B * T * log(Cl_ex / Cl_in);

	K_in = K_ex / exp(E_K / (k_B * T));
	Na_in = Na_ex / exp(E_Na / (k_B * T));
	P_K = v_PK * pow(x[11], npk) / (K_PK + pow(x[11], npk));
	a_v = 4 * P_Ca * x[16] + P_K * K_in + P_Na * Na_in + P_Cl * Cl_ex;
	b_v = P_K * K_in - P_K * K_ex + P_Na * Na_in - P_Na * Na_ex + P_Cl * Cl_ex - P_Cl * Cl_in;
	c_v = -1 * (P_K * K_ex + 4 * P_Ca * Ca_ex + P_Na * Na_ex + P_Cl * Cl_in);
	V_rest = k_B * T * log((-b_v + pow((pow(b_v, 2) - (4 * a_v) * c_v), 0.5)) / (2 * a_v));
	//////////The definition of g_ex is constrast with the profile of I_ex. The value and unit of K_ex1 is strange.
	g_ex = v_ex1 * pow(abs(g_Na * (V_rest - E_Na) * 10e-6), n_ex1) / (K_ex1 + pow(abs(g_Na * (V_rest - E_Na) * 10e-6), n_ex1)) + v_ex2 * pow(x[16], n_ex2) / (K_ex2 + pow(x[16], n_ex2));
	R_m = V_R * V_rest / (K__R + V_rest);
	g_L = 1 / R_m;
	V_reset = 4 + V_rest;
	I_s = g_Na * E_Na + g_Ca * E_Ca + g_K * E_K + g_KCa * E_K + g_L * E_L - g_ex * E_ex - g_GABA * E_GABA;
	R_s = 1 / (g_Na + g_Ca + g_K + g_KCa + g_L - g_ex - g_GABA);
	tau_m = C_m * R_s;
	theta = 20 + V_rest;
	f_r = -1 * pow((tau_m * log((theta - I_s * R_s) / (V_reset - I_s * R_s))), -1);


	beta = x[18] / (K_D + x[18]);
	v_k = V_MK * x[16] / (x[16] + K_MK) + V_beta * beta / (beta + K_beta);

	dxdt[0] = (v_sP + C_T * CB_T * x[19] / (K_C + CB_T * x[19])) * pow(x[13], n) / (pow(K_AP, n) + pow(x[13], n)) - v_mP * x[0] / (K_mP + x[0]) - k_dmp * x[0];
	dxdt[1] = v_sC * pow(x[13], n) / (pow(K_AC, n) + pow(x[13], n)) - v_mC * x[1] / (K_mC + x[1]) - k_dmc * x[1];
	dxdt[2] = v_sB * pow(K_IB, m) / (pow(K_IB, m) + pow(x[13], m)) - v_mB * x[2] / (K_mB + x[2]) - k_dmb * x[2];
	dxdt[3] = k_sP * x[0] - V_1P * x[3] / (K_p + x[3]) + V_2P * x[5] / (K_dp + x[5]) + k_4 * x[7] - k_3 * x[3] * x[4] - k_dn * x[3];
	dxdt[4] = k_sC * x[1] - V_1C * x[4] / (K_p + x[4]) + V_2C * x[6] / (K_dp + x[6]) + k_4 * x[7] - k_3 * x[3] * x[4] - k_dnc * x[4];
	dxdt[5] = V_1P * x[3] / (K_p + x[3]) - V_2P * x[5] / (K_dp + x[5]) - v_dPC * x[5] / (K_d + x[5]) - k_dn * x[5];
	dxdt[6] = V_1C * x[4] / (K_p + x[4]) - V_2C * x[6] / (K_dp + x[6]) - v_dCC * x[6] / (K_d + x[6]) - k_dn * x[6];
	dxdt[7] = -V_1PC * x[7] / (K_p + x[7]) + V_2PC * x[9] / (K_dp + x[9]) - k_4 * x[7] + k_3 * x[3] * x[4] + k_2 * x[8] - k_1 * x[7] - k_dn * x[7];
	dxdt[8] = -V_3PC * x[8] / (K_p + x[8]) + V_4PC * x[10] / (K_dp + x[10]) - k_2 * x[8] + k_1 * x[7] - k_7 * x[13] * x[8] + k_8 * x[15] - k_dn * x[8];
	dxdt[9] = V_1PC * x[7] / (K_p + x[7]) - V_2PC * x[9] / (K_dp + x[9]) - v_dPCC * x[9] / (K_d + x[9]) - k_dn * x[9];
	dxdt[10] = V_3PC * x[8] / (K_p + x[8]) - V_4PC * x[10] / (K_dp + x[10]) - v_dPCN * x[10] / (K_d + x[10]) - k_dn * x[10];
	dxdt[11] = k_sB * x[2] - V_1B * x[11] / (K_p + x[11]) + V_2B * x[12] / (K_dp + x[12]) - k_5 * x[11] + k_6 * x[13] - k_dn * x[11];
	dxdt[12] = V_1B * x[11] / (K_p + x[11]) - V_2B * x[12] / (K_dp + x[12]) - v_dBC * x[12] / (K_d + x[12]) - k_dn * x[12];
	dxdt[13] = -V_3B * x[13] / (K_p + x[13]) + V_4B * x[14] / (K_dp + x[14]) + k_5 * x[11] - k_6 * x[13] - k_7 * x[13] * x[8] + k_8 * x[15] - k_dn * x[13];
	dxdt[14] = V_3B * x[13] / (K_p + x[13]) - V_4B * x[14] / (K_dp + x[14]) - v_dBN * x[14] / (K_d + x[14]) - k_dn * x[14];
	dxdt[15] = -k_8 * x[15] + k_7 * x[13] * x[8] - v_dIN * x[15] / (K_d + x[15]) - k_dn * x[15];
	
	dxdt[16] = v_v0 * pow(x[11], nv0) / (K_v0 + pow(x[11], nv0)) + v_1 * beta_IP3 - (v_kk * pow(x[4], nkk) / (K_kk + pow(x[4], nkk))) * pow(x[16], v) - V_M2 * pow(x[16], n1) / (pow(K_2, n1) + pow(x[16], n1)) + V_M3 * (pow(x[17], m1) / (pow(K_R, m1) + pow(x[17], m1))) * pow(x[16], p1) / (pow(K_A, p1) + pow(x[16], p1)) + k_f * x[17];
	dxdt[17] = V_M2 * pow(x[16], n1) / (pow(K_2, n1) + pow(x[16], n1)) - V_M3 * (pow(x[17], m1) / (pow(K_R, m1) + pow(x[17], m1))) * pow(x[16], p1) / (pow(K_A, p1) + pow(x[16], p1)) - k_f * x[17];
	dxdt[18] = v_VIP * pow(f_r, n_VIP) / (K_VIP + pow(f_r, n_VIP)) - k_dVIP * pow(x[18], n_dVIP);
	dxdt[19] = (v_P / CB_T) * ((v_k / v_P) * (1 - x[19]) / (K_1 + 1 - x[19]) - x[19] / (K__2 + x[19]));
	
}


void write_SCN(const state_type& x, const double t)
{
	ofstream ofile;
	ofile.open("firing_model.csv", ios::app);
//	ofile << t << '\t' << x[0] << '\t' << g_K << '\t' << E_Ca << '\t' << g_Ca << '\t' << g_KCa << '\t' << E_GABA << '\t' << f_r << '\t' << P_K << '\t' << V_rest << '\t' << R_m << '\t' << x[19] << '\t' << x[11] << endl;
	ofile << t << '\t' << x[0] << '\t' << f_r << '\t' << g_ex << endl;
	ofile.close();
}

int main(int argc, char** argv)
{
	ofstream ofile;

	ofile.open("firing_model.csv");
	ofile.close();

	state_type x = { 2.80, 2.00, 7.94, 0.40, 12.0, 0.13, 9.00, 1.26, 0.16, 0.20, 0.091, 2.41, 0.48, 1.94, 0.32, 0.05, 0.10, 0.10, 0.00, 0.12 }; // initial conditions

	integrate(scn_firing, x, 0.0, 800.0, 0.1, write_SCN);
}
