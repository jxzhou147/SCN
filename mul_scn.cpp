#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <math.h>

#include <boost/numeric/odeint.hpp>

#include "mul_parameter.h"

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::matrix< double > state_type;
typedef boost::numeric::ublas::matrix< double > connect_type;
typedef boost::array< double, N> neuron_type;

#define N_rand  99 //random number

/*
x(i, 0) = M_P;
x(i, 1) = M_C;
x(i, 2) = M_B;
x(i, 3) = P_C;
x(i, 4) = C_C;
x(i, 5) = P_CP;
x(i, 6) = C_CP;
x(i, 7) = PC_C;
x(i, 8) = PC_N;
x(i, 9) = PC_CP;
x(i, 10) = PC_NP;
x(i, 11) = B_C;
x(i, 12) = B_CP;
x(i, 13) = B_N;
x(i, 14) = B_NP;
x(i, 15) = I_N;
x(i, 16) = Ca;
x(i, 17) = Ca_store;
x(i, 18) = VIP;
x(i, 19) = CB_a;
x(i, 20) = GABA_f
*/

connect_type a(N, N);

neuron_type g_K;
neuron_type g_L;
neuron_type g_Ca;
neuron_type E_Ca;
neuron_type g_KCa;
neuron_type g_ex;
neuron_type E_GABA;
neuron_type I_s;
neuron_type R_s;
neuron_type tau_m;

neuron_type a_v;
neuron_type b_v;
neuron_type c_v;
neuron_type P_K;
neuron_type V_reset;
neuron_type V_rest;
neuron_type R_m;
neuron_type theta;
neuron_type f_r;

neuron_type GABA;
neuron_type Cl_in;
neuron_type K_in;
neuron_type Na_in;

neuron_type beta;
neuron_type v_k;

neuron_type gamma;			//the VIP concentration sensed by a neuron
neuron_type delta;			//the GABA concentartion sensed by a neuron

pair<int, int> index_to_coordinate(int i)
{
	pair<int, int> coor(i / N0, i % N0);
	return coor;
}

int coordinate_to_index(int m, int n)
{
	return m * N0 + n;
}

pair<int, int> boundary_convert(int m)				//return (m - 1, m + 1)
{
	pair<int, int> m_n;
	if (m == 0)
		m_n = make_pair(N0 - 1, m + 1);
	else if (m == N0 - 1)
		m_n = make_pair(m - 1, 0);
	else
		m_n = make_pair(m - 1, m + 1);
	return m_n;
}

void construct_connection(double p_link, connect_type& a)
{
	int m, n;	//the ith neuron's coordinate (m, n)  i = m*N + n
	for (size_t i = 0; i < N; i++)
	{
		m = index_to_coordinate(i).first;
		n = index_to_coordinate(i).second;

		//////////////generate near connections, using periodic boundary conditions
		a(i, coordinate_to_index(m, boundary_convert(n).first)) = 1;
		a(i, coordinate_to_index(m, boundary_convert(n).second)) = 1;
		a(i, coordinate_to_index(boundary_convert(m).first, n)) = 1;
		a(i, coordinate_to_index(boundary_convert(m).second, n)) = 1;

		/////////////////////generate long range connections
		for (size_t j = 0; j < N; j++)
		{
			if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link)
				a(i, j) = 1;
		}
	}
}

void scn_firing(const state_type& x, state_type& dxdt, double t)
{
	for (size_t i = 0; i < N; i++)
	{
		int k_connect = 0;
		double VIP_connect = 0;
//		for (size_t j = 0; j < (N * p_vip); j++)
		for (size_t j = 0; j < N; j++)
		{
			if ((a(j, i) == 1)) // & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_vip)
			{
				k_connect += 1;
				VIP_connect += x(j, 18);
			}
		}
		if (k_connect != 0)
			gamma[i] = VIP_connect / k_connect;
		else
			gamma[i] = 0;

		GABA[i] = GABA_0 + x(i, 20);
		if (t == 0)
		{
			for (size_t ii = i; ii < N; ii++)
				GABA[ii] = GABA[i];
		}
		int mu_connect = 0;
		double GABA_connect = 0;
		for (size_t j = 0; j < N; j++)
		{
			if (a(j, i) == 1)
			{
				mu_connect += 1;
				GABA_connect += GABA[j];
			}
		}
		if (mu_connect != 0)
			delta[i] = GABA_connect / mu_connect;
		else
			delta[i] =0;
			
//		GABA[i] = GABA_0 + x(i, 20);
//		gamma[i] = 0;
//		delta[i] = GABA[i];

		g_K[i] = g_K0 + v_gk * x(i, 0) / (K_gk + x(i, 0));
		g_Ca[i] = v_Ca * pow(x(i, 0), n_Ca) / (K_Ca + pow(x(i, 0), n_Ca));
		g_KCa[i] = v_KCa * pow(x(i, 4), n_KCa) / (K_KCa + pow(x(i, 4), n_KCa));

		Cl_in[i] = Cl_0 + v_Cl1 * x(i, 0) / (K_Cl1 + x(i, 0)) + v_Cl2 * pow(delta[i], n_Cl) / (K_Cl2 + pow(delta[i], n_Cl));

		E_Ca[i] = k_B * T * log(Ca_ex / x(i, 16)) / 2;
		E_GABA[i] = k_B * T * log(Cl_ex / Cl_in[i]);

		K_in[i] = K_ex / exp(E_K / (k_B * T));
		Na_in[i] = Na_ex / exp(E_Na / (k_B * T));
		P_K[i] = v_PK * pow(x(i, 11), npk) / (K_PK + pow(x(i, 11), npk));
		a_v[i] = 4 * P_Ca * x(i, 16) + P_K[i] * K_in[i] + P_Na * Na_in[i] + P_Cl * Cl_ex;
		b_v[i] = P_K[i] * K_in[i] - P_K[i] * K_ex + P_Na * Na_in[i] - P_Na * Na_ex + P_Cl * Cl_ex - P_Cl * Cl_in[i];
		c_v[i] = -1 * (P_K[i] * K_ex + 4 * P_Ca * Ca_ex + P_Na * Na_ex + P_Cl * Cl_in[i]);
		V_rest[i] = k_B * T * log((-b_v[i] + pow((pow(b_v[i], 2) - (4 * a_v[i]) * c_v[i]), 0.5)) / (2 * a_v[i]));
		g_ex[i] = v_ex1 * pow(abs(g_Na * (V_rest[i] - E_Na) * 10e-6), n_ex1) / (K_ex1 + pow(abs(g_Na * (V_rest[i] - E_Na) * 10e-6), n_ex1)) + v_ex2 * pow(x(i, 16), n_ex2) / (K_ex2 + pow(x(i, 16), n_ex2));
		R_m[i] = V_R * V_rest[i] / (K__R + V_rest[i]);
		g_L[i] = 1 / R_m[i];
		V_reset[i] = 4 + V_rest[i];
		I_s[i] = g_Na * E_Na + g_Ca[i] * E_Ca[i] + g_K[i] * E_K + g_KCa[i] * E_K + g_L[i] * E_L - g_ex[i] * E_ex - g_GABA * E_GABA[i];
		R_s[i] = 1 / (g_Na + g_Ca[i] + g_K[i] + g_KCa[i] + g_L[i] - g_ex[i] - g_GABA);
		tau_m[i] = C_m * R_s[i];
		theta[i] = 20 + V_rest[i];
		f_r[i] = -1 * pow((tau_m[i] * log((theta[i] - I_s[i] * R_s[i]) / (V_reset[i] - I_s[i] * R_s[i]))), -1);

		beta[i] = gamma[i] / (K_D + gamma[i]);
		v_k[i] = V_MK * x(i, 16) / (x(i, 16) + K_MK) + V_beta * beta[i] / (beta[i] + K_beta);

		dxdt(i, 0) = (v_sP + C_T * CB_T * x(i, 19) / (K_C + CB_T * x(i, 19))) * pow(x(i, 13), n) / (pow(K_AP, n) + pow(x(i, 13), n)) - v_mP * x(i, 0) / (K_mP + x(i, 0)) - k_dmp * x(i, 0);
		dxdt(i, 1) = v_sC * pow(x(i, 13), n) / (pow(K_AC, n) + pow(x(i, 13), n)) - v_mC * x(i, 1) / (K_mC + x(i, 1)) - k_dmc * x(i, 1);
		dxdt(i, 2) = v_sB * pow(K_IB, m) / (pow(K_IB, m) + pow(x(i, 13), m)) - v_mB * x(i, 2) / (K_mB + x(i, 2)) - k_dmb * x(i, 2);
		dxdt(i, 3) = k_sP * x(i, 0) - V_1P * x(i, 3) / (K_p + x(i, 3)) + V_2P * x(i, 5) / (K_dp + x(i, 5)) + k_4 * x(i, 7) - k_3 * x(i, 3) * x(i, 4) - k_dn * x(i, 3);
		dxdt(i, 4) = k_sC * x(i, 1) - V_1C * x(i, 4) / (K_p + x(i, 4)) + V_2C * x(i, 6) / (K_dp + x(i, 6)) + k_4 * x(i, 7) - k_3 * x(i, 3) * x(i, 4) - k_dnc * x(i, 4);
		dxdt(i, 5) = V_1P * x(i, 3) / (K_p + x(i, 3)) - V_2P * x(i, 5) / (K_dp + x(i, 5)) - v_dPC * x(i, 5) / (K_d + x(i, 5)) - k_dn * x(i, 5);
		dxdt(i, 6) = V_1C * x(i, 4) / (K_p + x(i, 4)) - V_2C * x(i, 6) / (K_dp + x(i, 6)) - v_dCC * x(i, 6) / (K_d + x(i, 6)) - k_dn * x(i, 6);
		dxdt(i, 7) = -V_1PC * x(i, 7) / (K_p + x(i, 7)) + V_2PC * x(i, 9) / (K_dp + x(i, 9)) - k_4 * x(i, 7) + k_3 * x(i, 3) * x(i, 4) + k_2 * x(i, 8) - k_1 * x(i, 7) - k_dn * x(i, 7);
		dxdt(i, 8) = -V_3PC * x(i, 8) / (K_p + x(i, 8)) + V_4PC * x(i, 10) / (K_dp + x(i, 10)) - k_2 * x(i, 8) + k_1 * x(i, 7) - k_7 * x(i, 13) * x(i, 8) + k_8 * x(i, 15) - k_dn * x(i, 8);
		dxdt(i, 9) = V_1PC * x(i, 7) / (K_p + x(i, 7)) - V_2PC * x(i, 9) / (K_dp + x(i, 9)) - v_dPCC * x(i, 9) / (K_d + x(i, 9)) - k_dn * x(i, 9);
		dxdt(i, 10) = V_3PC * x(i, 8) / (K_p + x(i, 8)) - V_4PC * x(i, 10) / (K_dp + x(i, 10)) - v_dPCN * x(i, 10) / (K_d + x(i, 10)) - k_dn * x(i, 10);
		dxdt(i, 11) = k_sB * x(i, 2) - V_1B * x(i, 11) / (K_p + x(i, 11)) + V_2B * x(i, 12) / (K_dp + x(i, 12)) - k_5 * x(i, 11) + k_6 * x(i, 13) - k_dn * x(i, 11);
		dxdt(i, 12) = V_1B * x(i, 11) / (K_p + x(i, 11)) - V_2B * x(i, 12) / (K_dp + x(i, 12)) - v_dBC * x(i, 12) / (K_d + x(i, 12)) - k_dn * x(i, 12);
		dxdt(i, 13) = -V_3B * x(i, 13) / (K_p + x(i, 13)) + V_4B * x(i, 14) / (K_dp + x(i, 14)) + k_5 * x(i, 11) - k_6 * x(i, 13) - k_7 * x(i, 13) * x(i, 8) + k_8 * x(i, 15) - k_dn * x(i, 13);
		dxdt(i, 14) = V_3B * x(i, 13) / (K_p + x(i, 13)) - V_4B * x(i, 14) / (K_dp + x(i, 14)) - v_dBN * x(i, 14) / (K_d + x(i, 14)) - k_dn * x(i, 14);
		dxdt(i, 15) = -k_8 * x(i, 15) + k_7 * x(i, 13) * x(i, 8) - v_dIN * x(i, 15) / (K_d + x(i, 15)) - k_dn * x(i, 15);

		dxdt(i, 16) = v_v0 * pow(x(i, 11), nv0) / (K_v0 + pow(x(i, 11), nv0)) + v_1 * beta_IP3 - (v_kk * pow(x(i, 4), nkk) / (K_kk + pow(x(i, 4), nkk))) * pow(x(i, 16), v) - V_M2 * pow(x(i, 16), n1) / (pow(K_2, n1) + pow(x(i, 16), n1)) + V_M3 * (pow(x(i, 17), m1) / (pow(K_R, m1) + pow(x(i, 17), m1))) * pow(x(i, 16), p1) / (pow(K_A, p1) + pow(x(i, 16), p1)) + k_f * x(i, 17);
		dxdt(i, 17) = V_M2 * pow(x(i, 16), n1) / (pow(K_2, n1) + pow(x(i, 16), n1)) - V_M3 * (pow(x(i, 17), m1) / (pow(K_R, m1) + pow(x(i, 17), m1))) * pow(x(i, 16), p1) / (pow(K_A, p1) + pow(x(i, 16), p1)) - k_f * x(i, 17);
		dxdt(i, 18) = v_VIP * pow(f_r[i], n_VIP) / (K_VIP + pow(f_r[i], n_VIP)) - k_dVIP * pow(x(i, 18), n_dVIP);
		dxdt(i, 19) = (v_P / CB_T) * ((v_k[i] / v_P) * (1 - x(i, 19)) / (K_1 + 1 - x(i, 19)) - x(i, 19) / (K__2 + x(i, 19)));
		dxdt(i, 20) = v_GABA * pow(f_r[i], n_GABA) / (K_GABA + pow(f_r[i], n_GABA)) - k_dGABA * pow(GABA[i], n_dGABA);
	}
}


void write_SCN(const state_type& x, const double t)
{
	ofstream ofile;
	ofile.open("firing_model_network.csv", ios::app);
	ofile << t << '\t' << x(0, 0) << '\t' << f_r[0] << '\t' << x(5, 0) << '\t' << f_r[5] << '\t' << x(15, 0) << '\t' << f_r[15] << endl;
	ofile.close();
}

int main(int argc, char** argv)
{
	ofstream ofile;

	ofile.open("firing_model_network.csv");
	ofile.close();

	construct_connection(p_link, a);

	state_type x(N, 21);
	for (size_t i = 0; i < N; i++)
	{
	    x(i, 0) = 2.80;
		x(i, 1) = 2.00;
		x(i, 2) = 7.94;
		x(i, 3) = 0.40;
		x(i, 4) = 12.0;
		x(i, 5) = 0.13;
		x(i, 6) = 9.00;
		x(i, 7) = 1.26;
		x(i, 8) = 0.16;
		x(i, 9) = 0.20;
		x(i, 10) = 0.091;
		x(i, 11) = 2.41;
		x(i, 12) = 0,48;
		x(i, 13) = 1.94;
		x(i, 14) = 0.32;
		x(i, 15) = 0.05;
		x(i, 16) = 0.10;
		x(i, 17) = 0.10;
		x(i, 18) = 0.00;
		x(i, 19) = 0.12;
		x(i, 20) = 0.00;
	}

	integrate(scn_firing, x, 0.0, 500.0, 0.1, write_SCN);
//	runge_kutta4< state_type > rk;
//	integrate_const(rk, scn_firing, x, 0.0, 500.0, 0.1, write_SCN);

	return 0;
}
