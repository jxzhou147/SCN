//parameters in the model

const int N = 25;          //neurons' number
const int N0 = 5;			//N0^2 = N
const double p_vip = 0.2;			//VIP neurons' percentage
const double p_link = 0.1;	//probability of long range connections

const double k_1 = 0.45;
const double k_2 = 0.2;
const double k_3 = 0.4;
const double k_4 = 0.2;
const double k_5 = 0.4;
const double k_6 = 0.2;
const double k_7 = 0.5;
const double k_8 = 0.1;
const double K_AP = 0.6;
const double K_AC = 0.6;
const double K_IB = 2.2;
const double k_dmb = 0.01;
const double k_dmc = 0.01;
const double k_dmp = 0.01;
const double k_dn = 0.01;
const double k_dnc = 0.12;
const double K_d = 0.3;
const double K_dp = 0.1;
const double K_p = 0.1;
const double K_mB = 0.4;
const double K_mC = 0.4;
const double K_mP = 0.31;
const double k_sB = 0.12;
const double k_sC = 1.6;
const double k_sP = 0.6;
const int m = 2;
const int n = 4;
const double V_1B = 0.5;
const double V_1C = 0.6;
const double V_1P = 0.4;
const double V_1PC = 0.4;
const double V_2B = 0.1;
const double V_2C = 0.1;
const double V_2P = 0.3;
const double V_2PC = 0.1;
const double V_3B = 0.5;
const double V_3PC = 0.4;
const double V_4B = 0.2;
const double V_4PC = 0.1;
const double V_phos = 0.4;
const double v_dBC = 0.5;
const double v_dBN = 0.6;
const double v_dCC = 0.7;
const double v_dIN = 0.8;
const double v_dPC = 0.7;
const double v_dPCC = 0.7;
const double v_dPCN = 0.7;
const double v_mB = 0.8;
const double v_mC = 1.0;
const double v_mP = 1.1;
const double v_sB = 1.0;
const double v_sC = 1.1;
const double v_sP = 0.94;

const double n1 = 2.2;
const double m1 = 6;
const double p1 = 4.2;
const double E_K = -97;
const double k_B = 8.61585492e-2; // R/F * 10e3
const double T = 310.15;
const double g_K0 = 9.7;
const double v_gk = 10;
const double K_gk = 10;
const double g_Na = 36;
const double E_Na = 45;
const double v_kk = 3.3;
const double K_kk = 0.02;
const double nkk = 0.1;
const double v_v0 = 0.09;
const double K_v0 = 4.5;
const double nv0 = 4.5;
const double v = 2;
const double v_1 = 0.0003;
const double beta_IP3 = 0.5;
const double V_M2 = 149.5;
const double K_2 = 5;
const double V_M3 = 400;
const double K_R = 3;
const double K_A = 0.67;
const double k_f = 0.001;
const double Ca_ex = 5;
const double v_Ca = 12.3;
const double K_Ca = 22;
const double n_Ca = 2.2;
const double v_KCa = 3;
const double K_KCa = 0.16;
const double n_KCa = -1;
const double E_L = -29;
const double GABA_0 = 0.1;
const double v_GABA = 0.5;
const double K_GABA = 20;
const double g_GABA = 12.3;
const double n_GABA = 1.9;
const double k_dGABA = 0.5;
const double n_dGABA = 0.2;
const double Cl_0 = 1;
const double v_Cl1 = 15.5;
const double v_Cl2 = 19;
const double K_Cl1 = 4;
const double K_Cl2 = 1;
const double n_Cl = -0.2;
const double Cl_ex = 114.5;
const double v_ex1 = 101;
const double K_ex1 = 574.05;
const double n_ex1 = 2.5;
const double v_ex2 = 3.5;
const double K_ex2 = 1;
const double n_ex2 = -1;
const double E_ex = 0;
const double P_Ca = 0.05;
const double P_Na = 0.036;
const double P_Cl = 0.3;
const double K_ex = 1;
const double Na_ex = 145;
const double v_PK = 1.9;
const double K_PK = 1;
const double npk = -2;
const double V_R = 0.41;
const double K__R = 34;
const double v_VIP = 0.5;
const double K_VIP = 20;
const double n_VIP = 1.9;
const double k_dVIP = 0.5;
const double n_dVIP = 0.2;
const double V_MK = 5;
const double K_MK = 4.2;
const double V_beta = 2.1;
const double K_beta = 1.8;
const double C_T = 1.5;
const double K_C = 0.15;
const double K_D = 0.06;
const double C_m = 6;

const double CB_T = 1;
const double v_P = 1;
const double K_1 = 0.01;
//const double K__2 = 0.01;
const double K__2 = 0.1;
