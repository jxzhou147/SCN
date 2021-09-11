#ifndef FUNS_H
#define FUNS_H

#include "parameter.h"

struct Params
{
    double **a;
    double *ed;
    double *eta;
    double daylength;
    double K;
    double K_f;
};

struct Enrange
{
    double daylength;
    double K_f;
};

void write_curve(const double y[], const double t, int light);
int poincare(double t, const double y[], double f[], void *params);
double firingRate(double M_P, double C_C, double B_C, double Ca, double delta, double t);
void construct_connection(double (*a_vip)[N], double (*a_gaba)[N], double p_link, double d);
void construct_connection_vip(double (*a_vip)[N], double (*a_gaba)[N], double p_link, double d, double p_vd);
//void construct_connection_vd(double (*a_vip)[N], double (*a_gaba)[N], double p_link, double d, double p_vd, int tid);
void construct_connection_vd(double **a, double p_link, double d, double p_vd, double p_dv, int tid);
void connection_number(double **a_vip, double **a_gaba, int connect_num[]);
//void dynamicNet(double **a_vip, double **a_gaba, double **a_vip_static, double **a_gaba_static, const double t, double delay, double p_cut, double day);
void dynamicNet(double **a, double p_cut);
//void connection_number(const double (*a_vip)[N], const double (*a_gaba)[N], int connect_num[]);
int scn_firing(double t, const double y[], double f[], void *params);
void write_MP(const double y[], const double t, int light);
void write_connect(const double (*a)[N]);
void write_fr(const double f_r[], const double t, int light);
void write_input(const double p_intput, const double t, int light);
void write_results(const double P[], const double P_var[], const double R[], const double p_vd, const int connect_num[], const double beta_light, double p_cut, double day);
void write_peaks(const double peak_ave[], const int peak_ave_num);
void write_peaks_vd(const double peak_ave_v[], const double peak_ave_d[], const int peak_ave_num);
void write_peaks_vd_p(const double peak_ave_v[], const double peak_ave_d[], const int peak_ave_num, const double p_vd, const double p_dv, const double daylength);
void write_peaks_p(const double peak_ave[], const int peak_ave_num, const int tid, const double pulse);
void write_er(const double daylength, const double P[], int ifEntrain);
void write_rt(const double t, const double P[], const double P_var[], const double R[]);
void write_period_N(const double (*peak)[5]);
void find_peak_ave(const double t, double peak_ave[], double peak_now[], double peak_past[], double peak_old[], const double y[], int ind_ave[]);
void find_peak_ave_vd(const double t, double peak_ave_v[], double peak_now_v[], double peak_past_v[], double peak_old_v[], double peak_ave_d[], double peak_now_d[], double peak_past_d[], double peak_old_d[], const double y[], int ind_ave_v[], int ind_ave_d[]);
void find_peak(const double t, double (*peak)[5], double (*peak_now)[2], double (*peak_past)[2], double (*peak_old)[2], const double y[], int ind[]);
void period(double P[], double P_var[], const double (*peak)[5]);
void period_ave(double P[], double P_var[], const double (*peak)[5]);
void period_N1(double P[], double P_var[], const double (*peak)[5]);
void syn_degree(double R[], const double P[], const double (*peak)[5]);
void syn_degree_N1(double R[], const double P[], const double (*peak)[5]);
double real_root_no(double x, double n);
double real_root_transfer(double x, double n);
double real_root_zero(double x, double n);

#endif