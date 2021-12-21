#ifndef NFUNS_H
#define NFUNS_H

#include "parameter.h"

struct Params
{
    double **a;
    double *ed;
    double *eta;
    double daylength;
    double tcycle;
    double K;
    double K_f;
};

// odes
int poincare(double t, const double y[], double f[], void *params);

// network functions
void construct_connection_random(double **a, double p_random, int tid);
int connection_number(double **a);
void dynamicWireAndCut(double **a, const double p_cut, int tid);
void construct_connection_triVal(double **a, double p_random, double p_8, int tid);
void dynamic_triVal(double **a, const double p_cut, int tid);
void construct_connection_norm(double **a, double p_random, double p_sigma, int tid);
void dynamic_norm(double **a, const double p_cut, double p_sigma, int tid);
void construct_connection_uniform(double **a, double p_random, int tid);
void dynamic_uniform(double **a, const double p_cut, int tid);
void dynamic_uniform_random(double **a, int tid);
void construct_connection_powerLaw(double **a, double p_random, int tid);
void dynamic_uniform_powerLaw(double **a, int tid);
void construct_connection_vd(double **a, double p_link, double d, double p_vd, double p_dv, int tid);
void dynamicNet(double **a, double p_cut);
void construct_connection_sw_triVal(double **a, int m, double p_random, double p_8, int tid);
void construct_connection_sf_triVal(double **a, int m, double p_8, int tid);
void construct_connection_random_undirected(double **a, double p_random, int tid);
void dynamicWireAndCut_undirected(double **a, const double p_cut, int tid);
void dynamic_triVal_random(double **a, double p_8, int tid);
int ifConnected(double **a);
int part_num(double **a);
void net_component(double **a, int *component);

// resluts analysis
void find_peak_t(const double t, double (*peak)[90], double (*peak_now)[2], double (*peak_past)[2], double (*peak_old)[2], const double y[], int ind[]);
void find_peak(const double t, double (*peak)[5], double (*peak_now)[2], double (*peak_past)[2], double (*peak_old)[2], const double y[], int ind[]);
double period(const double (*peak)[5]);
double period_ave(const double (*peak)[5]);
double period_std(const double (*peak)[5]);
double period_one(const double (*peak)[5]);
double syn_degree(const double p, const double (*peak)[5]);
int ifEntrainFun(const double (*peak)[5], const double tcycle);

// write files
void write_curve(const double y[], const double t, int light);
void write_rt(const double t, const double p, const double r);
void write_results(const double r, const double p_cut, double t_re);
void write_results_period(const double period_a, const double period_s, const double r, const double K, int *component);
void write_er(const double tcycle, const double K_f, const double p, int ifEntrain);
void write_period_N(const double (*peak)[5]);
void write_connect(double **a);
void write_net(double p_random, int connect_num, int part_num);

#endif