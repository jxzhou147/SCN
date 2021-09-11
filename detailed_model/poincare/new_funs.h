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

// resluts analysis
void find_peak(const double t, double (*peak)[5], double (*peak_now)[2], double (*peak_past)[2], double (*peak_old)[2], const double y[], int ind[]);
double period(const double (*peak)[5]);
double period_one(const double (*peak)[5]);
double syn_degree(const double p, const double (*peak)[5]);
int ifEntrainFun(const double (*peak)[5], const double tcycle);

// write files
void write_curve(const double y[], const double t, int light);
void write_rt(const double t, const double p, const double r);
void write_results(const double p, const double r, const double p_random);
void write_er(const double tcycle, const double K_f, const double p, int ifEntrain);
void write_period_N(const double (*peak)[5]);

#endif