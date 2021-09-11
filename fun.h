#ifndef FUN_H
#define FUN_H

#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

//#include "kuramoto.h"

using namespace std;

class Kuramoto
{
public:
    double a[N][N];
    double peak[N][5];

    int ind[N];

    double e[N];
    double interaction[N];
    double eta[N];
    double theta[N];
    double peak1[N];
    double peak2[N];

    double R_[3] = {0, 0, 0};
    double R[3] = {0, 0, 0};
    double P[3] = {0, 0, 0}; //periods
    int times = 0;

    int t1;
    int t2;

    int light;

    pair<int, int> index_to_coordinate(int i);
    double dist(int i, int j);
    void construct_connection();
    void rewiring();
    void syn_degree();
    double delta_theta();
    void period();
    void find_peak(const double t);
    double phase(double pha_pi);
    void initialization();

    //    int  func (double t, const double y[], double f[], void *params);

    void write_kuramoto(const double y[], const double t);
    void write_peaks();
    void write_connect();
    void write_results();
};

int func(double t, const double y[], double f[], void *params);

#endif