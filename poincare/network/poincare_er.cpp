#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <pthread.h>

#include "poincare.h"

using namespace std;
using namespace std::chrono;

pthread_rwlock_t rwlock;
double compute_tau = 200;
int tnum = 0;
int numEntrain = 0;
int ifEntrain = 0;

struct Params
{
    double (*a)[N];
    double *e;
    double *interaction;
    double *eta;
    double *theta;
    int light;
    double daylength;
    double K;
    double K_f;
    double t0;
};

struct Enrange
{
    double daylength;
    double K_f;
};

pair<int, int> index_to_coordinate(int i)
{
    pair<int, int> coor(i / N0, i % N0);
    return coor;
}

double dist(int i, int j)
{
    double xi = (double)(index_to_coordinate(i).first);
    double yi = (double)(index_to_coordinate(i).second);
    double xj = (double)(index_to_coordinate(j).first);
    double yj = (double)(index_to_coordinate(j).second);

    double dis = sqrt(pow((xi - xj), 2) + pow((yi - yj), 2));
    return dis;
}

void random_wiring(int i, double (*a)[N], double d, int Nh, int Ne)
{
    //srand((unsigned)time(NULL));

    int wiring_node = 0;
    double wiring_max = 0.0;
    double wiring_rand = 0;

    for (size_t j = Nh; j < Ne; j++)
    {
        if ((a[i][j] == 0) & (dist(i, j) > d))
        {
            wiring_rand = (rand() % (N_rand + 1) / (float)(N_rand + 1));
            if (wiring_rand > wiring_max)
            {
                wiring_node = j;
                wiring_max = wiring_rand;
            }
        }
    }

    a[i][wiring_node] = 1;
}

void weighted_wiring(int i, const double y[], double (*a)[N], double tau)
{
    //srand((unsigned)time(NULL));

    double dif;
    //double tau = 0.3;
    int wiring_node = 0;
    double wiring_max = 0.0;
    double wiring_rand = 0;

    for (size_t j = 0; j < N; j++)
    {
        if ((a[i][j] == 0) & (i != j))
        {
            dif = exp(fmod(abs(y[i] - y[j]), 2 * PI) / tau);
            wiring_rand = (rand() % (N_rand + 1) / (float)(N_rand + 1)) / dif;
            if (wiring_rand > wiring_max)
            {
                wiring_node = j;
                wiring_max = wiring_rand;
            }
        }
    }

    a[i][wiring_node] = 1;
}

void construct_connection(double (*a)[N], double p_link, double d)
{
    //srand((unsigned)time(NULL));

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;

            if ((dist(i, j) <= d) & (i != j))
            {
                a[i][j] = 1;

                if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link)
                {
                    a[i][j] = 0;
                    random_wiring(i, a, d, 0, N);
                }
            }

            //if ((i < N_v) & (j > N_v)) //ventral to dorsal SCN's interaction is replusive (a_ij: i to j)
            //  a[i][j] = -1 * a[i][j];
        }
    }
}

void construct_connection_new(double (*a)[N], double p_link_v, double p_link_d, double d)
{
    //srand((unsigned)time(NULL));
    for (size_t i = 0; i < N_v; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;
            if (j < N_v)
            {
                if ((dist(i, j) <= d) & (i != j))
                {
                    a[i][j] = 1;
                }
            }
            if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link_v)
                a[i][j] = 1;
        }
    }

    for (size_t i = N_v; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;
            if (j >= N_v)
            {
                if ((dist(i, j) <= d) & (i != j))
                {
                    a[i][j] = 1;
                }
            }
            if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link_d)
                a[i][j] = 1;
        }
    }
}

void rewiring_new(const double y[], double (*a)[N], double tau, double p_re, double d)
{
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if ((a[i][j] == 1) & (dist(i, j) > d) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_re))
            {
                a[i][j] = 0;
                random_wiring(i, a, d, 0, N);
            }
        }
    }
}

void rewiring(const double y[], double (*a)[N], double tau, double p_re, double d)
{
    //srand((unsigned)time(NULL));

    double dif;
    //double tau = 0.3;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            //dif = exp(fmod(abs(y[i] - y[j]), 2 * PI) / tau); //weights of rewiring and cutting
            dif = 1;
            //if ((a[i][j] == 1) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < (p_re)))
            if ((a[i][j] == 1) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < (p_re * dif)))
            {
                a[i][j] = 0;
                //weighted_wiring(i, y, a, tau);
                random_wiring(i, a, d, 0, N);
            }

            //if ((i < N_v) & (j > N_v)) //ventral to dorsal SCN's interaction is replusive (a_ij: i to j)
            //  a[i][j] = -1 * a[i][j];

            //			if ((a[i][j] == 0) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < (p_re / dif)))
            //				a[i][j] = 1;
            /*
            if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_re)
            {
                if (a[i][j] == 1)
                {
                    a[i][j] = 0;
                    random_wiring(i, a);
                }
            }*/
        }
    }
}

int connection_number(const double (*a)[N])
{
    int num = 0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] == 1)
                num++;
        }
    }

    return num;
}

void syn_degree(double R_[], const double theta[])
{
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;

    for (size_t i = 0; i < N_v; i++)
    {
        sum1 += cos(theta[i]);
        sum2 += sin(theta[i]);
    }

    R_[0] += sqrt(pow(sum1, 2) + pow(sum2, 2));

    for (size_t i = N_v; i < N; i++)
    {
        sum3 += cos(theta[i]);
        sum4 += sin(theta[i]);
    }

    R_[1] += sqrt(pow(sum3, 2) + pow(sum4, 2));

    R_[2] += sqrt(pow(sum1 + sum3, 2) + pow(sum2 + sum4, 2));
}

double delta_theta_old(const double (*peak)[5])
{
    double tv[5] = {0, 0, 0, 0, 0};
    double ts[5] = {0, 0, 0, 0, 0};
    double delta_the = 0;

    //for (size_t i = 0; i < 5; i++)
    for (size_t i = 0; i < 1; i++)
    {
        for (size_t j = 0; j < N_v; j++)
            tv[i] += peak[j][i];
        for (size_t j = N_v; j < N; j++)
            ts[i] += peak[j][i];

        //delta_the += abs((tv[i] / N_v - ts[i] / (N - N_v)));
        delta_the += abs((peak[0][i] - peak[N_v][i]));
    }

    //delta_the = delta_the / 5;
    delta_the = delta_the / 1;

    return delta_the;
}

double delta_theta(const double (*peak)[5])
{
    double tv[5] = {0, 0, 0, 0, 0};
    double ts[5] = {0, 0, 0, 0, 0};
    double tvs[5] = {0, 0, 0, 0, 0};
    double tvsi;
    double delta_the = 0;

    for (size_t i = 0; i < 5; i++)
    {
        for (size_t j = 0; j < N_v; j++)
        {
            for (size_t k = N_v; k < N; k++)
            {
                tvsi = abs((peak[j][i] - peak[k][i]));
                if (tvsi > 12)
                    tvsi = abs(24 - tvsi);
                tvs[i] += tvsi;
            }
        }
        delta_the += (tvs[i] / (N_v * (N - N_v)));
    }

    delta_the = delta_the / 5;

    return delta_the;
}

void period(double P[], const double (*peak)[5])
{
    for (size_t i = 0; i < N_v; i++)
        P[0] += ((peak[i][4] - peak[i][0]) / 4);
    P[0] /= N_v;

    for (size_t i = N_v; i < N; i++)
        P[1] += (peak[i][4] - peak[i][0]) / 4;
    P[1] /= (N - N_v);

    P[2] = (P[0] * N_v + P[1] * (N - N_v)) / N;
}

void find_peak(const double t, double (*peak)[5], double peak1[], double peak2[], const double theta[], int ind[])
{
    int ind_min = 5;
    for (size_t i_min = 0; i_min < N; i_min++)
    {
        if (ind_min > ind[i_min])
            ind_min = ind[i_min];
    }

    if (ind_min < 5)
    {
        for (size_t i = 0; i < N; i++)
            peak1[i] = theta[i];

        for (size_t i_p = 0; i_p < N; i_p++)
        {
            if (((peak2[i_p] - peak1[i_p]) > 6) & (ind[i_p] < 5))
            {
                peak[i_p][ind[i_p]] = t;
                ind[i_p] = ind[i_p] + 1;
            }
        }

        for (size_t i = 0; i < N; i++)
            peak2[i] = theta[i];
    }
}

void find_peak_ave(const double t, double peak_ave[], double &peak1_ave, double &peak2_ave, const double theta_ave, int &ind_ave)
{
    peak1_ave = theta_ave;

    if (((peak2_ave - peak1_ave) > 5) & (ind_ave < 40))
    {
        peak_ave[ind_ave] = t;
        ind_ave++;
    }

    peak2_ave = theta_ave;
}

double phase(double pha_pi)
{
    //double pha = 2.0 * PI / 24.0; A HUGE MISTAKE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double pha = pha_pi * 24.0 / (2.0 * PI);
    if ((pha > 3.0) & (pha < 9.0))
        return 0.0;
    else
    {
        if (pha > 9.0)
            return -1.0 * sin((pha - 9) * 2 * PI / 18.0);
        else
            return -1.0 * sin((pha + 24 - 9) * 2 * PI / 18.0);
    }
}

double phase_light(double pha_pi, double light)
{
    double t0 = 2;
    double pha;
    if (pha_pi > 0)
        pha = pha_pi * 24.0 / (2.0 * PI);
    else
        pha = (pha_pi + 2 * PI) * 24.0 / (2.0 * PI);

    if ((pha > t0) & (pha < (light - t0))) //// keep the morning 3 hours and dusk 3 hours having effect of entrainment
        return 0.0;
    else
    {
        if (pha > (light - t0))
            return -1 * sin((pha - (light - t0)) * 2 * PI / (24 - (light - 2 * t0)));
        else
            return -1 * sin((pha + 24 - (light - t0)) * 2 * PI / (24 - (light - 2 * t0)));
    }
}

double peak_theta(const double theta[], int flag) ///// flag = 0: SCN; flag = 1: ventral SCN; flag = 2: dorsal SCN
{
    int num = 64;
    int count[num] = {0}; // -3.2, -3.1, -3.0, ..., 0, 0.1, ... , 3.0, 3.1
    int temp;
    int imin, imax;

    if (flag == 0)
    {
        imin = 0;
        imax = N;
    }
    else if (flag == 1)
    {
        imin = 0;
        imax = N_v;
    }
    else
    {
        imin = N_v;
        imax = N;
    }

    for (size_t i = imin; i < imax; i++)
    {
        temp = (int)floor(theta[i] * 10.0) + 32;
        count[temp]++;
    }

    ////// find the biggest number of count
    int big = 0;
    int big_ind = 0;
    for (size_t i = 0; i < num; i++)
    {
        if (count[i] > big)
        {
            big = count[i];
            big_ind = i;
        }
    }

    return ((double)big_ind / 10.0 - 3.14);
}

//////// 2N-dim, 0:N-1 -> x, N:2N-1 -> y
int kuramoto(double t, const double y[], double f[], void *params)
{
    Params *param;
    param = (struct Params *)params;
    for (size_t i = 0; i < N; i++)
    {
        param->interaction[i] = 0;
        if (param->e[i] != 0)
        {
            for (size_t j = 0; j < N; j++)
                param->interaction[i] += param->a[j][i] * y[j] / param->e[i];
        }
    }
    /*
    if ((t > param->t0) & (t <= (param->t0 + 10)))
        param->light = 1;
    else
        param->light = 0;
*/
    //if (fmod(t, 24.0) < param->daylength)
    //  param->light = 1;
    //else
    //  param->light = 0;

    double delay = 0;
    if (t > param->t0)
        delay = 0;
    if (fmod(t + delay, param->daylength) < (param->daylength / 2.0))
        param->light = 1;
    else
        param->light = 0;

    for (size_t i = 0; i < N_v; i++)
    {
        //dxdt[i] = omega * eta[i] + K * interaction[i] + K_f * sin(Omega * t - x[i]);
        //f[i] = omega_v * param->eta[i] + K * param->interaction[i] + K_f * param->light * phase(param->theta[i]);
        //f[i] = omega_v * param->eta[i] + K * param->interaction[i] + K_f * sin(Omega * t - y[i]);
        //f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_v * param->eta[i] + param->K * param->interaction[i] + param->K_f * param->light * phase_light(param->theta[i], param->daylength);
        f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_v * param->eta[i] + param->K * param->interaction[i] + param->K_f * param->light;
        //f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_v * param->eta[i];
        f[i + N] = gam * y[i + N] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) + y[i] * omega_v * param->eta[i];
        //f[i + N] = gam * y[i + N] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) + y[i] * omega_v * param->eta[i];

        param->theta[i] = atan2(y[i + N], y[i]);
    }

    for (size_t i = N_v; i < N; i++)
    {
        f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_d * param->eta[i] + param->K * param->interaction[i];

        //f[i] = gam * y[i] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) - y[i + N] * omega_d * param->eta[i];
        f[i + N] = gam * y[i + N] * (a - sqrt(pow(y[i], 2) + pow(y[i + N], 2))) + y[i] * omega_d * param->eta[i];
        param->theta[i] = atan2(y[i + N], y[i]);
    }

    return GSL_SUCCESS;
}

void write_kuramoto(const double y[], const double t, const double theta[], int light)
{
    //double phase_sun = fmod(Omega * t, 2 * PI);
    //double delta_phase[N];
    //double ave_phase = 0;

    //for (size_t i = 0; i < N; i++)
    //  ave_phase += theta[i];
    //ave_phase = ave_phase / double(N); // this average is not perfect matching each other, so that the delta phase has some big peaks
    /*
    for (size_t i = 0; i < N; i++)
    {
        delta_phase[i] = theta[i] - ave_phase;
        if (delta_phase[i] > PI)
            delta_phase[i] = 2 * PI - delta_phase[i];
        if (delta_phase[i] < (-1.0 * PI))
            delta_phase[i] = 2 * PI + delta_phase[i];
    }
*/
    ofstream ofile;
    ofile.open("kuramoto.csv", ios::app);
    //ofile << t << '\t' << light << '\t' << y[0] << '\t' << y[5] << '\t' << y[15] << '\t' << y[30] << '\t' << y[40] << '\t' << y[50] << '\t' << y[70] << '\t' << y[80] << '\t' << y[90] << '\t' << y[95] << endl;
    //ofile << t << '\t' << light << '\t' << theta[0] << '\t' << theta[5] << '\t' << theta[15] << '\t' << theta[30] << '\t' << theta[40] << '\t' << theta[50] << '\t' << theta[70] << '\t' << theta[80] << '\t' << theta[90] << '\t' << theta[95] << endl;
    /*
    ofile << t << '\t';
    for (int i = 0; i < N; i++)
    {
        ofile << theta[i] << '\t';
    }
    ofile << endl;
*/
    ofile << t << '\t' << light << '\t' << peak_theta(theta, 1) << '\t' << peak_theta(theta, 2) << '\t' << peak_theta(theta, 0) << endl;
    ofile.close();
}

void write_peaks(const double (*peak)[5])
{
    ofstream ofile;
    ofile.open("peaks.csv", ios::app);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 5; j++)
            ofile << peak[i][j] << '\t';
        ofile << endl;
    }
    ofile.close();
}

void write_connect(const double (*a)[N])
{
    ofstream ofile;
    ofile.open("connection.csv", ios::app);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
            ofile << a[i][j] << '\t';
        ofile << endl;
    }
    ofile.close();
}

void write_results(const double t, const double (*a)[N], const double R[], const double (*peak)[5], const double peak_ave[], const double peak_ave_v[], const double peak_ave_d[], const double on_sun[], const double P[], const double daylength, const double tau, const double t_re, double purt, double p_link, double p_link_v, double p_link_d, double p_re, double d, double sigma, double K, double K_f, int flag)
{
    int ifEntrainment = 0;
    if (ifEntrain > 24)
        ifEntrainment = 1;

    ofstream ofile;
    ofile.open("results_er.csv", ios::app);
    //ofile << t << '\t' << N << '\t' << connection_number(a) << '\t' << t_re << '\t' << d << '\t' << p_link << '\t' << p_re << '\t' << Omega << '\t' << daylength << '\t' << K_f << '\t' << R[0] << '\t' << R[1] << '\t' << R[2] << '\t' << delta_theta(peak) << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << '\t' << tau << '\t' << flag << endl;
    //ofile << t << '\t' << tau << '\t' << t_re << '\t' << daylength << '\t' << R[0] << '\t' << R[1] << '\t' << R[2] << '\t' << delta_theta(peak) << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << endl;
    //ofile << t << '\t' << p_re << '\t' << R[0] << '\t' << R[1] << '\t' << R[2] << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << endl;
    ofile << p_re << '\t' << K_f << '\t' << daylength << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << '\t' << ifEntrainment << endl;
    //ofile << sigma;
    //for (size_t i = 0; i < N; i++)
    //  ofile << '\t' << ((peak[i][4] - peak[i][0]) / 4);
    //ofile << endl;
    //for (size_t i = 0; i < 40; i++)
    //  ofile << on_sun[i] << '\t' << peak_ave_v[i] << '\t' << peak_ave_d[i] << '\t' << peak_ave[i] << endl;
    ofile.close();
}

void *threadfunc(void *arg)
{
    //double daylength = *(double *)arg;
    Enrange *enrange;
    enrange = (struct Enrange *)arg;
    double K_f = enrange->K_f;
    double daylength = enrange->daylength;
    cout << K_f << '\t' << daylength << endl;
    double p_re = 0.01;
    //double K_f = 0.15;
    int tid = tnum++;
    double t_re = 0.1;
    double sigma = 0.05;
    //double sigma = 0;
    double K = 0.15;
    //double daylength = 24;
    double d = 1;
    double p_link = 0;
    double purt = 0;
    double tau = 3.0;
    double a[N][N];
    double peak[N][5];
    int ind[N];
    double p_link_v = 0.001;
    double p_link_d = p_link_v;
    //double p_link_d = 0.01;

    double peak_ave[40];
    int ind_ave = 0;
    double peak1_ave, peak2_ave;

    double peak_ave_v[40];
    int ind_ave_v = 0;
    double peak1_ave_v, peak2_ave_v;

    double peak_ave_d[40];
    int ind_ave_d = 0;
    double peak1_ave_d, peak2_ave_d;

    double on_sun[40];
    int ind_sun = 0;
    int on_sun1, on_sun2;

    double e[N];
    double interaction[N];
    double eta[N];
    double y[2 * N];
    double theta[N];
    double peak1[N] = {0};
    double peak2[N] = {0};

    double R_[3] = {0, 0, 0};
    double R[3] = {0, 0, 0};
    double P[3] = {0, 0, 0}; //periods
    int times = 0;

    int t1;
    int t2;

    int light;
    //double daylength;
    int write_light;

    double daylength_0 = daylength;

    //for (int i = 0; i < 25; i++)
    //{
    //	daylength = double(i);

    srand(tid);

    for (size_t i = 0; i < 3; i++)
    {
        R_[i] = 0;
        R[i] = 0;
        P[i] = 0; //periods
    }

    //construct_connection(a, p_link, d);
    construct_connection_new(a, p_link_v, p_link_d, d);

    for (size_t i = 0; i < N; i++)
    {
        e[i] = 0;
        for (size_t j = 0; j < N; j++)
            e[i] += abs(a[j][i]);

        ind[i] = 0; // initialize the ind array

        for (size_t j = 0; j < 5; j++)
            peak[i][j] = 0;
    }

    normal_distribution<double> norm{1, sigma};
    //random_device rd;
    //default_random_engine rng{rd()};
    default_random_engine rng; ////// the norm distribution is same for different threads
    for (size_t i = 0; i < N; i++)
    {
        eta[i] = norm(rng);
        if (((rand() % (N_rand + 1) / (float)(N_rand + 1)) < purt) & (i < N_v))
            eta[i] = 0;
    }

    for (size_t i = 0; i < N; i++)
    {
        y[i] = (rand() % (N_rand + 1) / (float)(N_rand + 1)) * 2.0 - 1.0;
        y[i + N] = (rand() % (N_rand + 1) / (float)(N_rand + 1)) * 2.0 - 1.0;
        theta[i] = atan2(y[i + N], y[i]);
    }

    //double t = 0.0, t0 = 600.0, t_dark = 600; // t0 is the time of light off
    double t = 0.0, t0 = 3990.0, t01 = 3000;

    Params param =
        {
            a,
            e,
            interaction,
            eta,
            theta,
            light,
            daylength,
            K,
            K_f,
            t0};

    //const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
    const gsl_odeiv_step_type *T = gsl_odeiv_step_rk2;
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 2 * N);

    gsl_odeiv_system sys = {kuramoto, NULL, 2 * N, &(param)};

    double h = 0.1;
    double y_err[2 * N];
    double dydt_in[2 * N], dydt_out[2 * N];

    double compute_t1 = 0, compute_t2 = 0;

    GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
    while (t < t01)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        compute_t1 = floor(t / compute_tau);
        /*
        if (compute_t1 > 1)
        {
            if (times < 1200)
            {
                times++;
                syn_degree(R_, theta);
            }
        }
        
        if (t > (t01 - 120))
        {
            times++;
            syn_degree(R_, theta);
        }
*/
        ///////// rewiring
        //if (t > 100)
        //{
        t1 = floor(t / t_re);
        if (t1 != t2)
        {
            //rewiring(y, a, tau, p_re, d);
            rewiring_new(y, a, tau, p_re, d);

            for (size_t i = 0; i < N; i++)
            {
                e[i] = 0;
                for (size_t j = 0; j < N; j++)
                    e[i] += abs(a[j][i]);
            }
        }
        t2 = floor(t / t_re);
        //}

        //if (t > (compute_t1 * compute_tau))
        //  find_peak(t, peak, peak1, peak2, theta, ind);

        //if ((compute_t1 != compute_t2))
        //  write_kuramoto(y, t, theta, write_light);
        /*
        if ((compute_t1 != compute_t2) & (compute_t1 > 1))
        {
            R[0] = R_[0] / (N_v * times);
            R[1] = R_[1] / ((N - N_v) * times);
            R[2] = R_[2] / (N * times);

            //period(P, peak);

            pthread_rwlock_wrlock(&rwlock);
            //write_results(t, a, R, peak, peak_ave, peak_ave_v, peak_ave_d, on_sun, P, daylength, tau, t_re, purt, p_link, p_link_v, p_link_d, p_re, d, sigma, K, K_f, 0);
            pthread_rwlock_unlock(&rwlock);

            times = 0;
            for (size_t i = 0; i < 3; i++)
            {
                R_[i] = 0;
                R[i] = 0;
                P[i] = 0; //periods
            }

            for (size_t i = 0; i < N; i++)
            {
                ind[i] = 0; // initialize the ind array

                for (size_t j = 0; j < 5; j++)
                    peak[i][j] = 0;
            }
        }
*/
        compute_t2 = floor(t / compute_tau);
        /*
        if ((t > t0) & (t <= (t0 + 10)))
            write_light = 1;
        else
            write_light = 0;
        */
        double delay = 0;
        if (t > t0)
            delay = 0;
        if (fmod(t + delay, 24.0) < (daylength / 2.0))
            write_light = 1;
        else
            write_light = 0;

        /*
        if (t > (t0 - 300))
        {
            find_peak_ave(t, peak_ave, peak1_ave, peak2_ave, peak_theta(theta, 0), ind_ave);
            find_peak_ave(t, peak_ave_v, peak1_ave_v, peak2_ave_v, peak_theta(theta, 1), ind_ave_v);
            find_peak_ave(t, peak_ave_d, peak1_ave_d, peak2_ave_d, peak_theta(theta, 2), ind_ave_d);

            if (ind_sun < 40)
            {
                on_sun1 = write_light;
                if (on_sun1 < on_sun2)
                {
                    on_sun[ind_sun] = t;
                    ind_sun++;
                }
                on_sun2 = write_light;
            }
        }
        */
        //if ((t > 1800) || (t < 1000))
        //write_kuramoto(y, t, theta, write_light);

        //if ((t < 2000) || (t > t0))
        //{
        //   pthread_rwlock_wrlock(&rwlock);
        // write_kuramoto(y, t, theta, write_light);
        //pthread_rwlock_unlock(&rwlock);
        //}

        if (t > (t01 - 200))
            find_peak(t, peak, peak1, peak2, theta, ind);

        ////// A VERY HUGE MISTAKE !!!!! TO FIX THIS IN OTHER FILES!
        for (size_t i = 0; i < 2 * N; i++)
            dydt_in[i] = dydt_out[i];
        t += h;
    }

    //write_kuramoto(y, t, theta, write_light);

    period(P, peak);
    /*
    R[0] = R_[0] / (N_v * times);
    R[1] = R_[1] / ((N - N_v) * times);
    R[2] = R_[2] / (N * times);
    //cout << times << endl;

    pthread_rwlock_wrlock(&rwlock);
    //write_results(t, a, R, peak, peak_ave, peak_ave_v, peak_ave_d, on_sun, P, daylength, tau, t_re, purt, p_link, p_re, d, sigma, K, K_f, 0);
    write_results(t, a, R, peak, peak_ave, peak_ave_v, peak_ave_d, on_sun, P, daylength, tau, t_re, purt, p_link, p_link_v, p_link_d, p_re, d, sigma, K, K_f, 0);
    pthread_rwlock_unlock(&rwlock);
*/
    pthread_rwlock_wrlock(&rwlock);
    numEntrain++;
    if ((abs(P[0] - daylength) < 0.05) & (abs(P[1] - daylength) < 0.05) & (abs(P[2] - daylength) < 0.05))
        ifEntrain++;
    pthread_rwlock_unlock(&rwlock);

    if (numEntrain == 30)
    {
        write_results(t, a, R, peak, peak_ave, peak_ave_v, peak_ave_d, on_sun, P, daylength, tau, t_re, purt, p_link, p_link_v, p_link_d, p_re, d, sigma, K, K_f, 0);
    }

    //write_peaks(peak);
    //write_connect(a);
    /*
    daylength = 0;

    times = 0;

    for (size_t i = 0; i < 3; i++)
    {
        R_[i] = 0;
        R[i] = 0;
        P[i] = 0; //periods
    }
    for (size_t i = 0; i < N; i++)
    {
        ind[i] = 0; // initialize the ind array

        for (size_t j = 0; j < 5; j++)
            peak[i][j] = 0;
    }

    while (t < t_dark)
    {
        int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
        if (status != GSL_SUCCESS)
            break;

        compute_t1 = floor(t / compute_tau);

        if (t > (compute_t1 * compute_tau))
        {
            if (times < 100)
            {
                times++;
                syn_degree(R_, theta);
            }
        }

        ///////// rewiring
        if (t > 100)
        {
            t1 = floor(t / t_re);
            if (t1 != t2)
            {
                rewiring(y, a, tau, p_re);

                for (size_t i = 0; i < N; i++)
                {
                    e[i] = 0;
                    for (size_t j = 0; j < N; j++)
                        e[i] += abs(a[j][i]);
                }
            }
            t2 = floor(t / t_re);
        }

        if (t > (compute_t1 * compute_tau))
            find_peak(t, peak, peak1, peak2, theta, ind);

        if (compute_t1 != compute_t2)
        {
            R[0] = R_[0] / (N_v * times);
            R[1] = R_[1] / ((N - N_v) * times);
            R[2] = R_[2] / (N * times);

            period(P, peak);

            pthread_rwlock_wrlock(&rwlock);

            write_results(t, a, R, peak, P, daylength_0, tau, t_re, purt, p_link, p_re, sigma, K, K_f, 1);

            pthread_rwlock_unlock(&rwlock);

            times = 0;
            for (size_t i = 0; i < 3; i++)
            {
                R_[i] = 0;
                R[i] = 0;
                P[i] = 0; //periods
            }
            for (size_t i = 0; i < N; i++)
            {
                ind[i] = 0; // initialize the ind array

                for (size_t j = 0; j < 5; j++)
                    peak[i][j] = 0;
            }
        }

        compute_t2 = floor(t / compute_tau);
        //if ((t > (t_dark - 1000)) || (t < (t0 + 1000)))
        //  write_kuramoto(y, t, theta, light);

        dydt_in[0] = dydt_out[0];
        dydt_in[1] = dydt_out[1];
        t += h;
    }
*/
    gsl_odeiv_step_free(s);
    //}
    return (void *)0;
}

int main(int argc, char *argv[])
{
    auto start = high_resolution_clock::now();

    ofstream ofile;
    ofile.open("kuramoto.csv");
    ofile.close();

    //ofile.open("results.csv");
    //ofile.close();

    double K_f = stod(argv[1]);
    double daylength = stod(argv[2]);

    Enrange enrange =
        {
            daylength,
            K_f};

    int nthread = 30;
    //int same_num = 1;
    double values[nthread];
    pthread_t threads[nthread];

    pthread_rwlock_init(&rwlock, NULL);

    for (size_t i = 0; i < nthread; i++)
    {
        pthread_create(&(threads[i]), NULL, threadfunc, &(enrange));
    }

    for (size_t i = 0; i < nthread; i++)
        pthread_join(threads[i], NULL);

    pthread_rwlock_destroy(&rwlock);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << "Time taken by program: " << duration.count() << " seconds" << endl;

    return 0;
}
