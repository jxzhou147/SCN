#include <fstream>
#include <iostream>
#include <math.h>
#include "parameter.h"
#include "funs.h"

using namespace std;

void write_MP(const double y[], const double t, int light)
{
    ofstream ofile;
    ofile.open("mp.csv", ios::app);

    ofile << t << '\t';
    for (int i = 0; i < N; i++)
    {
        ofile << y[i * Ns] << '\t';
    }
    ofile << light << endl;

    ofile.close();
}

void write_n1(const double y[], const double t, int light)
{
    ofstream ofile;
    ofile.open("n1.csv", ios::app);

    ofile << t << '\t';
    for (int i = 0; i < Ns; i++)
    {
        ofile << y[i] << '\t';
    }
    ofile << light << endl;

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

void write_fr(const double f_r[], const double t, int light)
{
    ofstream ofile;
    ofile.open("firingRate.csv", ios::app);

    ofile << t << '\t';

    for (int i = 0; i < N; i++)
        ofile << f_r[i] << '\t';

    ofile << light << endl;

    ofile.close();
}

void write_input(const double p_intput, const double t, int light)
{
    ofstream ofile;
    ofile.open("p_input.csv", ios::app);

    ofile << t << '\t' << light << '\t' << p_intput << endl;

    ofile.close();
}

void write_results(const double P[], const double P_var[], const double R[], const double p_vd, const int connect_num[], const double beta_light, double p_cut, double day)
{
    ofstream ofile;
    ofile.open("results.csv", ios::app);

    ofile << p_vd << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << '\t' << R[0] << '\t' << R[1] << '\t' << R[2] << '\t' << connect_num[0] << '\t' << connect_num[1] << endl;

    ofile.close();
}

void write_rt(const double t, const double P[], const double P_var[], const double R[])
{
    ofstream ofile;
    ofile.open("rt.csv", ios::app);

    ofile << t << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << '\t' << P_var[0] << '\t' << P_var[1] << '\t' << P_var[2] << '\t' << R[0] << '\t' << R[1] << '\t' << R[2] << endl;

    ofile.close();
}

void write_er(const double daylength, const double P[], int ifEntrain)
{
    int ifEntrainment = 0;
    if (ifEntrain > 24)
        ifEntrainment = 1;

    ofstream ofile;
    ofile.open("results_er.csv", ios::app);
    ofile << daylength << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << '\t' << ifEntrain << '\t' << ifEntrainment << endl;

    ofile.close();
}

void write_peaks(const double peak_ave[], const int peak_ave_num)
{
    ofstream ofile;
    ofile.open("peaks.csv", ios::app);

    for (size_t i = 0; i < peak_ave_num; i++)
        ofile << peak_ave[i] << endl;
    ofile.close();
}

void write_peaks_vd(const double peak_ave_v[], const double peak_ave_d[], const int peak_ave_num)
{
    ofstream ofile;
    ofile.open("peaks.csv", ios::app);

    for (size_t i = 0; i < peak_ave_num; i++)
        ofile << peak_ave_v[i] << '\t' << peak_ave_d[i] << endl;
    ofile.close();
}

void write_peaks_p(const double peak_ave[], const int peak_ave_num, const int tid, const double pulse)
{
    ofstream ofile;
    ofile.open("peaks.csv", ios::app);

    ofile << pulse << '\t' << tid << '\t';
    for (size_t i = 0; i < peak_ave_num; i++)
        ofile << peak_ave[i] << '\t';
    ofile << endl;

    ofile.close();
}

void write_period_N(const double (*peak)[5])
{
    ofstream ofile;
    ofile.open("period_N.csv", ios::app);

    for (size_t i = 0; i < N; i++)
    {
        ofile << ((peak[i][4] - peak[i][0]) / 4) << endl;
    }

    ofile.close();
}

void find_peak(const double t, double (*peak)[5], double (*peak_now)[2], double (*peak_past)[2], double (*peak_old)[2], const double y[], int ind[])
{
    int ind_min = 5;
    for (size_t i_min = 0; i_min < (N + 3); i_min++)
    {
        if (ind_min > ind[i_min])
            ind_min = ind[i_min];
    }

    double average_v = 0;
    for (size_t i = 0; i < N_v; i++)
        average_v += y[i * Ns];
    average_v = average_v / N_v;

    double average_d = 0;
    for (size_t i = N_v; i < N; i++)
        average_d += y[i * Ns];
    average_d = average_d / (N - N_v);

    double average = 0;
    for (size_t i = 0; i < N; i++)
        average += y[i * Ns];
    average = average / N;

    if (ind_min < 5)
    {
        for (size_t i = 0; i < (N + 3); i++)
        {
            peak_old[i][0] = peak_past[i][0]; // [0] is M_P, [1] is t
            peak_old[i][1] = peak_past[i][1];
            peak_past[i][0] = peak_now[i][0];
            peak_past[i][1] = peak_now[i][1];
        }
        for (size_t i = 0; i < N; i++)
        {
            peak_now[i][0] = y[i * Ns];
            peak_now[i][1] = t;
        }
        peak_now[N][0] = average_v;
        peak_now[N][1] = t;
        peak_now[N + 1][0] = average_d;
        peak_now[N + 1][1] = t;
        peak_now[N + 2][0] = average;
        peak_now[N + 2][1] = t;

        for (size_t i_p = 0; i_p < (N + 3); i_p++)
        {
            if ((peak_past[i_p][0] > 0) & (peak_past[i_p][0] > peak_old[i_p][0]) & (peak_past[i_p][0] > peak_now[i_p][0]) & (ind[i_p] < 5))
            {
                if (ind[i_p] == 0)
                {
                    peak[i_p][ind[i_p]] = peak_past[i_p][1];
                    ind[i_p] = ind[i_p] + 1;
                }
                else
                {
                    if ((peak_past[i_p][1] - peak[i_p][ind[i_p] - 1]) > 10)
                    {
                        peak[i_p][ind[i_p]] = peak_past[i_p][1];
                        ind[i_p] = ind[i_p] + 1;
                    }
                }
            }
        }
    }
}

void find_peak_ave(const double t, double peak_ave[], double peak_now[], double peak_past[], double peak_old[], const double y[], int ind_ave[])
{
    double average = 0;
    for (size_t i = 0; i < N; i++)
        average += y[i * Ns];
    average = average / N;

    peak_old[0] = peak_past[0]; // [0] is M_P, [1] is t
    peak_old[1] = peak_past[1];
    peak_past[0] = peak_now[0];
    peak_past[1] = peak_now[1];
    peak_now[0] = average;
    peak_now[1] = t;

    if ((peak_past[0] > peak_old[0]) & (peak_past[0] > peak_now[0]))
    {
        peak_ave[ind_ave[0]] = peak_past[1];
        ind_ave[0]++;
    }
}

void find_peak_ave_vd(const double t, double peak_ave_v[], double peak_now_v[], double peak_past_v[], double peak_old_v[], double peak_ave_d[], double peak_now_d[], double peak_past_d[], double peak_old_d[], const double y[], int ind_ave_v[], int ind_ave_d[])
{
    double average_v = 0;
    double average_d = 0;
    for (size_t i = 0; i < N_v; i++)
        average_v += y[i * Ns];
    average_v = average_v / N_v;
    for (size_t i = N_v; i < N; i++)
        average_d += y[i * Ns];
    average_d = average_d / (N - N_v);

    peak_old_v[0] = peak_past_v[0]; // [0] is M_P, [1] is t
    peak_old_v[1] = peak_past_v[1];
    peak_past_v[0] = peak_now_v[0];
    peak_past_v[1] = peak_now_v[1];
    peak_now_v[0] = average_v;
    peak_now_v[1] = t;
    peak_old_d[0] = peak_past_d[0]; // [0] is M_P, [1] is t
    peak_old_d[1] = peak_past_d[1];
    peak_past_d[0] = peak_now_d[0];
    peak_past_d[1] = peak_now_d[1];
    peak_now_d[0] = average_d;
    peak_now_d[1] = t;

    if ((peak_past_v[0] > peak_old_v[0]) & (peak_past_v[0] > peak_now_v[0]))
    {
        peak_ave_v[ind_ave_v[0]] = peak_past_v[1];
        ind_ave_v[0]++;
    }
    if ((peak_past_d[0] > peak_old_d[0]) & (peak_past_d[0] > peak_now_d[0]))
    {
        peak_ave_d[ind_ave_d[0]] = peak_past_d[1];
        ind_ave_d[0]++;
    }
}

void period(double P[], double P_var[], const double (*peak)[5])
{
    for (size_t i = 0; i < N_v; i++)
        P[0] += ((peak[i][4] - peak[i][1]) / 3);
    P[0] /= N_v;
    for (size_t i = 0; i < N_v; i++)
        P_var[0] += pow((peak[i][4] - peak[i][1]) / 3 - P[0], 2);
    P_var[0] = pow(P_var[0] / N_v, 0.5); // Attention! pow(x, 1 /2) is wrong!

    for (size_t i = N_v; i < N; i++)
        P[1] += (peak[i][4] - peak[i][1]) / 3;
    P[1] /= (N - N_v);
    for (size_t i = N_v; i < N; i++)
        P_var[1] += pow((peak[i][4] - peak[i][1]) / 3 - P[1], 2);
    P_var[1] = pow(P_var[1] / (N - N_v), 0.5);

    P[2] = (P[0] * N_v + P[1] * (N - N_v)) / N;
    for (size_t i = 0; i < N; i++)
        P_var[2] += pow((peak[i][4] - peak[i][1]) / 3 - P[2], 2);
    P_var[2] = pow(P_var[2] / N, 0.5);
}

void period_N1(double P[], double P_var[], const double (*peak)[5])
{
    for (size_t i = 0; i < N_v; i++)
        P[0] += (peak[i][2] - peak[i][1]);
    P[0] /= N_v;
    for (size_t i = 0; i < N_v; i++)
        P_var[0] += pow((peak[i][2] - peak[i][1]) - P[0], 2);
    P_var[0] = pow(P_var[0] / N_v, 0.5); // Attention! pow(x, 1 /2) is wrong!

    for (size_t i = N_v; i < N; i++)
        P[1] += (peak[i][2] - peak[i][1]);
    P[1] /= (N - N_v);
    for (size_t i = N_v; i < N; i++)
        P_var[1] += pow((peak[i][2] - peak[i][1]) - P[1], 2);
    P_var[1] = pow(P_var[1] / (N - N_v), 0.5);

    P[2] = (P[0] * N_v + P[1] * (N - N_v)) / N;
    for (size_t i = 0; i < N; i++)
        P_var[2] += pow((peak[i][2] - peak[i][1]) - P[2], 2);
    P_var[2] = pow(P_var[2] / N, 0.5);
}

void period_ave(double P[], double P_var[], const double (*peak)[5])
{
    P[0] = ((peak[N][4] - peak[N][1]) / 3);
    P[1] = (peak[N + 1][4] - peak[N + 1][1]) / 3;
    P[2] = (peak[N + 2][4] - peak[N + 2][1]) / 3;
}

void syn_degree(double R[], const double P[], const double (*peak)[5])
{
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;

    for (size_t i = 0; i < N_v; i++)
    {
        sum1 += cos((peak[i][2] - peak[N][2]) * 2 * PI / P[0]);
        sum2 += sin((peak[i][2] - peak[N][2]) * 2 * PI / P[0]);
    }

    R[0] += sqrt(pow(sum1, 2) + pow(sum2, 2)) / N_v;

    for (size_t i = N_v; i < N; i++)
    {
        sum3 += cos((peak[i][2] - peak[N + 1][2]) * 2 * PI / P[1]);
        sum4 += sin((peak[i][2] - peak[N + 1][2]) * 2 * PI / P[1]);
    }

    R[1] += sqrt(pow(sum3, 2) + pow(sum4, 2)) / (N - N_v);

    for (size_t i = 0; i < N; i++)
    {
        sum5 += cos((peak[i][2] - peak[N + 2][2]) * 2 * PI / P[2]);
        sum6 += sin((peak[i][2] - peak[N + 2][2]) * 2 * PI / P[2]);
    }

    R[2] += sqrt(pow(sum5, 2) + pow(sum6, 2)) / N;
}

void syn_degree_24(double R[], const double P[], const double (*peak)[5])
{
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;

    for (size_t i = 0; i < N_v; i++)
    {
        sum1 += cos((peak[i][2] - peak[N][2]) * 2 * PI / 24.0);
        sum2 += sin((peak[i][2] - peak[N][2]) * 2 * PI / 24.0);
    }

    R[0] += sqrt(pow(sum1, 2) + pow(sum2, 2)) / N_v;

    for (size_t i = N_v; i < N; i++)
    {
        sum3 += cos((peak[i][2] - peak[N + 1][2]) * 2 * PI / 24.0);
        sum4 += sin((peak[i][2] - peak[N + 1][2]) * 2 * PI / 24.0);
    }

    R[1] += sqrt(pow(sum3, 2) + pow(sum4, 2)) / (N - N_v);

    for (size_t i = 0; i < N; i++)
    {
        sum5 += cos((peak[i][2] - peak[N + 2][2]) * 2 * PI / 24.0);
        sum6 += sin((peak[i][2] - peak[N + 2][2]) * 2 * PI / 24.0);
    }

    R[2] += sqrt(pow(sum5, 2) + pow(sum6, 2)) / N;
}

void syn_degree_N1(double R[], const double P[], const double (*peak)[5])
{
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;

    for (size_t i = 0; i < N_v; i++)
    {
        sum1 += cos((peak[i][0] - peak[N][0]) * 2 * PI / P[0]);
        sum2 += sin((peak[i][0] - peak[N][0]) * 2 * PI / P[0]);
    }

    R[0] += sqrt(pow(sum1, 2) + pow(sum2, 2)) / N_v;

    for (size_t i = N_v; i < N; i++)
    {
        sum3 += cos((peak[i][0] - peak[N + 1][0]) * 2 * PI / P[1]);
        sum4 += sin((peak[i][0] - peak[N + 1][0]) * 2 * PI / P[1]);
    }

    R[1] += sqrt(pow(sum3, 2) + pow(sum4, 2)) / (N - N_v);

    for (size_t i = 0; i < N; i++)
    {
        sum5 += cos((peak[i][0] - peak[N + 2][0]) * 2 * PI / P[2]);
        sum6 += sin((peak[i][0] - peak[N + 2][0]) * 2 * PI / P[2]);
    }

    R[2] += sqrt(pow(sum5, 2) + pow(sum6, 2)) / N;
}

double real_root_no(double x, double n) // If n is like 0.5, 1.5, when x < 0 ,there is no real root
{
    if (x >= 0)
        return pow(x, n);
    else
        return 0.0;
}

double real_root_transfer(double x, double n) // If n is like 0.2, 2.2, when x < 0 ,transfer to the real root
{
    if (x >= 0)
        return pow(x, n);
    else
        return -1.0 * pow(-1 * x, n);
}

double real_root_zero(double x, double n) // If n < 0, when x = 0 , cout a problem
{
    if (x == 0)
    {
        std::cout << "wrong!" << std::endl;
        return pow(0.1, n);
    }
    else
        return pow(x, n);
}