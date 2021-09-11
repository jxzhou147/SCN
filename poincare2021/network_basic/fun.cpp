#include <fstream>
#include <iostream>
#include <math.h>
#include "parameter.h"
//#include "funs.h"

using namespace std;

void write_curve(const double y[], const double t, int light)
{
    ofstream ofile;
    ofile.open("curve.csv", ios::app);

    ofile << t << '\t';
    for (int i = 0; i < N; i++)
    {
        ofile << y[i] << '\t';
    }
    ofile << light << endl;

    ofile.close();
}

void write_connect(double **a)
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

void write_results(const double r, const double p_cut, double t_re)
{
    ofstream ofile;
    ofile.open("results.csv", ios::app);

    ofile << p_cut << '\t' << t_re << '\t' << r << endl;

    ofile.close();
}

void write_results_period(const double period_a, const double period_s, const double r, const double K, int *component)
{
    ofstream ofile;
    ofile.open("results.csv", ios::app);

    ofile << K << '\t' << r << '\t' << period_a << '\t' << period_s << '\t' << component[0] << '\t' << component[1] << endl;

    ofile.close();
}

void write_rt(const double t, const double p, const double r)
{
    ofstream ofile;
    ofile.open("rt.csv", ios::app);

    ofile << t << '\t' << p << '\t' << r << endl;

    ofile.close();
}

void write_er(const double tcycle, const double K_f, const double p, int ifEntrain)
{
    int ifEntrainment = 0;
    if (ifEntrain > 24)
        ifEntrainment = 1;

    ofstream ofile;
    ofile.open("results_er.csv", ios::app);
    ofile << tcycle << '\t' << K_f << '\t' << p << '\t' << ifEntrain << '\t' << ifEntrainment << endl;

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

void write_peaks_vd_p(const double peak_ave_v[], const double peak_ave_d[], const int peak_ave_num, const double p_vd, const double p_dv, const double daylength)
{
    ofstream ofile;
    ofile.open("peaks.csv", ios::app);

    ofile << "k_vd = " << p_vd << ", k_dv = " << p_dv << ", daylength = " << daylength << endl;

    for (size_t i = 0; i < peak_ave_num; i++)
        ofile << peak_ave_v[i] << '\t' << peak_ave_d[i] << endl;
    ofile.close();
}

void write_peaks_p(const double k_vd, const double peak_v[], const double peak_d[], double period_v[], double period_d[], double phase_gap[], const int peak_num)
{
    ofstream ofile;
    ofile.open("peaks.csv", ios::app);

    ofile << "k_vd = " << k_vd << endl;

    for (size_t i = 0; i < (peak_num - 1); i++)
    {
        period_v[i] = peak_v[i + 1] - peak_v[i];
        period_d[i] = peak_d[i + 1] - peak_d[i];
        phase_gap[i] = peak_v[i] - peak_d[i];
        ofile << peak_v[i] << '\t' << peak_d[i] << '\t' << period_v[i] << '\t' << period_d[i] << '\t' << phase_gap[i] << endl;
    }
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
        ofile << ((peak[i][4] - peak[i][1]) / 3) << endl;
    }

    ofile.close();
}

void write_net(double p_random, int connect_num, int part_num)
{
    ofstream ofile;
    ofile.open("net.csv", ios::app);

    ofile << p_random << ',' << connect_num << ',' << part_num << endl;

    ofile.close();
}

void find_peak(const double t, double (*peak)[5], double (*peak_now)[2], double (*peak_past)[2], double (*peak_old)[2], const double y[], int ind[])
{
    int ind_min = 5;
    for (size_t i_min = 0; i_min < (N + 1); i_min++)
    {
        if (ind_min > ind[i_min])
            ind_min = ind[i_min];
    }

    double average = 0;
    for (size_t i = 0; i < N; i++)
        average += y[i];
    average = average / N;

    if (ind_min < 5)
    {
        for (size_t i = 0; i < (N + 1); i++)
        {
            peak_old[i][0] = peak_past[i][0]; // [0] is x, [1] is t
            peak_old[i][1] = peak_past[i][1];
            peak_past[i][0] = peak_now[i][0];
            peak_past[i][1] = peak_now[i][1];
        }
        for (size_t i = 0; i < N; i++)
        {
            peak_now[i][0] = y[i];
            peak_now[i][1] = t;
        }
        peak_now[N][0] = average;
        peak_now[N][1] = t;

        for (size_t i_p = 0; i_p < (N + 1); i_p++)
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

double period_ave(const double (*peak)[5])
{
    double p = 0;

    for (int i = 0; i < N; i++)
        p += ((peak[i][4] - peak[i][1]) / 3);
    p = p / N;

    return p;
}

double period_std(const double (*peak)[5])
{
    double sum = 0;
    double per = 0;
    double ave = period_ave(peak);
    for (int i = 0; i < N; i++)
    {
        per = ((peak[i][4] - peak[i][1]) / 3);
        sum += pow(per - ave, 2);
    }
    sum = sum / N;

    return sqrt(sum);
}

double period(const double (*peak)[5])
{
    double p = 0;

    p = (peak[N][4] - peak[N][1]) / 3;

    return p;
}

double period_one(const double (*peak)[5])
{
    double p = 0;

    p = (peak[0][4] - peak[0][1]) / 3;

    return p;
}

double syn_degree(const double p, const double (*peak)[5])
{
    double r = 0;
    double sum1 = 0;
    double sum2 = 0;

    for (size_t j = 1; j < 5; j++)
    {
        sum1 = 0;
        sum2 = 0;
        for (size_t i = 0; i < N; i++)
        {
            sum1 += cos((peak[i][j] - peak[N][j]) * 2 * PI / p);
            sum2 += sin((peak[i][j] - peak[N][j]) * 2 * PI / p);
        }

        r += sqrt(pow(sum1, 2) + pow(sum2, 2)) / N;
    }

    r /= 4.0;

    return r;
}

int ifEntrainFun(const double (*peak)[5], const double tcycle)
{
    int entrainNum = 0;

    for (int i = 0; i < N; i++)
    {
        if (abs((peak[i][4] - peak[i][1]) / 3 - tcycle) < 0.1)
            entrainNum++;
    }

    if (entrainNum > (N * 0.9))
        return 1;
    return 0;
}