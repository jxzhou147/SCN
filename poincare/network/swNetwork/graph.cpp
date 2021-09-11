#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <random>
#include <chrono>

#include "poincare.h"

using namespace std;
using namespace std::chrono;

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

void random_wiring(int i, double (*a)[N], double d)
{
    //srand((unsigned)time(NULL));

    int wiring_node = 0;
    double wiring_max = 0.0;
    double wiring_rand = 0;

    for (size_t j = 0; j < N; j++)
    {
        if ((a[i][j] == 0) & (dist(i, j) > d) & (i != j))
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

void random_wiring_random(int i, double (*a)[N], double d)
{
    //srand((unsigned)time(NULL));

    int wiring_node = 0;
    double wiring_max = 0.0;
    double wiring_rand = 0;

    for (size_t j = 0; j < N; j++)
    {
        if ((a[i][j] == 0) & (i != j))
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

void construct_connection(double (*a)[N], double p_link, double d, double p_d)
{
    //srand((unsigned)time(NULL));
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            a[i][j] = 0;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if ((dist(i, j) <= d) & (i != j) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_d))
            {
                a[i][j] = 1;

                if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link)
                {
                    a[i][j] = 0;
                    random_wiring(i, a, d);
                }
            }
        }
    }
}

void rewiring_random(double (*a)[N], double p_re, double d)
{
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if ((a[i][j] == 1) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_re))
            {
                a[i][j] = 0;
                random_wiring_random(i, a, d);
            }
        }
    }
}

void degree(const double (*a)[N], double deg[], double deg_in[], double deg_out[])
{
    for (size_t i = 0; i < N; i++)
    {
        deg_in[i] = 0;
        deg_out[i] = 0;
    }

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] == 1)
                deg_out[i]++;
            if (a[j][i] == 1)
                deg_in[i]++;
        }
        deg[i] = deg_in[i] + deg_out[i];
    }
}

void cluster(const double (*a)[N], const double deg[], double cluster_coe[])
{
    for (size_t i = 0; i < N; i++)
    {
        cluster_coe[i] = 0;
        if (deg[i] >= 2)
        {
            for (size_t j = 0; j < N; j++)
            {
                if ((a[i][j] == 1) || (a[j][i] == 1))
                {
                    for (size_t k = 0; k < N; k++)
                    {
                        if ((a[i][k] == 1) || (a[k][i] == 1))
                        {
                            if (a[j][k] == 1)
                                cluster_coe[i]++;
                        }
                    }
                }
            }
            cluster_coe[i] = cluster_coe[i] / (deg[i] * (deg[i] - 1));
        }
    }
}

double shortest_path(const double (*a)[N]) // Dijkstra method
{
    double sp = 0;
    double dist[N];
    bool s[N];
    const int INF = 1000000000;
    int ind = 0;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            dist[j] = INF;
            s[j] = false;
            if (a[i][j] == 1)
                dist[j] = 1;
        }

        dist[i] = 0;
        s[i] = true;

        for (size_t j = 0; j < N; j++)
        {
            int u = i;
            int mindist = INF;
            for (size_t k = 0; k < N; k++)
            {
                if (s[k] == false && dist[k] < mindist)
                {
                    u = k;
                    mindist = dist[k];
                }
            }
            s[u] = true;
            for (size_t k = 0; k < N; k++)
            {
                if (s[k] == false && a[u][k] == 1)
                {
                    if (dist[u] + 1 < dist[k])
                        dist[k] = dist[u] + 1;
                }
            }
        }

        for (size_t j = 0; j < N; j++)
        {
            if (dist[j] == INF)
            {
                dist[j] = 0;
                ind++;
            }
            sp += dist[j];
        }
    }

    sp = sp / ((double)N * (double)N - ind);
    return sp;
}

void write_graph(int n, const double (*a)[N], const double deg[], double cluster_coe[], double d, double p_d, double p_link, double p_re)
{
    ofstream ofile;
    ofile.open("graph.csv", ios::app);

    double ave_deg = 0;
    for (size_t i = 0; i < N; i++)
        ave_deg += deg[i];
    ave_deg /= (double)N;

    double ave_cluster = 0;
    for (size_t i = 0; i < N; i++)
        ave_cluster += cluster_coe[i];
    ave_cluster /= (double)N;

    ofile << n << '\t' << d << '\t' << p_d << '\t' << p_link << '\t' << p_re << '\t' << ave_deg << '\t' << ave_cluster << '\t' << shortest_path(a) << endl;
    ofile.close();
}

void write_degree_distribution(int n, const double deg[], double d, double p_d, double p_link, double p_re)
{
    ofstream ofile;
    ofile.open("degree.csv", ios::app);
    ofile << n << '\t' << d << '\t' << p_d << '\t' << p_link << '\t' << p_re << endl;
    for (size_t i = 0; i < N; i++)
        ofile << deg[i] << '\t';
    ofile << endl;
}

int main()
{
    double d = 1;
    double p_d = 0.6;
    double p_link = 0.1;
    double p_re = 0.1;
    int n0 = 100;

    double a[N][N];
    double deg[N];
    double deg_in[N];
    double deg_out[N];
    double cluster_coe[N];

    ofstream ofile;
    ofile.open("graph.csv");
    ofile.close();

    construct_connection(a, p_link, d, p_d);

    degree(a, deg, deg_in, deg_out);
    write_degree_distribution(0, deg, d, p_d, p_link, p_re);

    for (size_t n = 0; n < n0; n++)
    {
        //degree(a, deg, deg_in, deg_out);
        //cluster(a, deg, cluster_coe);
        //write_graph(n, a, deg, cluster_coe, d, p_d, p_link, p_re);
        rewiring_random(a, p_re, d);
    }

    degree(a, deg, deg_in, deg_out);
    write_degree_distribution(n0, deg, d, p_d, p_link, p_re);

    return 0;
}