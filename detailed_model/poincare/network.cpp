#include <iostream>
#include <math.h>
#include <random>
#include <chrono>
#include "parameter.h"
//#include "funs.h"

using namespace std;

pair<int, int> index_to_coordinate(int i, int N0)
{
    pair<int, int> coor(i / N0, i % N0);
    return coor;
}

double dist(int i, int j, int N0)
{
    double xi = (double)(index_to_coordinate(i, N0).first);
    double yi = (double)(index_to_coordinate(i, N0).second);
    double xj = (double)(index_to_coordinate(j, N0).first);
    double yj = (double)(index_to_coordinate(j, N0).second);

    double dis = sqrt(pow((xi - xj), 2) + pow((yi - yj), 2));
    return dis;
}

/*
void random_wiring(int i, double (*a)[N], double d, int Nh, int Ne, int N0)
{
    int wiring_node = 0;
    double wiring_max = 0.0;
    double wiring_rand = 0;

    for (size_t j = Nh; j < Ne; j++)
    {
        if ((a[i][j] == 0) & (dist(i, j, N0) > d) & (i != j))
        {
            wiring_rand = (rand() / double(RAND_MAX));
            if (wiring_rand > wiring_max)
            {
                wiring_node = j;
                wiring_max = wiring_rand;
            }
        }
    }

    a[i][wiring_node] = 1;
}
*/

void random_wiring(double **a, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    int wiring = 1;
    int wiring_rand_i;
    int wiring_rand_j;

    while (wiring == 1)
    {
        wiring_rand_i = (rand() % N);
        wiring_rand_j = (rand() % N);

        if (a[wiring_rand_i][wiring_rand_j] == 0)
            wiring = 0;
    }

    a[wiring_rand_i][wiring_rand_j] = 1;
}

void construct_connection_random(double **a, double p_random, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;
        }

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if ((rand() / double(RAND_MAX)) < p_random)
            {
                a[i][j] = 1;
            }
        }
    }
}

#if 0
void construct_connection_vd(double **a, double p_link, double d, double p_vd, double p_dv, int tid)
{
    //srand((unsigned int)(time(NULL)));
    srand(rand() * (unsigned int)(time(NULL)) * tid);
    //// Initialization
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;
        }

    //// Ventral to Ventral
    for (size_t i = 0; i < N_v; i++)
    {
        for (size_t j = 0; j < N_v; j++)
        {
            //// Near connections
            if ((dist(i, j, N1) <= d) & (i != j))
                a[i][j] = 1;

            //// periodic boundary
            if ((i % N1) == 0)
            {
                a[i][i + N1 - 1] = 1;
                a[i + N1 - 1][i] = 1;
            }
            if (i < N1)
            {
                a[i][N_v + i - N1] = 1;
                a[N_v + i - N1][i] = 1;
            }
            /*
            //// rewiring for small world network
            if ((a_gaba[i][j] == 1) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link))
            {
                a_gaba[i][j] = 0;
                random_wiring(i, a_gaba, d, 0, N_v, N1);
            }
*/
            //// add remote connections
            if ((a[i][j] == 0) & ((rand() / double(RAND_MAX)) < p_link))
            {
                a[i][j] = 1;
            }
        }
    }

    //// Dorsal to Dorsal
    for (size_t i = N_v; i < N; i++)
    {
        for (size_t j = N_v; j < N; j++)
        {
            //// Near connections
            if ((dist(i - N_v, j - N_v, N2) <= d) & (i != j))
                a[i][j] = 1;

            //// periodic boundary
            if (((i - N_v) % N2) == 0)
            {
                a[i][i + N2 - 1] = 1;
                a[i + N2 - 1][i] = 1;
            }
            if ((i - N_v) < N2)
            {
                a[i][N - N_v + i - N_v - N2] = 1;
                a[N - N_v + i - N_v - N2][i] = 1;
            }

            //// add remote connections
            if ((a[i][j] == 0) & ((rand() / double(RAND_MAX)) < p_link))
            {
                a[i][j] = 1;
            }
        }
    }

    //// Ventral to Dorsal
    for (size_t i = 0; i < N_v; i++)
    {
        for (size_t j = N_v; j < N; j++)
        {
            if ((rand() / double(RAND_MAX)) < p_vd)
            {
                a[i][j] = 1;
            }
        }
    }

    //// Dorsal to Ventral
    for (size_t i = N_v; i < N; i++)
    {
        for (size_t j = 0; j < N_v; j++)
        {
            if ((rand() / double(RAND_MAX)) < p_dv)
            {
                a[i][j] = 1;
            }
        }
    }
}
#endif

int connection_number(double **a)
{
    int connect_num = 0;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] == 1)
                connect_num++;
        }
    }

    return connect_num;
}

#if 0
void dynamicNet(double **a, double p_cut)
{
    for (int i = 0; i < N_v; i++)
    {
        for (int j = N_v; j < N; j++)
        {
            if ((a[i][j] == 1) && ((rand() / double(RAND_MAX)) < p_cut))
                a[i][j] = 0;
        }
    }
}
#endif

// Accumulation and consumption: only cutting or wiring at one time
void dynamicAccAndCon(double **a_vip, double **a_gaba, double **a_vip_static, double **a_gaba_static, const double t, double delay, double p_cut, double day)
{
}

// Preferential wiring and cutting under entrainment: wiring number = cutting number * h(t)
void dynamicWireAndCut(double **a_vip, double **a_gaba, double **a_vip_static, double **a_gaba_static, const double t, double delay, double p_cut, double day)
{
}

//  Random wiring and cutting without entrainment: cutting number = wiring number
void dynamicWireAndCut(double **a, const double p_cut, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if ((a[i][j] == 1) && ((rand() / double(RAND_MAX)) < p_cut))
            {
                random_wiring(a, tid);
                a[i][j] = 0;
            }
        }
    }
}