#include <iostream>
#include <math.h>
#include <random>
#include <chrono>
#include "parameter.h"
#include <algorithm>
//#include "funs.h"

using namespace std;
using namespace std::chrono;

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

void dfs(int *vis, double **a, int v)
{
    vis[v] = 1;
    for (int i = 0; i < N; i++)
    {
        if (((a[v][i] == 1) || (a[i][v] == 1)) && (vis[i] == 0))
            dfs(vis, a, i);
    }
}

int ifConnected(double **a)
{
    int vis[N] = {0};
    dfs(vis, a, 0);

    int acc = 0;
    for (int i = 0; i < N; i++)
        acc += vis[i];

    return acc;
}

void net_component(double **a, int *component) // component[0]: number of parts; component[1]: node number of the bigest component
{
    component[0] = 0;
    component[1] = 0;
    int vis[N] = {0};
    int acc = 0;
    int acc_last = 0;
    int acc_once = 0;
    int dfs_ini = 0;

    while (acc != N)
    {
        for (int i = 0; i < N; i++)
        {
            if (vis[i] == 0)
            {
                dfs_ini = i;
                break;
            }
        }

        dfs(vis, a, dfs_ini);

        acc = 0;
        for (int i = 0; i < N; i++)
        {
            acc += vis[i];
        }

        acc_once = acc - acc_last;
        acc_last = acc;
        if (acc_once > component[1])
            component[1] = acc_once;
        component[0]++;
    }
}

int part_num(double **a)
{
    int num = 0;
    int vis[N] = {0};
    int acc = 0;
    int dfs_ini = 0;

    while (acc != N)
    {
        for (int i = 0; i < N; i++)
        {
            if (vis[i] == 0)
            {
                dfs_ini = i;
                break;
            }
        }
        dfs(vis, a, dfs_ini);
        acc = 0;
        for (int i = 0; i < N; i++)
        {
            acc += vis[i];
        }
        num++;
    }

    return num;
}

void random_wiring_vd(int i, double **a, double d, int Nh, int Ne, int N0)
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

void random_wiring_undirected(double **a, int tid)
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
    a[wiring_rand_j][wiring_rand_i] = 1;
}

void random_wiring_tri(double **a, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    int wiring = 1;
    int wiring_rand_i;
    int wiring_rand_j;

    while (wiring == 1)
    {
        wiring_rand_i = (rand() % N);
        wiring_rand_j = (rand() % N);

        if (fabs(a[wiring_rand_i][wiring_rand_j] - 0.2) < 1e-6)
            wiring = 0;
    }

    a[wiring_rand_i][wiring_rand_j] = 0.8;
}

void random_wiring_uniform(double **a, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    int wiring = 1;
    int wiring_rand_i;
    int wiring_rand_j;

    while (wiring == 1)
    {
        wiring_rand_i = (rand() % N);
        wiring_rand_j = (rand() % N);

        if (a[wiring_rand_i][wiring_rand_j] < 0.2)
            wiring = 0;
    }

    a[wiring_rand_i][wiring_rand_j] = (rand() / double(RAND_MAX)) / 5.0 + 0.8;
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

void construct_connection_random_undirected(double **a, double p_random, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;
        }

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i + 1; j < N; j++)
        {
            if ((rand() / double(RAND_MAX)) < p_random)
            {
                a[i][j] = 1;
                a[j][i] = 1;
            }
        }
    }
}
/*
void construct_connection_sw_triVal(double **a, double d, double p_random, double p_8, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;
            if (dist(i, j, N0) <= d) //!!! wrong! N0 is sqrt(N) , NOT N
                a[i][j] = 1;
        }
    }
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if ((a[i][j] == 1) && ((rand() / double(RAND_MAX)) < p_random))
            {
                random_wiring(a, tid);
                a[i][j] = 0;
            }
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] == 1)
            {
                if ((rand() / double(RAND_MAX)) < p_8)
                    a[i][j] = 0.8;
                else
                    a[i][j] = 0.2;
            }
        }
    }
}
*/

void BubbleSort(double *p, int length, int *ind_diff)
{
    for (int m = 0; m < length; m++)
    {
        ind_diff[m] = m;
    }

    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length - i - 1; j++)
        {
            if (p[j] > p[j + 1])
            {
                float temp = p[j];
                p[j] = p[j + 1];
                p[j + 1] = temp;

                int ind_temp = ind_diff[j];
                ind_diff[j] = ind_diff[j + 1];
                ind_diff[j + 1] = ind_temp;
            }
        }
    }
}

void construct_connection_sw_triVal(double **a, int m, double p_random, double p_8, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        double dis_i[N];
        int dis_index[N] = {0};
        for (size_t j = 0; j < N; j++)
        {
            dis_i[j] = dist(i, j, N0);
        }
        BubbleSort(dis_i, N, dis_index);

        for (size_t j = 0; j < m; j++)
        {
            a[i][dis_index[j]] = 1;
        }
        a[i][i] = 0;
    }

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if ((a[i][j] == 1) && ((rand() / double(RAND_MAX)) < p_random))
            {
                random_wiring(a, tid);
                a[i][j] = 0;
            }
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] == 1)
            {
                if ((rand() / double(RAND_MAX)) < p_8)
                    a[i][j] = 0.8;
                else
                    a[i][j] = 0.2;
            }
        }
    }
}

void construct_connection_sf_triVal(double **a, int m, double p_8, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            a[i][j] = 0;
        }
    }

    for (size_t i = 0; i < (m + 1); i++)
    {
        for (size_t j = (i + 1); j < (m + 1); j++)
        {
            a[i][j] = 1;
            a[j][i] = 1;
        }
    }

    int k[N] = {0};
    int k_tot = 0;
    double p_i[N] = {0};
    for (size_t i = (m + 1); i < N; i++)
    {
        k_tot = 0;
        for (size_t j = 0; j < i; j++)
        {
            k[j] = 0;
            p_i[j] = 0;
        }
        for (size_t j = 0; j < i; j++)
        {
            for (size_t jj = 0; jj < i; jj++)
            {
                if ((a[j][jj] == 1) | (a[jj][j] == 1))
                {
                    k[j]++;
                    k_tot++;
                }
            }
        }
        for (size_t j = 0; j < i; j++)
            p_i[j] = double(k[j]) / double(k_tot);

        for (size_t e_num = 0; e_num < m; e_num++)
        {
            double rand01 = (rand() / double(RAND_MAX));
            double sum_pi = 0;
            for (size_t j = 0; j < i; j++)
            {
                sum_pi += p_i[j];
                if ((sum_pi > rand01) & (a[j][i] == 0))
                {
                    a[j][i] = 1;
                    break;
                }
            }
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] == 1)
            {
                if ((rand() / double(RAND_MAX)) < p_8)
                    a[i][j] = 0.8;
                else
                    a[i][j] = 0.2;
            }
        }
    }
}

void construct_connection_triVal(double **a, double p_random, double p_8, int tid)
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
                if ((rand() / double(RAND_MAX)) < p_8)
                    a[i][j] = 0.8 * 3;
                else
                    a[i][j] = 0.2;
            }
        }
    }
}

void construct_connection_norm(double **a, double p_random, double p_sigma, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    unsigned seed = system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    normal_distribution<double> norm{0.7, p_sigma};

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
                a[i][j] = norm(gen);
                if (a[i][j] < 0)
                    a[i][j] = 0;
            }
        }
    }
}

void construct_connection_uniform(double **a, double p_random, int tid)
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
                a[i][j] = (rand() / double(RAND_MAX));
            }
        }
    }
}

void construct_connection_powerLaw(double **a, double p_random, int tid)
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
                a[i][j] = pow((rand() / double(RAND_MAX)), 2);
            }
        }
    }
}

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
            /*
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
*/
            //// rewiring for small world network
            if ((a[i][j] == 1) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link))
            {
                a[i][j] = 0;
                random_wiring_vd(i, a, d, 0, N_v, N1);
            }

            /*
            //// add remote connections
            if ((a[i][j] == 0) & ((rand() / double(RAND_MAX)) < p_link))
            {
                a[i][j] = 1;
            }
            */
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
            /*
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
*/
            //// rewiring for small world network
            if ((a[i][j] == 1) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link))
            {
                a[i][j] = 0;
                random_wiring_vd(i, a, d, N_v, N, N2);
            }
            /*
            //// add remote connections
            if ((a[i][j] == 0) & ((rand() / double(RAND_MAX)) < p_link))
            {
                a[i][j] = 1;
            }
            */
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

void dynamicNet(double **a, double p_cut)
{
    for (int i = 0; i < N_v; i++)
    {
        for (int j = 0; j < N_v; j++)
        {
            if ((a[i][j] == 1) && ((rand() / double(RAND_MAX)) < p_cut))
                a[i][j] = 0;
        }
    }
}

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

void dynamic_triVal(double **a, const double p_cut, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if ((fabs(a[i][j] - 0.8) < 1e-6) && ((rand() / double(RAND_MAX)) < p_cut))
            {
                random_wiring_tri(a, tid);
                a[i][j] = 0.2;
            }
        }
    }
}

void dynamic_norm(double **a, const double p_cut, double p_sigma, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    unsigned seed = system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    normal_distribution<double> norm{0.7, p_sigma};

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if ((fabs(a[i][j] - 0.0) > 1e-3) && ((rand() / double(RAND_MAX)) < p_cut))
            {
                a[i][j] = norm(gen);
                if (a[i][j] < 0)
                    a[i][j] = -1 * a[i][j];
            }
        }
    }
}

void dynamic_uniform(double **a, const double p_cut, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if ((a[i][j] > 0.8) && ((rand() / double(RAND_MAX)) < p_cut))
            {
                random_wiring_uniform(a, tid);
                a[i][j] = (rand() / double(RAND_MAX)) / 5.0;
            }
        }
    }
}

void dynamic_uniform_random(double **a, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] > 0)
            {
                a[i][j] = (rand() / double(RAND_MAX));
            }
        }
    }
}

void dynamic_uniform_powerLaw(double **a, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] > 0)
            {
                a[i][j] = pow((rand() / double(RAND_MAX)), 2);
            }
        }
    }
}

void dynamicWireAndCut_undirected(double **a, const double p_cut, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            if ((a[i][j] == 1) && ((rand() / double(RAND_MAX)) < p_cut))
            {
                random_wiring_undirected(a, tid);
                a[i][j] = 0;
                a[j][i] = 0;
            }
        }
    }
}

void dynamic_triVal_random(double **a, double p_8, int tid)
{
    srand(rand() * (unsigned int)(time(NULL)) * tid);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (a[i][j] != 0)
            {
                if ((rand() / double(RAND_MAX)) < p_8)
                    a[i][j] = 0.8 * 3;
                else
                    a[i][j] = 0.2;
            }
        }
    }
}