#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <random>
#include <pthread.h>

#include "new_funs.h"

using namespace std;

pthread_rwlock_t rwlock;
int tnum = 1;

void *threadfunc(void *arg)
{
    double p_random = *(double *)arg;
    int tid = tnum++;
    double p_cut = 0.1;
    double t_re = 1;

    double **a;
    a = new double *[N];
    for (int i = 0; i < N; i++)
    {
        a[i] = new double[N];
    }

    construct_connection_random(a, p_random, tid);

    int component[2] = {0};
    net_component(a, component);

    pthread_rwlock_wrlock(&rwlock);
    write_net(p_random, component[0], component[1]);
    pthread_rwlock_unlock(&rwlock);

    for (int i = 0; i < N; i++)
    {
        delete[] a[i];
    }
    delete[] a;

    return (void *)0;
}

int main(int argc, char *argv[])
{
    int nthread = 480;
    int same_num = 30;
    int iv;
    double values[nthread];
    pthread_t threads[nthread];

    pthread_rwlock_init(&rwlock, NULL);

    for (size_t i = 0; i < nthread; i++)
    {
        iv = i / same_num;
        //if (i <= 20)
        values[i] = floor((double)i / (double)same_num) * 0.001;
        //else
        //   values[i] = (floor((double)i / (double)same_num) - 20) * 0.1;
        //values[i] = pow(10, i) * 0.1;
        pthread_create(&(threads[i]), NULL, threadfunc, values + i);
    }

    for (size_t i = 0; i < nthread; i++)
        pthread_join(threads[i], NULL);

    pthread_rwlock_destroy(&rwlock);

    return 0;
}