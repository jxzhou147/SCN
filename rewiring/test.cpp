#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_TNUM 256
int tnum = 0, nthread;
pthread_cond_t cond;
pthread_mutex_t mtx;

void *threadfunc(void *arg)
{
    int tid, varg = *(int *)arg;

    pthread_mutex_lock(&mtx);
    tid = tnum++;
    printf("hello from thread_%d: value-of-arg=%d\n", tid, varg);
    if (nthread == tnum)
        pthread_cond_broadcast(&cond);
    else
        pthread_cond_wait(&cond, &mtx);
    pthread_mutex_unlock(&mtx);

    printf("bye from thread_%d\n", tid);

    return (void *)0;
}

int main(int argc, char **argv)
{
    int i, values[MAX_TNUM];
    pthread_t threads[MAX_TNUM];

    if (argc != 2)
    {
        printf("please enter the number of therads to be created\n");
        return 0;
    }
    else
    {
        sscanf(argv[1], "%d", &nthread);
        if (nthread > MAX_TNUM)
            nthread = MAX_TNUM;
    }

    pthread_mutex_init(&mtx, NULL);
    pthread_cond_init(&cond, NULL);

    printf("creating %d threads...\n", nthread);
    for (i = 0; i < nthread; i++)
    {
        values[i] = i;
        pthread_create(&(threads[i]), NULL, threadfunc, values + i);
    }
    printf("all creating requests have been submitted...\n");

    printf("waiting the created threads to terminate...\n");

    for (i = 0; i < nthread; i++)
        pthread_join(threads[i], NULL);

    printf("all created threads have terminated...\n");

    pthread_mutex_destroy(&mtx);
    pthread_cond_destroy(&cond);

    return 0;
}