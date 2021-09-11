#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <random>

#include "parameter.h"
#include "funs.h"

using namespace std;
using namespace std::chrono;

int main()
{ /*
    double v_sP[100], v_sB[100], v_mB[100];

    unsigned seed = system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    normal_distribution<double> norm(0, 1);

    for (int i = 0; i < 100; i++)
    {
        v_sP[i] = v_sP0 + v_sP0 * 0.06 * norm(gen);
        v_sB[i] = v_sB0 + v_sB0 * 0.01 * norm(gen);
        v_mB[i] = v_mB0 + v_mB0 * 0.01 * norm(gen);
        //cout << v_sP[i] << '\t' << v_sB[i] << '\t' << v_mB[i] << endl;
    }

    ofstream ofile;
    ofile.open("test.csv");

    for (int i = 0; i < 100; i++)
    {
        ofile << v_sP[i] << '\t' << v_sB[i] << '\t' << v_mB[i] << endl;
    }

    ofile.close();

    for (int i = 0; i < 100; i++)
    {
        
        //srand(i);
        double p_vd = 0.01;

        double **a_vip;
        double **a_gaba;
        a_vip = new double *[N];
        a_gaba = new double *[N];
        for (int i = 0; i < N; i++)
        {
            a_vip[i] = new double[N];
            a_gaba[i] = new double[N];
        }
        //double a_vip[N][N];
        //double a_gaba[N][N];
        int connect_num[2] = {0, 0}; // connect_num[0]: vip; connect_num[1]:gaba

        construct_connection_vd(a_vip, a_gaba, p_link, d, p_vd, 1);
        connection_number(a_vip, a_gaba, connect_num);

        cout << connect_num[0] << '/' << connect_num[1] << '\t';
        cout << (rand() / double(RAND_MAX)) << endl;

        for (int i = 0; i < N; i++)
        {
            delete[] a_vip[i];
            delete[] a_gaba[i];
        }
        delete[] a_vip;
        delete[] a_gaba;
        
        
    }
    */

    for (int i = 0; i < 100; i++)
        if (fmod(i, 24.0) >= (12.0))
        {
            cout << i << endl;
        }

    return 0;
}