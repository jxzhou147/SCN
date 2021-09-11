#include<iostream> 
#include<fstream>
#include<vector> 
#include<gsl/gsl_errno.h> 
#include<gsl/gsl_odeiv.h>
#include<math.h>
#include<random>

#include "kuramoto.h"
#include "fun.h"

using namespace std;


int func (double t, const double y[], double f[], void *params)
 {
	Kuramoto kura_func = *(Kuramoto *)params;

    for (size_t i = 0; i < N; i++)
	{
		kura_func.interaction[i] = 0;
		if (kura_func.e[i] != 0)
		{
			for (size_t j = 0; j < N; j++)
				kura_func.interaction[i] += kura_func.a[j][i] * sin(y[j] - y[i]) / kura_func.e[i];
		}
	}

    if (fmod(t, 2 * PI) < (daylength * 2 * PI / 24.0))
		kura_func.light = 1;
	else
		kura_func.light = 0;

    for (size_t i = 0; i < N_v; i++)
	{
		kura_func.theta[i] = fmod(y[i], 2 * PI);
		f[i] = omega_v * kura_func.eta[i] + K * kura_func.interaction[i] + K_f * kura_func.light * kura_func.phase(kura_func.theta[i]);
	}

	for (size_t i = N_v; i < N; i++)
	{
		f[i] = omega_d * kura_func.eta[i] + K * kura_func.interaction[i];
		kura_func.theta[i] = fmod(y[i], 2 * PI);
	}
 
	return GSL_SUCCESS; 
}

pair<int, int> Kuramoto::index_to_coordinate(int i)
{
	pair<int, int> coor(i / N0, i % N0);
	return coor;
}

double Kuramoto::dist(int i, int j)
{
	double xi = (double)(index_to_coordinate(i).first);
	double yi = (double)(index_to_coordinate(i).second);
	double xj = (double)(index_to_coordinate(j).first);
	double yj = (double)(index_to_coordinate(j).second);

	double dis = sqrt(pow((xi - xj), 2) + pow((yi - yj), 2));
	return dis;
}

void Kuramoto::construct_connection()
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (dist(i, j) < d)
			{
				a[i][j] = 1;
				if ((i < N_v) & (j > N_v))			//ventral to dorsal SCN's interaction is replusive (a_ij: i to j)
					a[i][j] = -1;
			}
			else
			{
				if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link)
				{
					a[i][j] = 1;
					if ((i < N_v) & (j > N_v))			
						a[i][j] = -1;
				}
				else
					a[i][j] = 0;
			}
		}
	}
}

void Kuramoto::rewiring()
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if ((a[i][j] != 0) & (dist(i, j) >= d))
			{
					if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_re)
						a[i][j] = 0;
			}

			if ((a[i][j] == 0) &((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_re))
			{
					if ((i < N_v) & (j > N_v))			
						a[i][j] = -1;
					else
						a[i][j] = 1;
			}
		}
	}
}

void Kuramoto::syn_degree()
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

	R_[2] += sqrt(pow(sum1, 2) + pow(sum2, 2) + pow(sum3, 2) + pow(sum4, 2));
}

double Kuramoto::delta_theta()
{
	double tv[5] = { 0, 0, 0, 0, 0 };
	double ts[5] = { 0, 0, 0, 0, 0 };
	double delta_the = 0;

	for (size_t i = 0; i < 5; i++)
	{
		for (size_t j = 0; j < N_v; j++)
			tv[i] += peak[j][i];
		for (size_t j = N_v; j < N; j++)
			ts[i] += peak[j][i];

		delta_the += abs((tv[i] / N_v - ts[i] / (N - N_v)));
	}

	delta_the = delta_the / 5;

	return delta_the;
}

void Kuramoto::period()
{
	for (size_t i = 0; i < N_v; i++)
		P[0] += ((peak[i][4] - peak[i][0]) / 4);
	P[0] /= N_v;
	for (size_t i = N_v; i < N; i++)
		P[1] += (peak[i][4] - peak[i][0]) / 4;
	P[1] /= (N - N_v);

	P[2] = (P[0] * N_v + P[1] * (N - N_v)) / N;
}

void Kuramoto::find_peak(const double t)
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
				//peak(i_p, ind[i_p]) = peak2[i_p];
				ind[i_p] = ind[i_p] + 1;
			}
		}

		for (size_t i = 0; i < N; i++)
			peak2[i] = theta[i];
	}
}

double Kuramoto::phase(double pha_pi)
{
	double pha = 2.0 * PI / 24.0;
	if ((pha > 3.0) & (pha < 9.0))
		return 0.0;
	else
	{
		if(pha > 9.0)
			return sin((pha - 9) * 2 * PI / 18.0);
		else
			return sin((pha + 24 - 9) * 2 * PI / 18.0);
	}
}

void Kuramoto::initialization()
{
	for (size_t i = 0; i < N; i++)
    {
		e[i] = 0;
		for (size_t j = 0; j < N; j++)
			e[i] += abs(a[j][i]);
	}

	normal_distribution<> norm{ 1, 0.05 };
	random_device rd;
	default_random_engine rng{ rd() };
	for (size_t i = 0; i < N; i++)
		eta[i] = norm(rng); 

	for (size_t i = 0; i < N; i++)
	{
		ind[i] = 0;			// initialize the ind array

		for (size_t j = 0; j < 5; j++)
			peak[i][j] = 0;
	}	
}

void Kuramoto::write_kuramoto(const double y[], const double t)
{
	ofstream ofile;
	ofile.open("kuramoto.csv", ios::app);
	ofile << t << '\t' << theta[0] << '\t' << theta[10] << '\t' << theta[40] << '\t' << theta[80] << '\t' << R_[1] << '\t' << R_[2] << endl;
	ofile.close();
}

void Kuramoto::write_peaks()
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

void Kuramoto::write_connect()
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

void Kuramoto::write_results()
{
	ofstream ofile;
	ofile.open("results.csv", ios::app);
	ofile << N << '\t' << t_re << '\t' << p_link << '\t' << p_re << '\t' << Omega << '\t' << daylength << '\t' << R[0] << '\t' << R[1] << '\t' << R[2] << '\t' << delta_theta() << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << endl;
	ofile.close();
}