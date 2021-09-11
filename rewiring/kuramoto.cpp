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

#include "kuramoto.h"

using namespace std;
using namespace std::chrono;

double a[N][N];
double peak[N][5];

int ind[N];

double e[N];
double interaction[N];
double eta[N];
double theta[N];
double peak1[N];
double peak2[N];

double R_[3] = {0, 0, 0};
double R[3] = {0, 0, 0};
double P[3] = {0, 0, 0}; //periods
int times = 0;

int ton;
int toff;

int light;
double tau_off = 3.0;
double tau_on = 3.0;
double daylength = 3;

pair<int, int>
index_to_coordinate(int i)
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

void random_wiring(int i)
{
	srand((unsigned)time(NULL));

	int wiring_node = 0;
	double wiring_max = 0.0;
	double wiring_rand = 0;

	for (size_t j = 0; j < N; j++)
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

	//if ((i < N_v) & (wiring_node > N_v)) //ventral to dorsal SCN's interaction is replusive (a_ij: i to j)
	//	a[i][wiring_node] = -1;
}

void weighted_wiring(int i, const double y[])
{
	srand((unsigned)time(NULL));

	double dif;
	int wiring_node = 0;
	double wiring_max = 0.0;
	double wiring_rand = 0;

	for (size_t j = 0; j < N; j++)
	{
		if ((a[i][j] == 0) & (i != j))
		{
			//dif = exp(fmod(abs(y[i] - y[j]), 2 * PI) / tau_on);
			dif = fmod(abs(y[i] - y[j]), 2 * PI) / tau_on + 1;
			wiring_rand = (rand() % (N_rand + 1) / (float)(N_rand + 1)) / dif;
			if (wiring_rand > wiring_max)
			{
				wiring_node = j;
				wiring_max = wiring_rand;
			}
		}
	}

	a[i][wiring_node] = 1;

	//if ((i < N_v) & (wiring_node > N_v)) //ventral to dorsal SCN's interaction is replusive (a_ij: i to j)
	//	a[i][wiring_node] = -1;
}

void construct_connection()
{
	srand((unsigned)time(NULL));

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			a[i][j] = 0;

			if ((dist(i, j) < d) & (i != j))
			{
				a[i][j] = 1;

				if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link)
				{
					//a[i][j] = 0;
					random_wiring(i);
				}
			}

			//else if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_link)
			//	a[i][j] = 1;

			//if ((i < N_v) & (j > N_v)) //ventral to dorsal SCN's interaction is replusive (a_ij: i to j)
			//	a[i][j] = -1 * a[i][j];
		}
	}
}

void rewiring(const double y[])
{
	srand((unsigned)time(NULL));

	double dif;
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			//dif = exp(fmod(abs(y[i] - y[j]), 2 * PI) / tau_off); //weights of rewiring and cutting
			dif = fmod(abs(y[i] - y[j]), 2 * PI) / tau_on;

			if ((a[i][j] != 0) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < (p_re * dif)))
			{
				a[i][j] = 0;
				//weighted_wiring(i, y);
				random_wiring(i);
			}

			//			if ((a[i][j] == 0) & ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < (p_re / dif)))
			//				a[i][j] = 1;
			/*
			if ((rand() % (N_rand + 1) / (float)(N_rand + 1)) < p_re)
			{
				if (a[i][j] == 1)
				{
					a[i][j] = 0;
					random_wiring(i);
				}
			}*/
		}
	}
}

int connection_number()
{
	int num = 0;
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (a[i][j] != 0)
				num++;
		}
	}

	return num;
}

void syn_degree()
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

double delta_theta()
{
	double tv[5] = {0, 0, 0, 0, 0};
	double ts[5] = {0, 0, 0, 0, 0};
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

void period()
{
	for (size_t i = 0; i < N_v; i++)
		P[0] += ((peak[i][4] - peak[i][0]) / 4);
	P[0] /= N_v;

	for (size_t i = N_v; i < N; i++)
		P[1] += (peak[i][4] - peak[i][0]) / 4;
	P[1] /= (N - N_v);

	P[2] = (P[0] * N_v + P[1] * (N - N_v)) / N;
}

void find_peak(const double t)
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

int kuramoto(double t, const double y[], double f[], void *params)
{
	for (size_t i = 0; i < N; i++)
	{
		interaction[i] = 0;
		if (e[i] != 0)
		{
			for (size_t j = 0; j < N; j++)
				interaction[i] += a[j][i] * sin(y[j] - y[i]) / e[i];
		}
	}

	if (fmod(t, 24.0) < daylength)
		light = 1;
	else
		light = 0;

	for (size_t i = 0; i < N_v; i++)
	{
		//dxdt[i] = omega * eta[i] + K * interaction[i] + K_f * sin(Omega * t - x[i]);
		f[i] = omega_v * eta[i] + K * interaction[i] + K_f * light * phase(theta[i]);
		theta[i] = fmod(y[i], 2 * PI);
	}

	for (size_t i = N_v; i < N; i++)
	{
		f[i] = omega_d * eta[i] + K * interaction[i];
		theta[i] = fmod(y[i], 2 * PI);
	}

	return GSL_SUCCESS;
}

void write_kuramoto(const double y[], const double t)
{
	double ave_theta = 0;
	for (size_t i = 0; i < N; i++)
		ave_theta += theta[i];
	ave_theta /= (double)N;

	double phase_sun = fmod(Omega * t, 2 * PI);

	ofstream ofile;
	ofile.open("kuramoto.csv", ios::app);
	ofile << t << '\t' << ave_theta << '\t' << phase_sun << '\t' << theta[0] << '\t' << theta[10] << '\t' << theta[40] << '\t' << theta[80] << endl;
	ofile.close();
}

void write_peaks()
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

void write_connect()
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

void write_results()
{
	ofstream ofile;
	ofile.open("results.csv", ios::app);
	ofile << N << '\t' << connection_number() << '\t' << t_re << '\t' << d << '\t' << p_link << '\t' << p_re << '\t' << Omega << '\t' << daylength << '\t' << K_f << '\t' << R[0] << '\t' << R[1] << '\t' << R[2] << '\t' << delta_theta() << '\t' << P[0] << '\t' << P[1] << '\t' << P[2] << endl;
	//ofile << daylength << '\t' << connection_number() << '\t';
	ofile.close();
}

int main()
{
	auto start = high_resolution_clock::now();

	ofstream ofile;

	ofile.open("kuramoto.csv");
	ofile.close();

	ofile.open("peaks.csv");
	ofile.close();

	ofile.open("connection.csv");
	ofile.close();

	//for (int i = 0; i < 25; i++)
	//{
	//	daylength = double(i);

	srand((unsigned)time(NULL));

	for (size_t i = 0; i < 10; i++)
	{
		cout << (rand() % (N_rand + 1) / (float)(N_rand + 1)) << endl;
	}

	times = 0;

	for (size_t i = 0; i < 3; i++)
	{
		R_[i] = 0;
		R[i] = 0;
		P[i] = 0; //periods
	}

	construct_connection();

	for (size_t i = 0; i < N; i++)
	{
		e[i] = 0;
		for (size_t j = 0; j < N; j++)
			e[i] += abs(a[j][i]);
	}

	normal_distribution<> norm{1, 0.05};
	random_device rd;
	default_random_engine rng{rd()};
	for (size_t i = 0; i < N; i++)
		eta[i] = norm(rng);

	const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, N);

	gsl_odeiv_system sys = {kuramoto, NULL, N, NULL};

	double t = 0.0, t0 = 4000.0, t1 = 600, t2 = 900, t3 = 1200; // t0 is the time of light off
	double h = 0.1;
	double y_err[N];
	double dydt_in[N], dydt_out[N];
	double y[N];

	for (size_t i = 0; i < N; i++)
	{
		y[i] = 2 * PI * (rand() % (N_rand + 1) / (float)(N_rand + 1));
		theta[i] = fmod(y[i], 2 * PI);
	}

	for (size_t i = 0; i < N; i++)
	{
		ind[i] = 0; // initialize the ind array

		for (size_t j = 0; j < 5; j++)
			peak[i][j] = 0;
	}

	GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
	while (t < t0)
	{
		int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
		if (status != GSL_SUCCESS)
			break;

		if (t > (t0 - 150))
		{
			if (times < 100)
			{
				times++;
				syn_degree();
			}
		}

		///////// rewiring
		if (t > 100)
		{
			ton = floor(t / t_re);
			if (ton != toff)
			{
				rewiring(y);

				for (size_t i = 0; i < N; i++)
				{
					e[i] = 0;
					for (size_t j = 0; j < N; j++)
						e[i] += abs(a[j][i]);
				}
			}
			toff = floor(t / t_re);
		}

		if (t > (t0 - 200))
			find_peak(t);

		write_kuramoto(y, t);

		for (size_t i = 0; i < N; i++)
			dydt_in[i] = dydt_out[i];
		t += h;
	}

	R[0] = R_[0] / (N_v * times);
	R[1] = R_[1] / ((N - N_v) * times);
	R[2] = R_[2] / (N * times);

	period();

	write_peaks();
	//write_connect();
	write_results();
	/*
	daylength = 18;

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

	while (t < t1)
	{
		int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
		if (status != GSL_SUCCESS)
			break;

		if (t > (t1 - 150))
		{
			if (times < 100)
			{
				times++;
				syn_degree();
			}
		}

		///////// rewiring
		if (t > 400)
		{
			ton = floor(t / t_re);
			if (ton != toff)
			{
				rewiring(y);

				for (size_t i = 0; i < N; i++)
				{
					e[i] = 0;
					for (size_t j = 0; j < N; j++)
						e[i] += abs(a[j][i]);
				}
			}
			toff = floor(t / t_re);
		}

		if (t > (t1 - 200))
			find_peak(t);

		//write_kuramoto(y, t);

		dydt_in[0] = dydt_out[0];
		dydt_in[1] = dydt_out[1];
		t += h;
	}

	R[0] = R_[0] / (N_v * times);
	R[1] = R_[1] / ((N - N_v) * times);
	R[2] = R_[2] / (N * times);

	period();

	//write_peaks();
	//write_connect();
	write_results();

	daylength = 12;

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

	while (t < t2)
	{
		int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
		if (status != GSL_SUCCESS)
			break;

		if (t > (t2 - 150))
		{
			if (times < 100)
			{
				times++;
				syn_degree();
			}
		}

		///////// rewiring
		if (t > 400)
		{
			ton = floor(t / t_re);
			if (ton != toff)
			{
				rewiring(y);

				for (size_t i = 0; i < N; i++)
				{
					e[i] = 0;
					for (size_t j = 0; j < N; j++)
						e[i] += abs(a[j][i]);
				}
			}
			toff = floor(t / t_re);
		}

		if (t > (t2 - 200))
			find_peak(t);

		//write_kuramoto(y, t);

		dydt_in[0] = dydt_out[0];
		dydt_in[1] = dydt_out[1];
		t += h;
	}

	R[0] = R_[0] / (N_v * times);
	R[1] = R_[1] / ((N - N_v) * times);
	R[2] = R_[2] / (N * times);

	period();

	//write_peaks();
	//write_connect();
	write_results();

	daylength = 6;

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

	while (t < t3)
	{
		int status = gsl_odeiv_step_apply(s, t, h, y, y_err, dydt_in, dydt_out, &sys);
		if (status != GSL_SUCCESS)
			break;

		if (t > (t3 - 150))
		{
			if (times < 100)
			{
				times++;
				syn_degree();
			}
		}

		///////// rewiring
		if (t > 400)
		{
			ton = floor(t / t_re);
			if (ton != toff)
			{
				rewiring(y);

				for (size_t i = 0; i < N; i++)
				{
					e[i] = 0;
					for (size_t j = 0; j < N; j++)
						e[i] += abs(a[j][i]);
				}
			}
			toff = floor(t / t_re);
		}

		if (t > (t3 - 200))
			find_peak(t);

		//write_kuramoto(y, t);

		dydt_in[0] = dydt_out[0];
		dydt_in[1] = dydt_out[1];
		t += h;
	}

	R[0] = R_[0] / (N_v * times);
	R[1] = R_[1] / ((N - N_v) * times);
	R[2] = R_[2] / (N * times);

	period();

	//write_peaks();
	//write_connect();
	write_results();
*/
	gsl_odeiv_step_free(s);
	//}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cout << "Time taken by program: " << duration.count() << " seconds" << endl;

	return 0;
}
