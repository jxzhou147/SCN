#ifndef KURAMOTO_H
#define KURAMOTO_H

//parameters  in the kuramoto model
#define PI 3.1415926
#define N_rand  99 //random number
const int N = 100;
const int N0 = 10;		//N0^2 = N
const int N_v = 20;		//VIP neurons

const double p_link = 0.05;		//probability of long range connections
const double d = 1.1;			// the distance threshold of short-range links

const double t_re = 5;		//rewiring time gap
const double p_re = 0.02;		//rewiring probability

const double K = 0.05;
const double omega_v = 2 * PI / 23.6;
const double omega_d = 2 * PI / 23.2;
const double Omega = 2 * PI / 24;
const double K_f = 0.05;		//K_f can be 0.5K, K, 2K
const double daylength = 0.0;

#endif