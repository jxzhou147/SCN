#ifndef POINCARE_H
#define POINCARE_H

//parameters  in the poincare model
#define PI 3.1415926
#define N_rand 9999 //random number
const int N = 400;
const int N0 = 20;   //N0^2 = N
const int N_v = 140; //VIP neurons

//const double p_link = 0.05; //probability of long range connections
//const double d = 1.5; // the distance threshold of short-range links
//const double d = 2;

//const double t_re = 1;    //rewiring time gap
//const double p_re = 0.02; //rewiring probability

//const double K = 0.05;
const double omega_v = 2 * PI / 23.4;
const double omega_d = 2 * PI / 23;
//const double omega_v = 2 * PI / 23.6;
//const double omega_d = 2 * PI / 23.2;
const double Omega = 2 * PI / 24;
//const double K_f = 0.1; //K_f can be 0.5K, K, 2K
//const double K_f = 0.025;
//const double daylength = 12.0;

//const double omega_v = 2 * PI / 24.2;
//const double omega_d = 2 * PI / 23.8;

const double gam = 0.1;
const double a = 1;

#endif