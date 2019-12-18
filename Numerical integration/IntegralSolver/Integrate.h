#pragma once
#include <math.h>
#include <omp.h>
#include "Timer.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>

#define ZERO 1.0E-4
#define INF 5.0
#define MAXIT 10
#define EPS 3.0e-14
#define PI 3.14159265359

using namespace std;

double gaussLegendre(int N);
double gaussLaguerre(int N);
double gaussLaguerrePara(int N);
double* bruteMonteCarlo(int N);
double* impsampMonteCarlo(int N);
