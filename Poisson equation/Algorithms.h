#pragma once
#include "armadillo"

inline double f(double x) { return 100.0 * exp(-10.0 * x); }
inline double analytical(double x) { return 1 - (1 - exp(-10)) * x - exp(-10 * x); }
double* algo_9n(int exponent);
double* algo_4n(int exponent);
arma::vec arma_lu(int exponent);
double errorAnalysis(double* numeric, double* analytic, int exponent);
