#include"Algorithms.h"
#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include"time.h"
#include"armadillo"

using namespace arma;
using namespace std;


double* algo_9n(int exponent) {
	// Algorithm for solving Dirichlet boundary value problem with
	// a total of 9 FLOPS.  N = 10^exponent steps.
	// Writes solution to file u_data_slowN.dat and returns pointer to solution.
		
	int N = pow(10, exponent);
		
	// Opening file for result u
	ofstream u_data;
	stringstream filename;
	filename << "u_data_slow" << N << "n.dat";
	u_data.open(filename.str());

	// Initializing stepsize as function of integration points
	double h = 1.0 / N;

	// Allocating memory for arrays
	double* u, * d, * d_tilde, * g, * g_tilde, * a, * b, * x;
	u = new double[N + 1];
	d = new double[N + 1];
	d_tilde = new double[N + 1];
	g = new double[N + 1];
	g_tilde = new double[N + 1];
	a = new double[N+1];
	b = new double[N+1];
	x = new double[N + 1];

	// Boundary conditions
	u[0] = u[N] = 0.0;

	// Discretizing x and initializing RHS g(x_i) and diagonal-arrays d, a, b
	for (int i = 0; i <= N ; i++) {
		x[i] = i * h;
		g[i] = pow(h, 2) * f(x[i]);
		d[i] = 2.0;
		a[i] = -1.0;
		b[i] = -1.0;
	}

	// Starting timer
	clock_t start, finish;
	start = clock();
	
	// Forward substitution (6 FLOPS, 8 mermory reads, 2 memory writes)
	d_tilde[0] = d[0];
	g_tilde[0] = g[0];
		
	for (int i = 1; i < N + 1; i++) {
		d_tilde[i] = d[i] - a[i] * b[i - 1] / d_tilde[i - 1];
		g_tilde[i] = g[i] - a[i] * g_tilde[i - 1] / d_tilde[i - 1];
	}

	// Backward substitution (3 FLOPS, 4 memory reads, 1 memory write)
	u_data << u[N] << "\n";
	u[N] = g_tilde[N] / d_tilde[N];
	for (int i = N - 1; i > 0; i--) {
		u[i] = (g_tilde[i] - b[i] * u[i+1]) / d_tilde[i];
		u_data << u[i] << '\n';
	}
		
	// Stop timer
	finish = clock();
	cout << "Total time spent on algoritghm with n = " << N << " was " << ((double)(finish - start) / CLOCKS_PER_SEC) * 1000 << "ms" << endl;
	u_data << u[0];
	u_data.close();

	// Freeing up memory
	delete[] d; delete[] d_tilde; delete[] g; delete[] g_tilde;
	delete[] a; delete[] b; delete[] x;

	return u;
}


double* algo_4n(int exponent) {
	// Algorithm for solving Dirichlet boundary value problem with
	// a total of 4 FLOPS.  N = 10^exponent steps.
	// Writes solution to file u_data_fastN.dat and returns pointer to solution.

	int N = pow(10, exponent);

	// Opening file for result u
	ofstream u_data;
	stringstream filename;
	filename << "u_data_fast" << N << "n.dat";
	u_data.open(filename.str());

	// Initializing stepsize as a function of n
	double h = 1.0 / N;

	// Allocating memory for arrays
	double* u, * d_tilde, * g, * g_tilde, * x;
	u = new double[N + 1];
	g = new double[N + 1];
	g_tilde = new double[N + 1];
	d_tilde = new double[N + 1];
	x = new double[N + 1];

	// Boundary conditions
	u[0] = u[N] = 0.0;

	// Discretizing x and initializing RHS g(x_i) and d
	// Contains precalculations to minimize FLOPS
	for (int i = 0; i <= N; i++) {
		x[i] = i * h;
		g[i] = pow(h, 2) * f(x[i]);
	}

	g_tilde[0] = g[0];
	d_tilde[0] = 2;
	for (int i = 1; i < N + 1; i++) {
		d_tilde[i] = 2.0- 1.0 / d_tilde[i-1];
		g_tilde[i] = g[i] + g_tilde[i - 1] /d_tilde[i-1];
	}

	// Starting timer
	clock_t start, finish;
	start = clock();

	// Backward substitution (4 FLOPS, 2 memory reads, 1 memory write)
	u_data << u[N] << "\n";
	u[N] = g_tilde[N] / d_tilde[N];
	for (int i = N - 1; i > 0; i--) {
		u[i] = g_tilde[i] /d_tilde[i] + u[i+1] / d_tilde[i];
		u_data << u[i] << "\n";
	}

	// Stop timer
	finish = clock();
	cout << "Total time spent on algoritghm with n = " << N << " was " << ((double)(finish - start) / CLOCKS_PER_SEC) * 1000 << "ms" << endl;
	u_data << u[0];
	u_data.close();

	// Freeing up memory
	delete[] g; delete[] g_tilde; delete[] x;

	return u;
}

vec arma_lu(int exponent) {
	// Algorithm for solving Dirichlet boundary value problem with
	// Armadillos LU decomposition functionality. N = 10^exponent steps.
	// Writes solution to file u_data_armaN.dat and returns vector solution.

	int N = pow(10, exponent);

	// Opening file for result u
	ofstream u_data;
	stringstream filename;
	filename << "u_data_arma" << N << "n.dat";
	u_data.open(filename.str());

	// Initializing stepsize as a function of n
	double h = 1.0 / N;

	// Discretizing x and initializing RHS g(x_i) and d
	vec x = linspace<vec>(0, 1, N);
	vec g = pow(h, 2) * 100.0 * exp(-10.0 * x);

	mat L, U, P, A;

	vec row(N);
	row.zeros(); row(0) = 2.0; row(1) = -1.0;
	A = toeplitz(row, row);

	// Starting timer
	clock_t start, finish;
	start = clock();

	lu(L, U, P, A);

	vec y = solve(L, g);
	vec u = solve(U, y);

	//cout << u << endl;

	// Stop timer
	finish = clock();
	cout << "Total time spent on algoritghm with n = " << N << " was " << ((double)(finish - start) / CLOCKS_PER_SEC) * 1000 << "ms" << endl;
	u_data << u;
	u_data.close();

	return u;

}

double errorAnalysis(double* numeric, double* analytic, int exponent) {
	// Takes two pointers numeric and analytic of same size N = 10^exponent as input
	// Returns the maximum relative error between them

	double epsilonMax = -1.0;
	double epsilon;
	for (int i = 1; i < pow(10, exponent)-1; i++) {
		epsilon = fabs((numeric[i] - analytic[i]) / analytic[i]);
		if (epsilon > epsilonMax) {
			epsilonMax = epsilon;
		}	
	}
	return log10(epsilonMax);
}
