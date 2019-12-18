#include "Functions.h"
#include "Timer.h"
#define _USE_MATH_DEFINES
#include <math.h>

int* offDiag(arma::mat A, int* p, int* q, int n) {
	/*
	Function for finding index off max-element of a tridiagonal matrix
	*/

	// Initialization
	int index[2];
	double max = 0.0;

	// Searching for max element
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			double aij = fabs(A(i, j));
			if (aij > max)
			{
				max = aij;  *p = i; *q = j;
			}
		}
	}

	// Setting and returning max-index
	index[0] = *p;
	index[1] = *q;
	return index;
}

arma::mat Rotate(arma::mat A, arma::mat R, int k, int l, int n) {
	/*
	Function for rotating a matrix A with dimensions (n x n) so that element (k,l) becomes 0.
	Returns the matrix A after rotation.
	*/

	// Initialization
	double s, c;

	// Initalization of parameters when element (k,l) is not 0
	if (A(k, l) != 0.0) {
		double t, tau;
		tau = (A(l, l) - A(k, k)) / (2 * A(k, l));

		if (tau >= 0) {
			t = 1.0 / (tau + sqrt(1.0 + tau * tau));
		}
		else {
			t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
		}

		c = 1 / sqrt(1 + t * t);
		s = c * t;
	}
	// Initialization of parameters when element (k,l) is 0
	else {
		c = 1.0;
		s = 0.0;
	}
	// Initialization for rotation
	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = A(k, k);
	a_ll = A(l, l);
	A(k, k) = c * c * a_kk - 2.0 * c * s * A(k, l) + s * s * a_ll;
	A(l, l) = s * s * a_kk + 2.0 * c * s * A(k, l) + c * c * a_ll;
	A(k, l) = 0.0;
	A(l, k) = 0.0;

	// Executing rotation
	for (int i = 0; i < n; i++) {
		if (i != k && i != l) {
			a_ik = A(i, k);
			a_il = A(i, l);
			A(i, k) = c * a_ik - s * a_il;
			A(k, i) = A(i, k);
			A(i, l) = c * a_il + s * a_ik;
			A(l, i) = A(i, l);
		}
		r_ik = R(i, k);
		r_il = R(i, l);

		// Setting eigenvectors
		R(i, k) = c * r_ik - s * r_il;
		R(i, l) = c * r_il + s * r_ik;
	}
	// Returning rotated matrix A
	return A;
}

arma::vec analytic_beam_eigenvalues(double diagonal, double offDiagonal, int n) {
	/*
	Function for calculating analytic eigenvalues of a tridiagonal matrix with
	dimensions (n x n). Returns vector of eigenvalues
	*/

	// Initialization
	arma::vec eigenvals;
	arma::vec j = arma::regspace(1, n);

	// Starting timer (see Timer-class)
	Timer timer;

	// Calculating eigenvalues
	eigenvals = diagonal + 2 * offDiagonal * cos(j * M_PI / (n + 1.0));

	// Returning eigenvalues
	return eigenvals;
}

arma::vec armadillo_eigenvalues(int n, double dia, double offDia, arma::vec addDia) {
	/*
	Function for finding eigenvalues of tridiagonal matrix using armadillo functions.
	Matrix has dimensions (n x n) and diagonal elements = dia + addDia. Returns vector of eigenvalues
	*/

	// Initialization
	arma::mat A = setupTriDiag(n, dia, offDia, addDia);
	arma::vec eigenvalues;
	arma::mat eigenvectors;

	// Starting timer (see Timer-class)
	Timer timer;

	 // Calculating eigenvalues/eigenvectors
	arma::eig_sym(eigenvalues, eigenvectors, A);

	// Returning eigenvalues
	return eigenvalues;
}

arma::vec JacobiSolver(double tol, int n, int maxiter, double dia, double offDia, arma::vec addDia) {
	/*
	Function for finding eigenvalues of a tridiagonal matrix with dimensions (n x n).
	Matrix has diagonal elements = dia + addDia. The function returns vector of eigenvalues
	*/

	// Initialization
	arma::mat B = setupTriDiag(n, dia, offDia, addDia);
	arma::mat R(n, n);
	arma::vec eigenvalues(n);
	R.eye();
	int maxNonDiag = 1.0;
	int iterations = 0;

	// Starting timer (see Timer-class)
	Timer timer;

	// Rotating matrix while maxNonDiag is greater than a given error margin (tol)
	// Also capped loop for a max amount of iterations
	while (maxNonDiag > tol && iterations <= maxiter) {
		int p, q;
		offDiag(B, &p, &q, n);
		B = Rotate(B, R, p, q, n);
		iterations++;
	}

	// Setting and returning eigenvalues
	for (int i = 0; i < n; i++) {
		eigenvalues[i] = B(i, i);
	}
	return eigenvalues;
}

arma::mat setupTriDiag(int n, double diagonal, double offDiagonal, arma::vec addDiagonal) {
	/*
	Function for setting up a tridiagonal matrix with dimensions (n x n), diagonal elements = dia + addDiagonal
	and off-diagonal elements offDiagonal. Returns tridiagonal matrix
	*/

	// Initialization
	arma::mat A;
	arma::mat I;
	arma::vec row(n);
	row.zeros();
	row(0) = diagonal;
	row(1) = offDiagonal;

	// Setting up matrix
	A = arma::toeplitz(row, row);
	I = arma::diagmat(addDiagonal);

	// Returning tridiagonal matrix
	return A + I;
}
