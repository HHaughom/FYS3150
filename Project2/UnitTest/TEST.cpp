#include "pch.h"
#include "CppUnitTest.h"
#include "../testTest/Functions.h"
#include "../testTest/Functions.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest {

	TEST_CLASS(TriDiagonal) {
		/*
		Test-class for setting up tridiagonal matrices and finding
		max values of off-diagonal elements.
		*/
public:
	TEST_METHOD(initTriDiagonal) {
		// Initialization
		int n = 5;
		double dia = 2.0;
		double offDia = -1.0;
		arma::vec addDia = arma::zeros<arma::vec>(n);

		// Generating matrix
		arma::mat A = setupTriDiag(n, dia, offDia, addDia);
		arma::vec row(n);
		row.zeros();
		row(0) = dia;
		row(1) = offDia;

		// Expected matrix
		arma::mat expected = arma::toeplitz(row, row);

		// Checking for equality
		bool success = true;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (A(i, j) != expected(i, j)) {
					success = false;
				}
			}
		}
		// Assert result
		Assert::IsTrue(success);
	}

	TEST_METHOD(maxOffNegativeDia) {
		// Initialization
		int n = 5;
		double dia = 2.0;
		double offDia = -1.0;
		arma::vec addDia = arma::zeros<arma::vec>(n);
		arma::mat A = setupTriDiag(n, dia, offDia, addDia);

		// Inserting a max value at index (1,2)
		A(1, 2) = -3.0;

		// Generating index-result
		int p, q;
		int* generated = offDiag(A, &p, &q, n);

		// Expected result
		int expected[2];
		expected[0] = 1;
		expected[1] = 2;

		// Checking for equality
		bool success = true;
		for (int i = 0; i < 2; i++) {
			if (generated[i] != expected[i]) {
				success = false;
			}
		}
		// Assert result
		Assert::IsTrue(success);
	}

	TEST_METHOD(maxOffPositiveDia) {
		// Initialization
		int n = 5;
		double dia = 2.0;
		double offDia = 1.0;
		arma::vec addDia = arma::zeros<arma::vec>(n);
		arma::mat A = setupTriDiag(n, dia, offDia, addDia);

		// Inserting max value at index (1,2)
		A(1, 2) = 3.0;

		// Generating index-result
		int p, q;
		int* generated = offDiag(A, &p, &q, n);

		// Expected result
		int expected[2];
		expected[0] = 1;
		expected[1] = 2;

		// Checking for equality
		bool success = true;
		for (int i = 0; i < 2; i++) {
			if (generated[i] != expected[i]) {
				success = false;
			}
		}
		// Assert result
		Assert::IsTrue(success);
	}
	};

	TEST_CLASS(Eigenvalues) {
		/*
		Test-class for finding eigenvalues of tridiagonal matrices
		using different algorithms.
		*/
public:
	TEST_METHOD(Analytic) {
		// Initialization
		int n = 5;
		double dia = 2.0;
		double offDia = -1.0;

		// Expected result
		arma::vec expected(n);
		expected[0] = 0.2679;
		expected[1] = 1.0000;
		expected[2] = 2.0000;
		expected[3] = 3.0000;
		expected[4] = 3.7321;

		// Generating result using analytic-solver
		arma::vec generated = sort(analytic_beam_eigenvalues(dia, offDia, n));

		// Checking for equality with a given error margin (tol)
		bool success = true;
		double tol = 1E-4;
		for (int i = 0; i < n; i++) {
			if (fabs(expected[i] - generated[i]) > tol) {
				success = false;
			}
		}
		// Assert result
		Assert::IsTrue(success);
	}

	TEST_METHOD(Armadillo) {
		// Initialization
		int n = 5;
		double dia = 2.0;
		double offDia = -1.0;
		arma::vec addDia = arma::zeros<arma::vec>(n);
		arma::mat A = setupTriDiag(n, dia, offDia, addDia);

		// Expected result
		arma::vec expected(n);
		expected[0] = 0.2679;
		expected[1] = 1.0000;
		expected[2] = 2.0000;
		expected[3] = 3.0000;
		expected[4] = 3.7321;

		// Generating result using armadillo-solver
		arma::vec generated = sort(armadillo_eigenvalues(n, dia, offDia, addDia));

		// Checking for equality with given error margin (tol)
		bool success = true;
		double tol = 1E-4;
		for (int i = 0; i < n; i++) {
			if (fabs(expected[i] - generated[i]) > tol) {
				success = false;
			}
		}
		// Assert result
		Assert::IsTrue(success);
	}

	TEST_METHOD(Jacobi) {
		// Initialization
		int n = 5;
		double dia = 2.0;
		double offDia = -1.0;
		arma::vec addDia = arma::zeros<arma::vec>(n);

		// Expected result
		arma::vec expected(n);
		expected[0] = 0.2679;
		expected[1] = 1.0000;
		expected[2] = 2.0000;
		expected[3] = 3.0000;
		expected[4] = 3.7321;

		// Generating result using jacobi-solver
		arma::vec generated = sort(JacobiSolver(1E-10, n, 10000, dia, offDia, addDia));

		// Checking for equality with given error margin (tol)
		bool success = true;
		double tol = 1E-4;
		for (int i = 0; i < n; i++) {
			if (fabs(expected[i] - generated[i]) > tol) {
				success = false;
			}
		}
		// Assert result
		Assert::IsTrue(success);
	}
	};
}
