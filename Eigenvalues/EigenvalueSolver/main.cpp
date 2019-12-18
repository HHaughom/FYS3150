#include "Functions.h"
#include "Timer.h"
#include <iostream>
#include <armadillo>

using namespace std;

int main() {

	// Initialization
	double dia, offDia, rhoMax;
	arma::vec addDia, rho;
	arma::vec omega(4);
	int n;
	arma::vec analyticSolve, armadilloSolve, jacobiSolve;
	arma::mat A;

	int ans = -1;
	while (ans != 0) {
		cout << "0. Exit" << endl;
		cout << "1. Analytic solver" << endl;
		cout << "2. Armadillo solver" << endl;
		cout << "3. Jacobi's method" << endl;
		cout << "4. Quantum dots in three dimensions, one electron" << endl;
		cout << "5. Quantum dots in three dimensions, two electrons" << endl;
		cout << "Choose method: ";
		cin >> ans;
		cout << endl;

		switch (ans) {
		case 0:
			break;

		case 1:
			cout << "This tests the analytic solver for finding the eigenvalues of the Toeplitz" << endl;
			cout << "matrix with diagonal elements = 2.0/h^2 and off-diagonal elements = -1.0/h^2" << endl;
			cout << "Choose dimension of Toeplitz matrix (n x n) n = ";
			cin >> n;

			dia = 2.0 * pow(n, 2);
			offDia = -1.0 * pow(n, 2);
			analyticSolve = analytic_beam_eigenvalues(dia, offDia, n);

			// Only printing out the first 10 eigenvalues
			cout << "The eigenvalues of the matrix are:" << endl;
			for (int i = 0; i < 10; i++) {
				cout << analyticSolve[i] << endl;
			}
			cout << endl;
			break;

		case 2:
			cout << "This tests the armadillo solver for finding the eigenvalues of the Toeplitz" << endl;
			cout << "matrix with diagonal elements = 2.0/h^2 and off-diagonal elements = -1.0/h^2" << endl;
			cout << "Choose dimension of Toeplitz matrix (n x n) n = ";
			cin >> n;

			addDia = arma::zeros<arma::vec>(n);
			dia = 2.0 * pow(n, 2);
			offDia = -1.0 * pow(n, 2);
			armadilloSolve = armadillo_eigenvalues(n, dia, offDia, addDia);

			// Only printing out the first 10 eigenvalues
			cout << "The eigenvalues of the matrix are:" << endl;
			for (int i = 0; i < 10; i++) {
				cout << armadilloSolve[i] << endl;
			}
			cout << endl;
			break;

		case 3:
			cout << "This tests the implementation of the Jacobi's method and finds the eigenvalues" << endl;
			cout << "of the Toeplitz matrix with diagonal elements = 2.0/h^2 and off-diagonal elements = -1.0/h^2" << endl;
			cout << "Choose dimension of the Toeplitz matrix (n x n) n = ";
			cin >> n;

			addDia = arma::zeros<arma::vec>(n);
			dia = 2.0 * pow(n, 2);
			offDia = -1.0 * pow(n, 2);
			jacobiSolve = sort(JacobiSolver(1.0E-10, n, 100000, dia, offDia, addDia));

			// Only printing out the first 10 eigenvalues
			cout << "The eigenvalues of the matrix are:" << endl;
			for (int i = 0; i < 10; i++) {
				cout << jacobiSolve[i] << endl;
			}
			cout << endl;
			break;

		case 4:
			cout << "This tests the implementation of the Jacobi's method and finds the eigenvalues" << endl;
			cout << "of the Tridiagonal matrix with found in quantum dots in 3 dimensions with 1 electron" << endl;
			cout << "Choose dimension of the tridiagonal matrix (n x n) n = ";
			cin >> n;
			cout << "Choose maximum rho, rhoMax = ";
			cin >> rhoMax;

			rho = arma::regspace(1, n) * rhoMax / n;
			addDia = pow(rho, 2);
			dia = 2.0 * pow(n / rhoMax, 2);
			offDia = -1.0 * pow(n / rhoMax, 2);
			jacobiSolve = sort(JacobiSolver(1.0E-10, n, 100000, dia, offDia, addDia));

			// Only printing out first 10 eigenvalues
			cout << "The eigenvalues of the matrix are:" << endl;
			for (int i = 0; i < 10; i++) {
				cout << jacobiSolve[i] << endl;
			}
			cout << endl;
			break;
			
		case 5:

			cout << "This tests the implementation of the Jacobi's method and finds the eigenvalues" << endl;
			cout << "of the Tridiagonal matrix with found in quantum dots in 3 dimensions with 2 electrons" << endl;
			cout << "Choose dimension of the tridiagonal matrix (n x n) n = ";
			cin >> n;
			cout << "Choose maximum rho, rhoMax = ";
			cin >> rhoMax;

			omega[0] = 0.01;
			omega[1] = 0.5;
			omega[2] = 1;
			omega[3] = 5;
			rho = arma::regspace(1, n) * rhoMax / n;

			// Solving for all values of omega
			for (int i = 0; i < 4; i++) {
				cout << "Omega = " << omega[i] << ":" << endl;
				addDia = pow(omega[i] * rho, 2) + 1 / rho;
				dia = 2.0 * pow(n / rhoMax, 2);
				offDia = -1.0 * pow(n / rhoMax, 2);
				jacobiSolve = sort(JacobiSolver(1.0E-10, n, 100000, dia, offDia, addDia));

				// Only printing out first 3 eigenvalues
				for (int j = 0; j < 3; j++) {
					cout << jacobiSolve[j] << endl;
				}
			}
			break;

		default:
			cout << endl;
			cout << "Bad usage, id = " << ans << endl;
			break;
		}
	}
	return 0;
}
