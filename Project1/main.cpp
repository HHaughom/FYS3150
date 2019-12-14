#include<iostream>
#include<cstdlib>
#include"Algorithms.h"
#include "armadillo"

using namespace arma;
using namespace std;

int main(int argc, char* argv[]) {
	int n = atof(argv[1]);
	int id = -1;
	while (id != 0) {
		cout << "0. Exit" << endl;
		cout << "1. 9n algorithm" << endl;
		cout << "2. 4n algorithm" << endl;
		cout << "3. LU armadillo" << endl;
		cout << "4. Error analysis" << endl;
		cout << "Choose what you want to test: ";
		cin >> id;
		cout << endl;

		switch (id) {
		case 0:
			break;
		case 1: {
			for (int i = 1; i <= n; i++) {
				int N = pow(10, i);
				double* solutionSlow = new double[N];
				solutionSlow = algo_9n(i);
				delete[] solutionSlow;
			}
			break;
		}
		case 2: {
			for (int i = 1; i <= n; i++) {
				int N = pow(10, i);
				double* solutionFast = new double[N];
				solutionFast = algo_4n(i);
				delete[] solutionFast;
			}
			break;
		}
		case 3: {
			for (int i = 1; i <= n; i++) {
				vec solutionArma;
				solutionArma = arma_lu(i);
			}
			break;
		}
		case 4: {

			ofstream errorFile1, errorFile2, errorFile3;
			errorFile1.open("relativeError_slow.dat");
			errorFile2.open("relativeError_fast.dat");
			errorFile3.open("relativeError_arma.dat");

			for (int j = 1; j < 8; j++) {

				double* analytic = new double[pow(10, j)];
				double* numeric_slow = new double[pow(10, j)];
				double* numeric_fast = new double[pow(10, j)];

				numeric_slow = algo_9n(j);
				numeric_fast = algo_4n(j);

				double x;

				for (int i = 0; i < pow(10, j); i++) {
					x = i * (1.0 / pow(10, j));
					//int n = pow(10, j) - i - 1;
					analytic[i] = analytical(x);
				}

				errorFile1 << errorAnalysis(numeric_slow, analytic, j) << endl;
				errorFile2 << errorAnalysis(numeric_fast, analytic, j) << endl;


				if (j < 5) {
					vec numeric_arma = arma_lu(j);
					double* numeric_arma_arr = new double[pow(10, j)];

					for (int i = 0; i < pow(10, j); i++) {
						//int n = pow(10, j) - i - 1;
						numeric_arma_arr[i] = numeric_arma(i);
					}

					errorFile3 << errorAnalysis(numeric_arma_arr, analytic, j) << endl;
					delete[] numeric_arma_arr;
				}	

				delete[] analytic; delete[] numeric_slow; delete[] numeric_fast;
			}

			errorFile1.close(); errorFile1.close(); errorFile3.close();
			break;
		}
		default:
			cout << "Bad usage: id = " << id << endl;
			break;
		}
	}

	return 0;
}
