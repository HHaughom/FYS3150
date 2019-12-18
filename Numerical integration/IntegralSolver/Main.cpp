/*
   Program to solve the integral arising from the expectation value of
   the correlation energy of the two atoms in a Helium atom.
   Method choices 
       1. Gauss-Legendre Quadrature
	   2. Gauss-Laguerre Quadrature
	          * Added option of parallelisation
	   3. Monte-Carlo, uniform sampling
	   4. Monte-Carlo, importance sampling 

   Takes number of steps N as input.
   Can choose to run for one N, or write a sequence of N to file.
*/

#include "Integrate.h"

int main() {

	int N;
	double integral, error, variance;
	double* result;

	double analytic = 5 * pow(PI, 2) / pow(16, 2);
	int choice = -1;

	char write = 'N';
	char parallel = 'N';

	while (choice != 0) {
		cout << "0. Exit" << endl;
		cout << "1. Gauss-Legendre" << endl;
		cout << "2. Gauss-Laguerre" << endl;
		cout << "3. Monte Carlo, uniform sampling" << endl;
		cout << "4. Monte Carlo, importance sampling" << endl;
		cout << "Choose integration method: ";
		cin >> choice;
		cout << endl;

		switch (choice) {
		case 0:
			break;

		case 1: {
			cout << "Integrate with Gauss-Legrende" << endl;
			cout << "Choose number of integration points, N = ";
			cin >> N;
			cout << "Write to file [Y/N]? ";
			cin >> write;
			cout << endl;
			
			if (write == 'Y') {
				ofstream file;
				file.open("gaussLeg.dat");
				file << "N  " << "Integral  " << "Error  " << endl;

				for (int i = 1; i <= N / 5; i++) {
					integral = gaussLegendre(5 * i);
					error = abs(integral - analytic) / analytic;
					file << 5 * i << "  " << integral << "  " << error << endl;
				}

				file.close();
			}

			else {
				integral = gaussLegendre(N);
				error = abs(integral - analytic) / analytic;
				cout << "Analytic: " << analytic << endl;
				cout << "Numeric: " << integral << endl;
				cout << "Error: " << error << endl;
			}
			break;
		}	
		case 2: {
			cout << "Integrate with Gauss-Laguerre" << endl;
			cout << "Choose number of integration points, N = ";
			cin >> N;
			cout << "Parallelize [Y/N]?";
			cin >> parallel;
			cout << "Write to file [Y/N]?";
			cin >> write;
			cout << endl;

			if (write == 'Y') {
				ofstream file;
				file.open("gaussLag.dat");
				file << "N  " << "Integral  " << "Error  " << endl;

				for (int i = 1; i <= N / 5; i++) {
					if (parallel == 'Y') {
						integral = gaussLaguerrePara(5 * i);
					}
					else {
						integral = gaussLaguerre(5 * i);
					}
					error = abs(integral - analytic) / analytic;
					file << 5 * i << "  " << integral << "  " << error << endl;
				}
				file.close();
			}

			else {
				if (parallel == 'Y') {
					integral = gaussLaguerrePara(N);
				}
				else {
					integral = gaussLaguerre(N);
				}
				error = abs(integral - analytic) / analytic;
				cout << "Analytic: " << analytic << endl;
				cout << "Numeric: " << integral << endl;
				cout << "Error: " << error << endl;
			}
			break;
		}
		case 3: {
			cout << "Integrate with brute force Monte Carlo" << endl;
			cout << "Choose number of integration points, N = ";
			cin >> N;
			cout << "Write to file [Y/N]? ";
			cin >> write;
			cout << endl;

			if (write == 'Y') {
				ofstream file;
				file.open("bruteMC.dat");
				file << "N  " << "Integral  " <<"Error  "<< "Variance" << endl;

				for (int i = 1; i < log10(N)+1; i++) {
					result = bruteMonteCarlo(pow(10, i)); 
					integral = result[0];
					variance = result[1];
					error = abs(integral - analytic) / analytic;
					file << pow(10, i) << "  " << integral << "  " << error << "  " << variance << endl;
				}

				file.close();
			}

			else {
				result = bruteMonteCarlo(N);
				integral = result[0];
				variance = result[1];
				error = abs(integral - analytic) / analytic;

				cout << "Analytic: " << analytic << endl;
				cout << "Numeric: " << integral << endl;
				cout << "Error: " << error << endl;
				cout << "Variance: " << variance << endl;
			}

			break;
		}
		case 4: {
			cout << "Integrate with importance sampling Monte Carlo" << endl;
			cout << "Choose number of integration points, N = ";
			cin >> N;
			cout << "Write to file [Y/N]? ";
			cin >> write;
			cout << endl;

			if (write == 'Y') {
				ofstream file;
				file.open("impsampMC.dat");
				file << "N  " << "Integral  " << "Error " << "Variance" << endl;

				for (int i = 1; i < log10(N)+1; i++) {
					result = impsampMonteCarlo(pow(10, i));
					integral = result[0];
					variance = result[1];
					error = abs(integral - analytic) / analytic;
					file << pow(10, i) << "  " << integral << "  " << error << "  " << variance << endl;
				}

				file.close();
			}

			else {
				result = impsampMonteCarlo(N);
				result = impsampMonteCarlo(N);
				integral = result[0];
				variance = result[1];
				error = abs(integral - analytic) / analytic;

				cout << "Analytic: " << analytic << endl;
				cout << "Numeric: " << integral << endl;
				cout << "Error: " << error << endl;
				cout << "Variance: " << variance << endl;
			}
			break;
		}
		}
	}

	return 0;
}
