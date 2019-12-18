// Main project 4, made as executable

#include "../include/Ising.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

int main() {
	/*
	Testing implementation and running Ising model for different parameters
	*/

	// Initialization for case 2 (test for different temperatures and random/ordered)
	int nmax;
	ofstream L20T1randomfile;
	ofstream L20T24randomfile;
	ofstream L20T1orderedfile;
	ofstream L20T24orderedfile;
	Ising L20random = Ising(20, false);
	Ising L20ordered = Ising(20, true);
	ofstream ET1file;
	ofstream ET24file;
	ofstream ET1acceptRateFile;
	ofstream ET2acceptRateFile;
	arma::vec ET1;
	arma::vec ET24;

	// Initialization for case 3 (Ising on different temperatures)
	double T1, T2, dT;
	int L;
	string ans = "Y";
	int carlo;
	string filename;
	ofstream ofile;

	int id = -1;
	while (id != 0) {
		cout << "0. Exit" << endl;
		cout << "1. Testing against benchmark values with lattice size L = 2" << endl;
		cout << "2. Lattice size L = 20 for different Monte Carlo cycles" << endl;
		cout << "3. Run for different temperatures, Monte Carlo cycles and lattice sizes" << endl;
		cout << "4. Energy acceptance rate in Metropolis for different Monte Carlos cycles" << endl;
		cout << "Choose id: ";
		cin >> id;
		cout << "\n";
		switch (id) {
		default:
			// Error in input id
			cout << "Bad usage, id = " << id << endl;
			break;
		case 0:
			// End program
			break;
		case 1: {
			// Testing for L = 2 agains benchmark values

      // Starting timer
      auto start = std::chrono::high_resolution_clock::now();

			// Initializing Ising model and running Metropolis algorithm with
			// 100000 Monte Carlo cycles
			Ising test = Ising(2, false);
			ExpectationValues ev = test.Metropolis(1.0, 1000000);

      // Stopping timer
      auto end = std::chrono::high_resolution_clock::now();

			// Printing out result to console
			cout << "<E>=" << ev.expecE << endl;
			cout << "<|M|>=" << ev.expecMabs << endl;
			cout << "<Cv>=" << ev.expecCv << endl;
			cout << "<|Chi|>=" << ev.expecChiAbs << endl;

      // Writing time spent on algorithm to console
      std::chrono::duration<double> duration = end-start;
      cout << "Time: " << duration.count() << "s" << endl;

			break;
		}
		case 2: {
			/*
			Finding expectation values for lattice of size L = 20 and different
			Monte Carlo cylces and writing results to file
			*/

			// Opening resultfiles
			L20T1randomfile.open("L20T1random.txt");
			L20T24randomfile.open("L20T24random.txt");
			L20T1orderedfile.open("L20T1ordered.txt");
			L20T24orderedfile.open("L20T24ordered.txt");
			ET1file.open("ET1data.txt");
			ET24file.open("ET24data.txt");

			// Writing header on resultfiles
			L20T1randomfile << left << setw(15) << "MC-cycles" << left << setw(15)
				<< "<E>"  << left << setw(15)
				<< "<|M|>" << left << setw(15) << "Cv" << left << setw(15) << "ChiAbs" << endl;
			L20T24randomfile << left << setw(15) << "MC-cycles" << left << setw(15)
				<< "<E>"  << left << setw(15)
				<< "<|M|>" << left << setw(15) << "Cv" << left << setw(15) << "ChiAbs" << endl;
			L20T1orderedfile << left << setw(30) << "MC-cycles" << left << setw(30)
				<< "<E>" << left << setw(30) << "<M>" << left << setw(30)
				<< "<|M|>" << left << setw(30) << "Cv" << left
				<< setw(30) << "Chi" << left << setw(30) << "ChiAbs" << endl;
			L20T24orderedfile << left << setw(30) << "MC-cycles" << left << setw(30)
				<< "<E>" << left << setw(12) << "<M>" << left << setw(12)
				<< "<|M|>" << left << setw(12) << "Cv" << left
				<< setw(12) << "Chi" << left << setw(12) << "ChiAbs" << endl;
			ET1file << left << setw(10) << "Energy" << endl;
			ET24file << left << setw(10) << "Energy" << endl;

			// Getting highest Monte Carlo cycles as input
			cout << "Highest power of 10**n, n = ";
			cin >> nmax;

      // Starting timer
      auto start = std::chrono::high_resolution_clock::now();

			// Looping over all exponents of 10 Monte Carlo cycles
			for (int n = 1; n < nmax + 1; n++) {
				int mcs = pow(10, n);
				cout << "Calculating for mcs = " << mcs << " ..." << endl;

				// Metropolis on all Ising models with L = 20
				ExpectationValues out1 = L20random.Metropolis(1.0, mcs);
				ExpectationValues out2 = L20random.Metropolis(2.4, mcs);
				ExpectationValues ordered1 = L20random.Metropolis(1.0, mcs);
				ExpectationValues ordered2 = L20random.Metropolis(2.4, mcs);

				// Writing to resultfiles
				L20T1randomfile << left << setw(15) << mcs << left << setw(15)
					<< out1.expecE  << left << setw(15)
					<< out1.expecMabs << left << setw(15) << out1.expecCv << left
					<< setw(15) << out1.expecChiAbs << endl;
				L20T24randomfile << left << setw(15) << mcs << left << setw(15)
					<< out2.expecE << left << setw(15)
					<< out2.expecMabs << left << setw(15) << out2.expecCv << left
					<< setw(15) << out2.expecChiAbs << endl;
				L20T1orderedfile << left << setw(15) << mcs << left << setw(15)
					<< ordered1.expecE  << left << setw(15)
					<< ordered1.expecMabs << left << setw(15) << ordered1.expecCv << left
					<< setw(15)  << ordered1.expecChiAbs << endl;
				L20T24orderedfile << left << setw(15) << mcs << left << setw(15)
					<< ordered2.expecE  << left << setw(15)
					<< ordered2.expecMabs << left << setw(15) << ordered2.expecCv << left
					<< setw(15)  << ordered2.expecChiAbs << endl;

				// Energy data vectors
				ET1 = out1.Evec;
				ET24 = out2.Evec;
			}

      // Stopping timer and writing duration to console
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = end-start;
      cout << "Time: " << duration.count() << "s" << endl;

			// Writing energy data to file
			for (int i = 0; i < ET1.n_elem; i++) {
				ET1file << left << setw(10) << ET1(i) << endl;
				ET24file << left << setw(10) << ET24(i) << endl;
			}

			// Closing all files
			L20T1randomfile.close();
			L20T24randomfile.close();
			L20T1orderedfile.close();
			L20T24orderedfile.close();
			ET1file.close();
			ET24file.close();
			break;
		}
		case 3: {
			/*
			Running Ising model at different temperatures and finding
			expectation values and writing them to file
			*/

			// Loop as long as long as user want to keep going
			while (ans == "Y" || ans == "y") {

				// Getting parameters to run Metropolis
				cout << "Insert parameters:\nT1 = ";
				cin >> T1;
				cout << "T2 = ";
				cin >> T2;
				cout << "dT = ";
				cin >> dT;
				cout << "Monte Carlo cylces, mcs = ";
				cin >> carlo;
				cout << "L = ";
				cin >> L;
				cout << "Filename: ";
				cin >> filename;

				// Initializing temperature vector and iteration variable
				arma::vec Tvec((T2 - T1) / dT + 3);
				Tvec(0) = T1;
				int i = 0;

				// Opening resultfile and writing header
				ofile.open(filename.c_str());
				ofile << "L = " << L << endl;
				ofile << left << setw(12) << "T" << left << setw(12) << "<E>" << left
					<< setw(12) <<  "<|M|>" << left
					<< setw(12) << "Cv" << left
					<< setw(12) << "ChiAbs" << endl;

        //Starting timer
        auto start = std::chrono::high_resolution_clock::now();

				// Looping over all temperatures. Console output to keep track during simulation
				cout << "Total temperature iterations = " << (int)floor(((T2 - T1) / dT))+1 << endl;
				for (double T = T1; T <= T2; T += dT) {
					cout << "Iteration " << i + 1 << " ..." << endl;

					// Running Metropolis on model
					Ising transition(L, false);
					ExpectationValues phaseTransition = transition.Metropolis(T, carlo);

					// Updating Tvec and iteration variable
					Tvec(i + 1) = Tvec(i) + dT;
					i++;

					// Printing results to file
					ofile << left << setw(15) << Tvec(i - 1) << left << setw(15) << phaseTransition.expecE << left << setw(15)
						<< phaseTransition.expecMabs
						<< left << setw(15) << phaseTransition.expecCv << left << setw(15)
						<< phaseTransition.expecChiAbs << endl;
				}

        // Stopping timer and writing duration to console
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end-start;
        cout << "Time: " << duration.count() << "s" << endl;

				// Close result-file
				ofile.close();

				cout << "\nDo you want to run again for different parameters? [Y/N]: ";
				cin >> ans;
			}
			break;
		}
		case 4: {
			/*
			Counting number of accepted energies in Metropolis for L = 20
			and different Monte Carlo cylces. T = 1, 2.4. Writing results to file.
			*/

			// Opening resultfiles
			ET1acceptRateFile.open("ET1accepted.txt");
			ET2acceptRateFile.open("ET24accepted.txt");

			// Writing header on resultfiles
			ET1acceptRateFile << left << setw(12) << "MC-cycles" << left << setw(12)
				<< "Accepted" << endl;
			ET2acceptRateFile << left << setw(12) << "MC-cycles" << left << setw(12)
				<< "Accepted" << endl;

			// Getting highest Monte Carlo cycles as input
			cout << "Highest power of 10**n, n = ";
			cin >> nmax;

      // Starting timer
      auto start = std::chrono::high_resolution_clock::now();

			// Looping over all exponents of 10 Monte Carlo cycles
			for (int n = 1; n < nmax + 1; n++) {
				int mcs = pow(10, n);
				cout << "Calculating for mcs = " << mcs << " ..." << endl;

				// Metropolis on all Ising models with L = 20
				ExpectationValues out1 = L20random.Metropolis(1.0, mcs);
				ExpectationValues out2 = L20random.Metropolis(2.4, mcs);

				// Writing to resultfiles
				ET1acceptRateFile << left << setw(12) << mcs << left << setw(12)
					<< out1.accepted << endl;
				ET2acceptRateFile << left << setw(12) << mcs << left << setw(12)
					<< out2.accepted << endl;
			}

      // Stopping timer and writing duration to console
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = end - start;
      cout << "Time: " << duration.count() << "s" << endl;

      // Closing files
			ET1acceptRateFile.close();
			ET1acceptRateFile.close();
			break;
		}
		}
	}
  return 0;
}
