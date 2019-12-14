// Class cpp

#include "../include/Ising.h"
#include <iostream>


using std::cout;
using std::endl;

Ising::Ising(int L, bool ordered) : L(L), ordered(ordered) {
  /*
  Function (constructor) to initialize a lattice of size (LxL) with
  ordered/random spin at every position
  */
	if (ordered) {
    // Setting up ordered spin (every spin up)
		for(int i = 0; i < L; i++) {
			for(int j = 0; j < L; j++) {
				spinMatrix[i][j] = 1;
			}
		}
	}
	else {
    // Setting up random spins
		long idum = -1;
		for (int i = 0; i < L; i++) {
			for (int j = 0; j < L; j++) {
				spinMatrix[i][j] = ran1(&idum);
				if (spinMatrix[i][j] > 0.5) {
					spinMatrix[i][j] = 1;
				}
				else {
					spinMatrix[i][j] = -1;
				}
			}
		}
	}
  // Calculating initial magnetization and energy of configuration
	M = calcM();
	E = calcE();

}


double Ising::calcE() {
	/*
	Function to calculate energy of configuration
	*/
	double energy = 0;
  // Looping over all spins
	for (int y = 0; y < L; y++) {
		for (int x = 0; x < L; x++) {
			energy -= spinMatrix[y][x] * (spinMatrix[(y + L - 1) % L][x] + spinMatrix[y][(x + L - 1) % L]);
		}
	}
	return energy;
}

double Ising::calcM() {
  /*
  Function to calculate magnetization of configuration
  */
	double magnet = 0;
  // Looping over all spins
	for (int y = 0; y < L; y++) {
		for (int x = 0; x < L; x++) {
			magnet += spinMatrix[y][x];
		}
	}
	return magnet;
}

ExpectationValues Ising::Metropolis(double T, int mcs) {
  /*
  Function to perform the Metropolis algorithm on the lattice at temperature
  T with 'mcs' Monte Carlo cycles
  */
  // Initialize variables
	double sumE = 0.0;
	double sumE2 = 0.0;
	double sumM2 = 0.0;
	double sumMabs = 0.0;
	double dE, dM;
	int acceptedEnergies = 0;
	long idum = -1;

	// Energi-vector for finding distribution
	arma::vec Etemp(mcs+1);
	Etemp(0) = calcE();

  // Initialize index and probabilities
	int x, y;
	double w, r;

  // Total amount of spins in lattice
	int spins = L * L;

	for (int cycle = 0; cycle < mcs; cycle++) {
		// Loop over all monte-carlo cycles

		// Generate random index in lattice
		x = (int) (ran1(&idum)*(double)L);
		y = (int) (ran1(&idum)*(double)L);

		// Flip spin at index (y,x)
		spinMatrix[y][x] *= -1.0;

		// Energy-change from flip
		dE = -(double)2 * spinMatrix[y][x] *
			(spinMatrix[y][(x + L - 1) % L] + spinMatrix[y][(x + L + 1) % L]
				+ spinMatrix[(y + L - 1) % L][x] + spinMatrix[(y + L + 1) % L][x]);

		dM = (double)2 * spinMatrix[y][x];

		// Accept new configuration if energy is lower than initial state
		if (dE <= 0) {
			E += dE;
			M += dM;
			acceptedEnergies++;
		}
		else {
			// Accept with probability given by Boltzmann distribution
			w = exp(-dE / T);
			r = ran1(&idum);
			if (r <= w) {
				// Accept new configuration
				E += dE;
				M += dM;
				acceptedEnergies++;
			}
			else {
				// Flip back, unchanged energy (not accepted flip)
				spinMatrix[y][x] *= -1;
			}
		}
		Etemp(cycle+1) = E;
		sumE += E;
		sumE2 += E * E;
		sumM2 += M * M;
		sumMabs += fabs(M);
	}

	// Calculating expectation values
	sumE /= mcs;
	sumMabs /= mcs;
	sumE2 /= mcs;
	sumM2 /= mcs;

	double varE = (sumE2 - (sumE * sumE));
	double varMabs = (sumM2 - (sumMabs* sumMabs ));

	ExpectationValues out;

	out.expecE = sumE / spins;
	out.expecMabs = sumMabs / spins;
	out.expecCv = varE / T / T / spins;
	out.expecChiAbs = varMabs / T / spins;

	out.accepted = acceptedEnergies;
  out.Evec = Etemp;

	return out;
}
