// Main class

#pragma once
#include <armadillo>
#include "ExpectationValues.h"
#include "lib.h"

class Ising {
public:
	bool ordered;
	int L;
	int** spinMatrix = (int**) matrix(L,L,sizeof(int));
	double M, E;
	Ising(int L, bool ordered);
	ExpectationValues Metropolis(double T, int mcs);
	double calcE();
	double calcM();
};
