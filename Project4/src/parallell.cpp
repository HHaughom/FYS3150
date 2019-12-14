// Parallellization with MPI

#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <chrono>
#include "../include/Ising.h"
#include "../include/lib.h"
#include "Ising.cpp"
#include "lib.cpp"


int main(int argc, char* argv[]) {

  int numProcs, procId;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Initialize different Ising-models with different lattice sizes
  // We have commented out the L100 case since this is testet on a
  // computer with 4 cores. This is to see maximum speedup
  Ising L20 = Ising(20, false);
  Ising L40 = Ising(40, false);
  Ising L60 = Ising(60, false);
  Ising L80 = Ising(80, false);
  //Ising L100 = Ising(100, false);

  // Setting temperature interval and Monte Carlo cylces
  double T1 = 2.0;
  double T2 = 3.0;
  double dT = 0.005;
  int mcs = 1000000;

  // Starting timer
  auto start = std::chrono::high_resolution_clock::now();

  // Looping over all T
  for(double T = T1; T <= T2; T+=dT) {

    std::cout << T << " " << procId << std::endl;

    L20.Metropolis(T, mcs);
    L40.Metropolis(T, mcs);
    L60.Metropolis(T, mcs);
    L80.Metropolis(T, mcs);
    //L100.Metropolis(T, mcs);
  }

  // Stopping timer and writing duration to console
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  std::cout << "Time: " << duration.count() << "s" << std::endl;

  MPI_Finalize();
  return 0;
}
