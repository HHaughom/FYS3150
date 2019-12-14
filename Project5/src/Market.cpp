// File for implementation of market functions

#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../include/Market.h"
#include "../include/tools.h"

// For debugging
#include <iostream>

// Base class (Market)
Market::Market(int N, double m0, double dm, double maxMoney): N(N), m0(m0), dm(dm), maxMoney(maxMoney), bins((int)maxMoney/dm) {
  /*
  Constructor initializing the Market
  */
  // Allocating memory for data in Market object
  this->agents = new double[this->N];
  this->dist = new double[this->bins];
  this->meanDist = new double[this->bins];

  // Filling agents with initial income m0
  for(int i = 0; i < N; i++) {
    this->agents[i] = m0;
  }

  // Initializing distribution
  int* temp = distribution(this->agents, this->N, this->maxMoney, this->dm);
  this->dist = normalize(temp, this->bins, this->N);

  // Initializing meanDistribution
  for(int i = 0; i < this->bins; i++) {
    this->meanDist[i] = 0.0;
  }

  // Freeing up memory
  delete [] temp;
}

int Market::transactions(double diffLimit, double saving) {
  /*
  Function to perform transactions on a market until equilibrium
  has been achieved. Returns transactions needed.
  */
  // Initializing agents and transaction variables
  int agent1, agent2;
  double epsilon, mi, mj, delta;

  // Initializing amount of transaction needed to reach equilibrium
  int eq = 0;
  int maxTransactions = 10000000;

  // Initializing difference
  double diff = diffLimit + 1;

  // Seed for pseudo-random number generator
  std::srand(std::time(0));

  // Perform transactions as long as equilibrium is not achieved
  while(diff >= diffLimit && eq < maxTransactions) {
    // Draw 2 random agents to make a transaction
    agent1 = std::rand() % this->N;
    agent2 = std::rand() % this->N;

    // If both agents are equal, draw a new agent until they are different
    while(agent1 == agent2) {
      agent2 = std::rand() % this->N;
    }

    // Draw epsilon from uniform distribution [0,1]
    epsilon = (double) std::rand() / RAND_MAX;

    // Calculating variables for transaction
    mi = this->agents[agent1];
    mj = this->agents[agent2];
    delta = (1 - saving) * (epsilon*mj - (1 - epsilon)*mi);

    // Update income of agents and incrementing k
    this->agents[agent1] += delta;
    this->agents[agent2] -= delta;
    eq++;

    // Finding new distribution and calculating variance between old and new
    // Only do this with a frequency of 100000 transactions
    if(eq % 100000 == 0) {
      // Initializing new distributions
      int* newDist = new int[this->bins];
      double* newNormalized = new double[this->bins];

      newDist = distribution(this->agents, this->N, this->maxMoney, this->dm);
      newNormalized = normalize(newDist, this->bins, this->N);
      diff = difference(this->dist, newNormalized, this->bins);

      // Updating distribution
      this->dist = newNormalized;

      // Freeing up memory
      delete [] newDist;
      delete [] newNormalized;
    }
  }

  // Returning transactions needed to reach equilibrium
  return eq;
}

void Market::simulateMarket(int runs, double saving, double diffLimit) {
    /*
    Function for simulating market with 'runs' runs and saving
    quotient 'saving'
    */
    // Initializing transactions variable
    int n;

    // Looping for all runs
    for(int i = 0; i < runs; i++) {
      // Print to console to keep track during execution
      std::cout << "Run nr. " << i+1 << std::endl;

      // Performing transactions
      n = this->transactions(diffLimit, saving);

      // Updating meanDist
      for(int j = 0; j < this->bins; j++) {
        this->meanDist[j] += this->dist[j];
      }
    }

    // Taking mean
    for(int i = 0; i < this->bins; i++) {
      this->meanDist[i] /= runs;
    }
}

void Market::writeToFile(std::string filename) {
  /*
  Function to write data in distribution to file 'filename'
  */
  // Initializing and opening file
  std::ofstream file;
  file.open(filename);

  // Writing data in file
  for(int i = 0; i < this->bins; i++) {
    file << std::left << std::setw(10) << this->meanDist[i] << std::endl;
  }

  // Closing file
  file.close();
}

// Neighbour class
int NeighboursMarket::transactions(double diffLimit, double saving, double alfa) {
  /*
  Function to perform transactions until equilibrium is achieved
  Here the probability of two agents to perform an exchange is dependent
  on the alfa parameter
  */
  int agent1, agent2;
  double mi, mj, epsilon, dm, a, prob, delta;

  int eq = 0;
  int maxTransactions = 10000000;
  double diff = diffLimit + 1;

  // Seed for pseudo-random number generator
  std::srand(std::time(0));

  // Perform transactions as long as equilibrium is not achieved
  while(diff >= diffLimit && eq < maxTransactions) {
    // Draw two random agents
    agent1 = std::rand() % this->N;
    agent2 = std::rand() % this->N;

    // Draw again if agent1 == agent2 until they are different
    while(agent1 == agent2) {
      agent2 = std::rand() % this->N;
    }
    mi = this->agents[agent1];
    mj = this->agents[agent2];

    // Draw random nubmer in unifrom distribution to decide whether
    // transaction happens
    a = (double) std::rand() / RAND_MAX;
    // Calculate probability of transaction to happen
    prob = pow(fabs(mi - mj),-alfa);

    // Perform transaction if a <= prob. If not, transaction doesn't happen
    if(a <= prob) {
      // Draw random number epsilon from uniform distribution
      epsilon = (double) std::rand() / RAND_MAX;

      // Calculate transaction change
      delta = (1 - saving) * (epsilon*mj - (1 - epsilon)*mi);

      // Update income of agents and incrementing eq
      this->agents[agent1] += delta;
      this->agents[agent2] -= delta;
      eq++;

      // Finding new distribution and calculating variance between old and new
      // Only do this with a frequency of 100000 transactions
      if(eq % 100000 == 0) {
        // Initializing new distributions
        int* newDist = new int[this->bins];
        double* newNormalized = new double[this->bins];

        newDist = distribution(this->agents, this->N, this->maxMoney, this->dm);
        newNormalized = normalize(newDist, this->bins, this->N);
        diff = difference(this->dist, newNormalized, this->bins);

        // Updating distribution
        this->dist = newNormalized;

        // Freeing up memory
        delete [] newDist;
        delete [] newNormalized;
      }
    }
  }

  // Returning transactions needed to reach equilibrium
  return eq;
}

void NeighboursMarket::simulateMarket(int runs, double saving, double diffLimit, double alfa) {
  /*
  Function to simulate neighbours market in 'runs' runs with saving quotient
  'saving'
  */
  // Initializing transactions variable
  int n;

  // Looping over all runs
  for(int i = 0; i < runs; i++) {
    // Print to console to keep track during execution
    std::cout << "Run nr. " << i+1 << std::endl;

    // Performing transaction
    n = this->transactions(diffLimit, saving, alfa);

    // Updating meanDist
    for(int j = 0; j < this->bins; j++) {
      this->meanDist[j] += this->dist[j];
    }
  }

  // Taking mean
  for(int i = 0; i < this->bins; i++) {
    this->meanDist[i] /= runs;
  }
}

// FormeTransactions class
FormerTransactionsMarket::FormerTransactionsMarket(int N, double m0, double dm, double maxMoney)
  : Market(N, m0, dm, maxMoney)
  {
  // Allocate memory for former transactions matrix
  this->formerTransactions = new int*[this->N];
  for(int i = 0; i < this->N; i++) {
    this->formerTransactions[i] = new int[this->N];
  }
  // Initialize former transactions matrix with zeros
  for(int i = 0; i < this->N; i++) {
    for(int j = 0; j < this->N; j++) {
      this->formerTransactions[i][j] = 0;
    }
  }
}

int FormerTransactionsMarket::transactions(double diffLimit, double saving, double alfa, double gamma) {
  /*
  Function to perform transactions until equilibrium is achieved.
  Both the financial difference and number of previous transactions
  between two agents is accounted for
  */
  int agent1, agent2, c;
  double mi, mj, epsilon, dm, a, prob, delta;
  int eq = 0;
  int maxTransactions = 10000000;
  double diff = diffLimit + 1;

  // Seed for pseudo-random number generator
  std::srand(std::time(0));

  // Perform transactions as long as equilibrium is not achieved
  while(diff > diffLimit && eq < maxTransactions) {
    // Draw two random agents
    agent1 = std::rand() % this->N;
    agent2 = std::rand() % this->N;

    // Find previous transactions between agents
    c = this->formerTransactions[agent1][agent2] +
        this->formerTransactions[agent2][agent1];

    // Draw again if agent1 == agent2 until they are different
    while(agent1 == agent2) {
      agent2 = std::rand() % this->N;
    }
    mi = this->agents[agent1];
    mj = this->agents[agent2];

    // Draw random nubmer in unifrom distribution to decide whether
    // transaction happens
    a = (double) std::rand() / RAND_MAX;
    // Calculate probability of transaction to happen
    prob = pow(fabs(mi - mj),-alfa) * pow((c + 1),gamma);

    // Perform transaction if a <= prob
    if(a <= prob) {
      // Increment number of transactions between agents and transactions
      this->formerTransactions[agent1][agent2]++;
      eq++;

      // Draw random number epsilon from uniform distribution
      epsilon = (double) std::rand() / RAND_MAX;

      // Calculate transaction change
      delta = (1 - saving) * (epsilon*mj - (1 - epsilon)*mi);

      // Update income of agents and incrementing eq
      this->agents[agent1] += delta;
      this->agents[agent2] -= delta;

      // Finding new distribution and calculating variance between old and new
      // Only do this with a frequency of 100000 transactions
      if(eq % 100000 == 0) {
        // Initializing new distributions
        int* newDist = new int[this->bins];
        double* newNormalized = new double[this->bins];

        newDist = distribution(this->agents, this->N, this->maxMoney, this->dm);
        newNormalized = normalize(newDist, this->bins, this->N);
        diff = difference(this->dist, newNormalized, this->bins);

        // Updating distribution
        this->dist = newNormalized;

        // Freeing up memory
        delete [] newDist;
        delete [] newNormalized;
      }
    }
  }

  // Returning transactions needed to reach equilibrium
  return eq;
}

void FormerTransactionsMarket::simulateMarket(int runs, double saving, double diffLimit, double alfa, double gamma) {
  /*
  Function to simulate market with 'runs' runs
  */
  // Initialize transactions variable
  int n;

  // Looping over all runs
  for(int i = 0; i < runs; i++) {
    // Print to console to keep track
    std::cout << "Run nr. " << i+1 << std::endl;

    // Performing transactions
    n = this->transactions(diffLimit, saving, alfa, gamma);

    // Updating meanDist
    for(int j = 0; j < this->bins; j++) {
      this->meanDist[j] += this->dist[j];
    }
  }

  // Taking mean
  for(int i = 0; i < this->bins; i++) {
    this->meanDist[i] /= runs;
  }
}
