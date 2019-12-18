// Tools/helping functions used

#include <cmath>
#include "../include/tools.h"

// For debugging
#include <iostream>

double* normalize(int* dist, int length, int n) {
  /*
  Function to normalize a distribution with length
  'length' with factor n
  */
  // Initialize result distribution
  double* result = new double[length];
  for(int i = 0; i < length; i++) {
    result[i] = (double)dist[i] / n;
  }

  // Returning normalized distribution
  return result;
}

int* distribution(double* array, int N, double limit, double dx) {
  /*
  Function that returns the distribution (unnormalized) of array
  */
  // Initialize channels and distribution array
  int bins = (int) limit/dx;
  int* dist = new int[bins];

  // Fill dist with data from array
  int index;
  for(int i = 0; i < N; i++) {
    if(array[i] > limit) {
      // If this happens often, it would affect the results
      // and it would show in the plot as a spike at the last point
      // To fix it, increase the limit
      index = bins-1;
    }
    else {
      index = (int) floor(std::fmod(array[i],limit) / dx);
    }
    dist[index]++;
  }

  // Return distribution
  return dist;
}

double difference(double* array1, double* array2, int n) {
  /*
  Fucntion to calculate the total difference between 2
  normalized distribution arrays of length 'n'
  */
  // Initialize variance variable and stepsize
  double diff = 0;

  // Calculating the variance
  for(int i = 0; i < n; i++) {
    diff += fabs(array1[i] - array2[i]);
  }

  // Returning mean difference
  return diff/n;
}
