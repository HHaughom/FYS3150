#pragma once
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long* idum) {
	/*
	 ** The function
	 **           ran0()
	 ** is an "Minimal" random number generator of Park and Miller
	 ** (see Numerical recipe page 279). Set or reset the input value
	 ** idum to any integer value (except the unlikely value MASK)
	 ** to initialize the sequence; idum must not be altered between
	 ** calls for sucessive deviates in a sequence.
	 ** The function returns a uniform deviate between 0.0 and 1.0.
	 */
	long     k;
	double   ans;

	*idum ^= MASK;
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0)* idum += IM;
	ans = AM * (*idum);
	*idum ^= MASK;
	return ans;
}

double gammln(double xx)
{	/*Support function for Gauss-Laguerre method*/
	double x, y, tmp, ser;
	static double cof[6] = { 76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5 };
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005 * ser / x);
}
