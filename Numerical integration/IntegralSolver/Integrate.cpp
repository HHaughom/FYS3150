#include "Integrate.h"
#include "lib.h"

//INTEGRAND

double integrandCartesian(double r1[], double r2[]) {
	/*Integrand given in cartesian coordinates. 
	  Returns integrand value for given electron positions r1 and r2.
	  Returns 0 if denominator close to zero.*/


	double functionValue;
	double alpha = 2;

	// evaluate the different terms of the exponential
	double R1 = sqrt(pow(r1[0], 2) + pow(r1[1], 2) + pow(r1[2], 2));
	double R2 = sqrt(pow(r2[0], 2) + pow(r2[1], 2) + pow(r2[2], 2));
	double distance = sqrt(pow(r1[0] - r2[0], 2) + pow(r1[1] - r2[1], 2) + pow(r1[2] - r2[2], 2));

	if (distance < ZERO) {
		functionValue = 0;
	}
	else {
		functionValue = exp(-2 * alpha * (R1 + R2)) / distance;
	}
	return functionValue;
}

double integrandSpherical(double r[], double theta[], double phi[]) {
	/*Integrand given in spherical coordinates.
	  Returns integrand value for given electron positions (ri, thetai, psii).
	  Returns 0 if denominator close to zero.*/

	double functionValue, enumerator, denominator;
	double alpha = 2;
	double cosBeta = cos(theta[0]) * cos(theta[1]) + sin(theta[0]) * sin(theta[1]) * cos(phi[0] - phi[1]);

	enumerator = sin(theta[0]) * sin(theta[1]) * pow(r[0], 2) * pow(r[1], 2) * exp(-4 * (r[0] + r[1]));
	denominator = sqrt(r[0] * r[0] + r[1] * r[1] - 2 * r[0] * r[1] * cosBeta);

	//handeling of singular points
	if (denominator < ZERO || isnan(denominator)) {
		functionValue = 0;
	}
	else {
		functionValue = enumerator / denominator;
	}
	return functionValue;
}

//GAUSS LEGENDRE 

void pointsGaussLegendre(double a, double b, double* x, double* w, int N) {
	/*Function copied from repository
	
		github.com/CompPhysics/ComputationalPhysics/blob/master/doc/Projects/2019/Project3/CodeExamples/exampleprogram.cpp

	  Finds the N meshpoints x and weights w from the legendre polynomials,
	  given integral limits [a,b]. Uses that integrand is symmetric.
		**Calculates the N Legendre polynomials through the recurrance relation.
		**Meshpoints found from zeros of polynomials, using Newtons method
		**Weights found from inverse of polynomials*/


	int m;
	double z1, z, xm, xl, pp, p3, p2, p1;
	double* x_low, * x_high, * w_low, * w_high;

	m = (N + 1) / 2;

	xm = 0.5 * (b + a);
	xl = 0.5 * (b - a);

	x_low = x;                          
	x_high = x + N - 1;
	w_low = w;
	w_high = w + N - 1;

	for (int i = 1; i <= m; i++) {                             // loops over desired roots
		z = cos(PI * (i - 0.25) / (N + 0.5));

		/*
	** Starting with the above approximation to the ith root
		** we enter the mani loop of refinement bt Newtons method.
		*/

		do {
			p1 = 1.0;
			p2 = 0.0;

			/* ** loop up recurrence relation to get the
			   ** Legendre polynomial evaluated at x*/

			for (int j = 1; j <= N; j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
			}

			/*  ** p1 is now the desired Legrendre polynomial. Next compute
				** ppp its derivative by standard relation involving also p2,
				** polynomial of one lower order.*/

			pp = N * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp;                   // Newton's method
		} while (fabs(z - z1) > ZERO);

		/*
		** Scale the root to the desired interval and put in its symmetric
		** counterpart. Compute the weight and its symmetric counterpart
		*/

		*(x_low++) = xm - xl * z;
		*(x_high--) = xm + xl * z;
		*w_low = 2.0 * xl / ((1.0 - z * z) * pp * pp);
		*(w_high--) = *(w_low++);
	}

}

double gaussLegendre(int N) {
	/*Function closely based on gauleg() function in

		github.com/CompPhysics/ComputationalPhysics/blob/master/doc/Projects/2019/Project3/CodeExamples/exampleprogram.cpp

	  Solves the integral using the cartesian integrand and
	  meshpoints and weights found in pointsGaussLegendre()*/

	double r1[3], r2[3], wgl[6];

	double* x = new double[N];
	double* w = new double[N];
	
	double integral = 0;
	double a = -INF;
	double b = INF;

	Timer timer;

	pointsGaussLegendre(a, b, x, w, N);

	for (int i = 0; i < N; i++) {

		r1[0] = x[i];
		wgl[0] = w[i];

		for (int j = 0; j < N; j++) {

			r1[1] = x[j];
			wgl[1] = w[j];

			for (int k = 0; k < N; k++) {

				r1[2] = x[k];
				wgl[2] = w[k];

				for (int l = 0; l < N; l++) {

					r2[0] = x[l];
					wgl[3] = w[l];

					for (int m = 0; m < N; m++) {

						r2[1] = x[m];
						wgl[4] = w[m];

						for (int n = 0; n < N; n++) {

							r2[2] = x[n];
							wgl[5] = w[n];

							integral += wgl[0] * wgl[1] * wgl[2] * wgl[3] * wgl[4] * wgl[5] * integrandCartesian(r1, r2);

						}
					}
				}
			}
		}
	}

	delete[] x;
	delete[] w;

	return integral;
}
 
//GAUSS LAGUERRE

void pointsGaussLaguerre(double* x, double* w, int N, double alpha) {
	/*Function copied from repository

		github.com/CompPhysics/ComputationalPhysics/blob/master/doc/Projects/2019/Project3/CodeExamples/exampleprogram.cpp

	  Finds the N meshpoints x and weights w from the laguerre polynomials,
	  given parameter alpha.
		**Calculates the N Laguerre polynomials through the recurrance relation.
		**Meshpoints found from zeros of polynomials, using Newtons method
		**Weights found from inverse of polynomials
		**MAXITERATION limits in case of non-convergence*/

	int its;
	double ai, p1, p2, p3, pp, z, z1;

	for (int i = 1; i <= N; i++) {
		if (i == 1) {
			z = (1.0 + alpha) * (3.0 + 0.92 * alpha) / (1.0 + 2.4 * N + 1.8 * alpha);
		}
		else if (i == 2) {
			z += (15.0 + 6.25 * alpha) / (1.0 + 0.9 * alpha + 2.5 * N);
		}
		else {
			ai = i - 2.0;
			z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alpha /
				(1.0 + 3.5 * ai)) * (z - x[i - 2]) / (1.0 + 0.3 * alpha);
		}
		for (its = 1; its <= MAXIT; its++) {
			p1 = 1.0;
			p2 = 0.0;
			for (int j = 1; j <= N; j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j - 1.0 + alpha - z) * p2 - (j - 1.0 + alpha) * p3) / j;
			}
			pp = (N * p1 - (N + alpha) * p2) / z;
			z1 = z;
			z = z1 - p1 / pp;
			if (fabs(z - z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations" << endl;
		x[i] = z;
		w[i] = -exp(gammln(alpha + N) - gammln((double)N)) / (pp * N * p2);
	}
}

double gaussLaguerre(int N) {
	/*Function closely based on gaussLaguerre() in

		github.com/CompPhysics/ComputationalPhysics/blob/master/doc/Projects/2019/Project3/CodeExamples/exampleprogram.cpp

	  Solves the integral using the spherical integrand and
	  meshpoints and weights found in pointsGaussLaguerre()*/

	double r[2], theta[2], phi[2], w[6], W;

	double* R = new double[N + 1];
	double* wgla = new double[N + 1];

	double* THETA = new double[N];
	double* wgle1 = new double[N];

	double* PHI = new double[N];
	double* wgle2 = new double[N];

	double integral = 0;
	double alpha = 2.0;

	double theta1 = 0;
	double theta2 = PI;

	double phi1 = 0;
	double phi2 = 2 * PI;

	Timer timer;

	pointsGaussLaguerre(R, wgla, N, alpha);
	pointsGaussLegendre(theta1, theta2, THETA, wgle1, N);
	pointsGaussLegendre(phi1, phi2, PHI, wgle2, N);

	for (int i = 1; i < N + 1; i++) {

		r[0] = R[i];
		w[0] = wgla[i];

		for (int j = 1; j < N + 1; j++) {

			r[1] = R[j];
			w[1] = wgla[j];

			W = pow(r[0], 2) * pow(r[1], 2) * exp(-(r[0] + r[1]));

			for (int k = 0; k < N; k++) {

				theta[0] = THETA[k];
				w[2] = wgle1[k];

				for (int l = 0; l < N; l++) {

					theta[1] = THETA[l];
					w[3] = wgle1[l];

					for (int m = 0; m < N; m++) {

						phi[0] = PHI[m];
						w[4] = wgle2[m];

						for (int n = 0; n < N; n++) {

							phi[1] = PHI[n];
							w[5] = wgle2[n];

							integral += w[0] * w[1] * w[2] * w[3] * w[4] * w[5] * integrandSpherical(r, theta, phi)/W;
						}
					}
				}
			}
		}
	}

	delete[] R;
	delete[] wgla;
	delete[] THETA;
	delete[] wgle1;
	delete[] PHI;
	delete[] wgle2;

	return integral;
}

//GAUSS LAGUERRE PARALLEL

double gaussLaguerrePara(int N) {
	/*Function closely based on gaussLaguerre() in
		github.com/CompPhysics/ComputationalPhysics/blob/master/doc/Projects/2019/Project3/CodeExamples/exampleprogram.cpp
	  Solves the integral using the spherical integrand and
	  meshpoints and weights found in pointsGaussLaguerre().

	  Includes parallelization, but is otherwise identical to gaussLaguerre*/


	double r[2], theta[2], phi[2], w[6], W;

	double* R = new double[N + 1];
	double* wgla = new double[N + 1];

	double* THETA = new double[N];
	double* wgle1 = new double[N];

	double* PHI = new double[N];
	double* wgle2 = new double[N];

	double integral = 0;
	double alpha = 2.0;

	double theta1 = 0;
	double theta2 = PI;

	double phi1 = 0;
	double phi2 = 2 * PI;

	Timer timer;

	pointsGaussLaguerre(R, wgla, N, alpha);
	pointsGaussLegendre(theta1, theta2, THETA, wgle1, N);
	pointsGaussLegendre(phi1, phi2, PHI, wgle2, N);

	int j, k, l, m, n;
	#pragma omp parallel for default(shared) private (j,k,l,m,n) reduction(+:integral)
	for (int i = 1; i < N + 1; i++) {

		r[0] = R[i];
		w[0] = wgla[i];

		for ( j = 1; j < N + 1; j++) {

			r[1] = R[j];
			w[1] = wgla[j];

			W = pow(r[0], 2) * pow(r[1], 2) * exp(-(r[0] + r[1]));

			for (k = 0; k < N; k++) {

				theta[0] = THETA[k];
				w[2] = wgle1[k];

				for ( l = 0; l < N; l++) {

					theta[1] = THETA[l];
					w[3] = wgle1[l];

					for ( m = 0; m < N; m++) {

						phi[0] = PHI[m];
						w[4] = wgle2[m];

						for (n = 0; n < N; n++) {

							phi[1] = PHI[n];
							w[5] = wgle2[n];

							integral += w[0] * w[1] * w[2] * w[3] * w[4] * w[5] * integrandSpherical(r, theta, phi) / W;
						}
					}
				}
			}
		}
	}

	delete[] R;
	delete[] wgla;
	delete[] THETA;
	delete[] wgle1;
	delete[] PHI;
	delete[] wgle2;

	return integral;
}

//BRUTE MONTE CARLO

double* bruteMonteCarlo(int N){
	/*Function closely based on example 

		github.com/CompPhysics/ComputationalPhysics/blob/master/doc/Projects/2019/Project3/CodeExamples/exampleprogram.cpp

	  Solves the integral using the spherical integrand and
	  meshpoints and weights found in pointsGaussLaguerre()*/

	double r1[3], r2[3], f, res[2];
	double variance;
	double* result;

	double integral = 0.;  
	double fsquared = 0.; 

	double lambda = INF;
	long idum = -1;
	double PDF = 1.0/pow((2 * lambda), 6);

	Timer timer;

	for (int i = 1; i <= N; i++) {
		for (int j = 0; j < 3; j++) {
			r1[j] = -lambda + 2 * lambda * ran0(&idum);
			r2[j] = -lambda + 2 * lambda * ran0(&idum);
		}
		f = integrandCartesian(r1, r2)/PDF;
		integral += f;
		fsquared += pow(f, 2); 
	}

	integral /= ((double)N); 
	fsquared /= ((double)N); 
	variance = (fsquared - pow(integral, 2))/((double)N);

	res[0] = integral;
	res[1] = variance;
	result = res;

	return result;
}

//IMPORTANCE SAMPLING MONTE CARLO

double* impsampMonteCarlo(int N){

	double f, variance, res[2];
	double r[2], theta[2], phi[2];
	double* result;

	double fsquared = 0;
	double integral = 0.;  
	long idum = -1;

	double lambda = INF;
	double thetaMax = PI;
	double phiMax = 2* PI;
	double PDF;

	Timer timer;

	//evaluate the integral with importance sampling    
	for (int i = 1; i <= N; i++) {
		for (int j = 0; j < 2; j++) {
			r[j] = -log(1 - ran0(&idum))/4; 
			theta[j] = thetaMax * ran0(&idum); 
			phi[j] = phiMax * ran0(&idum); 
		}
		PDF = 16*exp(-4*r[0]) * exp(-4*r[1]) * pow(phiMax, -2) * pow(thetaMax, -2);
		f = integrandSpherical(r, theta, phi)/PDF;
		integral += f;
		fsquared += pow(f, 2);
	}

	integral /= ((double)N);
	fsquared /= ((double)N);
	variance = (fsquared - pow(integral,2))/((double)N);

	res[0] = integral;
	res[1] = variance;
	result = res;

	return result;
}

