#include "pch.h"
#include "CppUnitTest.h"
#include "../P3/Integrate.h"
#include "../P3/Integrate.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest {

	TEST_CLASS(Integrand) {
		/*
		Testing functions returning integrand values.
		*/
public:
	TEST_METHOD(integrandCart) {
		// Initialization
		double valCalc1, valExp1, valCalc2, valExp2;

		double r1[] = { 1.,2.,3. };
		double r2[] = { 4.,5.,6. };
		double r3[] = { 4.,5.,6. };

		// Calculate functionvalue 
		valCalc1 = integrandCartesian(r1, r2);
		valCalc2 = integrandCartesian(r2, r3);

		// Expecteded functionvalue
		valExp1 = 3.4731140393369526E-23;
		valExp2 = 0;

		// Checking for equality
		bool success = true;
		if ((abs(valCalc1 - valExp1) > EPS) || (abs(valCalc2 - valExp2) > EPS)) {
			success = false;
		}
		// Assert result
		Assert::IsTrue(success);
	}

	TEST_METHOD(integrandSpher) {
		// Initialization
		double valCalc1, valExp1, valCalc2, valExp2;

		double r1[] = { 2.,4. };
		double theta1[] = { 1.,2. };
		double phi1[] = { 3.,4. };

		double r2[] = { 2.,2. };
		double theta2[] = { PI / 2.0,PI / 2.0 };
		double phi2[] = { PI,PI };

		// Calculate functionvalue 
		valCalc1 = integrandSpherical(r1, theta1, phi1);
		valCalc2 = integrandSpherical(r2, theta2, phi2);

		// Expecteded functionvalue
		valExp1 = 3.85329853514885E-10;
		valExp2 = 0;

		// Checking for equality
		bool success = true;
		if ((abs(valCalc1 - valExp1) > EPS * 1E4) || (abs(valCalc2 - valExp2) > EPS)) {
			success = false;
		}
		// Assert result
		Assert::IsTrue(success);
	}

	};

	TEST_CLASS(pointFinder) {
		/*
		Testing function for finding mesh points and weights.
		*/
public:
	TEST_METHOD(gaussLeg) {
		// Initialization

		double xExp, wExp;

		double a = -1;
		double b = 1;
		int N = 2;
		
		double* x = new double[N];
		double* w = new double[N];

		// Calculate points and weights 

		pointsGaussLegendre(a, b, x, w, N);

		// Expecteded points and weights

		xExp = 1.0 / sqrt(3);
		wExp = 1;

		// Checking for equality
		bool success = true;
		if ((abs(xExp - abs(x[0])) > 1E-3) || (abs(wExp - w[0]) > 1E-3)) {
			success = false;
		}

		// Assert result
		Assert::IsTrue(success);
	}
	};

}
