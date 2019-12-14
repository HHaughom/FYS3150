#include "pch.h"
#include "CppUnitTest.h"
#include "..\Ising\Ising.h"
#include "..\Ising\ExpectationValues.h"
#include "..\Ising\Ising.cpp"

//#include "\Users\formIdabel\source\repos\Ising\Ising\Ising.cpp"
//#include <armadillo>
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest
{
	TEST_CLASS(isingTest)
	{
	public:
		
		TEST_METHOD(ExpectationVal)

		{
			double analE, analMabs, analCv, analChi;
			double compE, compMabs, compCv, compChi;
			double diffE, diffM, diffCv, diffChi; 
			bool success = false;

			Ising test = Ising(2, true);
			ExpectationValues ev = test.Metropolis(1.0, 1000000);
			double epsilon = 1e-3;

			analE = -1.99598;
			analMabs = 0.998661;
			analCv = 0.0320823;
			analChi = 0.00401074;

			compE = ev.expecE;
			compMabs = ev.expecMabs;
			compCv = ev.expecCv;
			compChi = ev.expecChiAbs/4.0;
			
			diffE = fabs(analE - compE);
			diffM = fabs(analMabs - compMabs);
			diffCv = fabs(analCv - compCv);
			diffChi = fabs(analChi - compChi);

			if ((diffE < epsilon) || (diffE < epsilon) || (diffE < epsilon) || (diffE < epsilon)) {
				success = true;
			}

			Assert::IsTrue(success);
		}

		TEST_METHOD(StateVal) 
		
		{
			double expecE, expecM, compE, compM;
			double diffE, diffM;
			bool success = false;

			Ising test = Ising(2, true);
			ExpectationValues ev = test.Metropolis(1.0, 1000000);
			double epsilon = 1e-10;

			compE = test.calcE();
			compM = test.calcM();

			expecM = 4.0;
			expecE = -4.0;

			diffE = fabs(compE - expecE);
			diffM = fabs(compM - expecM);

			if ((diffE < epsilon) || (diffM < epsilon)) {
				success = true;
			}

			Assert::IsTrue(success);

		}

	};
}
