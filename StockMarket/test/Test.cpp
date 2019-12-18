#include "pch.h"
#include "CppUnitTest.h"
#include "../Market/Market.h"
#include "../Market/tools.h"
#include "../Market/tools.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Test
{
	TEST_CLASS(TestTools)
	{
	public:
		
		TEST_METHOD(TestDistribution)
		{
			double* arr = new double[5];
			
			arr[0] = 0.5;
			arr[1] = 1.5;
			arr[2] = 2.5;
			arr[3] = 3.5;
			arr[4] = 3.5;

			int N = 5;
			double limit = 5.0;
			double dx = 1.0;

			double* trueDistrib = new double[5];

			trueDistrib[0] = 1.0;
			trueDistrib[1] = 1.0;
			trueDistrib[2] = 1.0;
			trueDistrib[3] = 2.0;
			trueDistrib[4] = 0.0;

			int* distrib = distribution(arr, N, limit, dx);

			bool success = true;

			for (int i = 0; i < limit; i++) {
				if (fabs(trueDistrib[i] - distrib[i]) > 0.0001) {
					success = false;
					break;
				}
			}

			Assert::IsTrue(success);

		}
		TEST_METHOD(TestNormalize)
		{

			int* arr = new int[5];

			arr[0] = 1.0;
			arr[1] = 2.0;
			arr[2] = 3.0;
			arr[3] = 4.0;
			arr[4] = 0.0;

			int length = 5;
			int n = 10;

			double* trueNorm = new double[5];

			trueNorm[0] = 0.1;
			trueNorm[1] = 0.2;
			trueNorm[2] = 0.3;
			trueNorm[3] = 0.4;
			trueNorm[4] = 0.0;

			double* norm = normalize(arr, length, n);

			bool success = true;

			for (int i = 0; i < length; i++) {
				if (fabs(trueNorm[i] - norm[i]) > 0.0001) {
					success = false;
					break;
				}
			}
			Assert::IsTrue(success);

		}

		TEST_METHOD(TestDifference)
		{

			double* arr1 = new double[5];
			double* arr2 = new double[5];

			arr1[0] = 0.1;
			arr1[1] = 0.2;
			arr1[2] = 0.3;
			arr1[3] = 0.4;
			arr1[4] = 0.0;

			arr2[0] = 0.1;
			arr2[1] = 0.1;
			arr2[2] = 0.5;
			arr2[3] = 0.2;
			arr2[4] = 0.2;

			int N = 5;

			double trueDiff = 0.7/N;

			double diff = difference(arr1, arr2, N);

			bool success = true;


			if (fabs(trueDiff - diff) > 0.0001) {
				success = false;
			}
			Assert::IsTrue(success);
		}

	};
}
