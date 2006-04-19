#ifndef _TESTFOURTHORDERTENSOR_HPP_
#define _TESTFOURTHORDERTENSOR_HPP_

#include <cxxtest/TestSuite.h>
#include "FourthOrderTensor.hpp"

class TestFourthOrderTensor : public CxxTest::TestSuite
{
public:
	void testFourthOrderTensor()
	{
		int dim;
		
		dim=1;
		FourthOrderTensor<1> tensor1;
		
		TS_ASSERT_EQUALS( tensor1.rGetVal().size(), dim );		
		for(int M=0; M<dim; M++)
		{
			TS_ASSERT_EQUALS( tensor1.rGetVal()[M].size(), dim );		
			for(int N=0; N<dim; N++)
			{	
				TS_ASSERT_EQUALS( tensor1.rGetVal()[M][N].size(), dim );		
				for(int P=0; P<dim; P++)
				{
					TS_ASSERT_EQUALS( tensor1.rGetVal()[M][N][P].size(), dim );		
					for(int Q=0; Q<dim; Q++)
					{
						TS_ASSERT_DELTA( tensor1.rGetVal()[M][N][P][Q], 0, 1e-12 );		
					}
				}
			}
		}		


		dim=2;
		FourthOrderTensor<2> tensor2;
		
		TS_ASSERT_EQUALS( tensor2.rGetVal().size(), dim );		
		for(int M=0; M<dim; M++)
		{
			TS_ASSERT_EQUALS( tensor2.rGetVal()[M].size(), dim );		
			for(int N=0; N<dim; N++)
			{	
				TS_ASSERT_EQUALS( tensor2.rGetVal()[M][N].size(), dim );		
				for(int P=0; P<dim; P++)
				{
					TS_ASSERT_EQUALS( tensor2.rGetVal()[M][N][P].size(), dim );		
					for(int Q=0; Q<dim; Q++)
					{
						TS_ASSERT_DELTA( tensor2.rGetVal()[M][N][P][Q], 0, 1e-12 );		
					}
				}
			}
		}	



		dim=3;
		FourthOrderTensor<3> tensor3;
		
		TS_ASSERT_EQUALS( tensor3.rGetVal().size(), dim );		
		for(int M=0; M<dim; M++)
		{
			TS_ASSERT_EQUALS( tensor3.rGetVal()[M].size(), dim );		
			for(int N=0; N<dim; N++)
			{	
				TS_ASSERT_EQUALS( tensor3.rGetVal()[M][N].size(), dim );		
				for(int P=0; P<dim; P++)
				{
					TS_ASSERT_EQUALS( tensor3.rGetVal()[M][N][P].size(), dim );		
					for(int Q=0; Q<dim; Q++)
					{
						TS_ASSERT_DELTA( tensor3.rGetVal()[M][N][P][Q], 0, 1e-12 );		
					}
				}
			}
		}	
		
		
		tensor3.SetAsTensorProductOfIdentities();
		for(int M=0; M<dim; M++)
		{
			for(int N=0; N<dim; N++)
			{	
				for(int P=0; P<dim; P++)
				{
					for(int Q=0; Q<dim; Q++)
					{
						if((M==N) && (P==Q))
						{
							TS_ASSERT_DELTA( tensor3.rGetVal()[M][N][P][Q], 1.0, 1e-12 );		
						}
						else
						{
							TS_ASSERT_DELTA( tensor3.rGetVal()[M][N][P][Q], 0.0, 1e-12 );		
						}
					}
				}
			}
		}	
	}

	
};
#endif //_TESTFOURTHORDERTENSOR_HPP_
