#ifndef _FOURTHORDERTENSOR_HPP_
#define _FOURTHORDERTENSOR_HPP_

#include <cassert>
#include <vector>

/** FourthOrderTensor
 * 
 *  A class of fourth order tensors (ie tensors with four indices, eg. dTdE[M][N][P][Q]),
 *  over arbitrary dimension
 * 
 *  \todo overload [][][][] or use access via (,,,)
 *  \todo overload matrix multiplication
 * 
 */ 
template <int SPACE_DIM> 
class FourthOrderTensor
{
	
public:
	std::vector< std::vector < std::vector< std::vector<double> > > >   mVal;
	
	FourthOrderTensor()
	{
		// check dim>0 but <4, we don't want to create fourth order tensors over hyperspace
		assert(SPACE_DIM>0);
		assert(SPACE_DIM<4);
		
				
		// allocate memory and zero entries
		mVal.resize(SPACE_DIM);		
		for(int M=0; M<SPACE_DIM; M++)
		{
			mVal[M].resize(SPACE_DIM);		
			for(int N=0; N<SPACE_DIM; N++)
			{	
				mVal[M][N].resize(SPACE_DIM);		
				for(int P=0; P<SPACE_DIM; P++)
				{
					mVal[M][N][P].resize(SPACE_DIM);		
					for(int Q=0; Q<SPACE_DIM; Q++)
					{
						mVal[M][N][P][Q] = 0;
					}
				}
			}
		}				
	}
	
	
	void SetAsTensorProductOfIdentities()
	{
		for(int M=0; M<SPACE_DIM; M++)
		{
			for(int P=0; P<SPACE_DIM; P++)
			{	
				mVal[M][M][P][P] = 1;
			}
		}			
	}	
	
};



				

#endif //_FOURTHORDERTENSOR_HPP_
