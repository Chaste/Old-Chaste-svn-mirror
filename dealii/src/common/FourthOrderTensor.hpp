#ifndef _FOURTHORDERTENSOR_HPP_
#define _FOURTHORDERTENSOR_HPP_

#include <cassert>
#include <vector>

/** 
 *  FourthOrderTensor
 * 
 *  A class of fourth order tensors (ie tensors with four indices), over arbitrary dimension
 * 
 */ 
template <int DIM> 
class FourthOrderTensor
{
private:
    std::vector<double> mData;
    unsigned mDimSqd;
    unsigned mDimCubed;
    unsigned mDimToFour;
    
public:
    FourthOrderTensor()
    {
        // check dim>0 but <4
        assert(DIM>0);
        assert(DIM<4);
        
        mDimSqd = DIM*DIM;
        mDimCubed = DIM*DIM*DIM;
        mDimToFour = DIM*DIM*DIM*DIM;
        
        // allocate memory and zero entries
        mData.resize(mDimToFour, 0.0);
    }
    
    double& operator()(unsigned M, unsigned N, unsigned P, unsigned Q) 
    {
        assert(M<DIM);
        assert(N<DIM);
        assert(P<DIM);
        assert(Q<DIM);
        
        unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
        return mData[index];
    }
    
    void Zero()
    {
        for(unsigned i=0; i<mDimToFour; i++)
        {
            mData[i] = 0.0;
        }
    }
};

#endif //_FOURTHORDERTENSOR_HPP_
