#ifndef _FOURTHORDERTENSOR_HPP_
#define _FOURTHORDERTENSOR_HPP_

#include <cassert>
#include <vector>

#include <base/tensor.h>
#include "Exception.hpp"

using namespace dealii;

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
    
    /**
     *  Set to be the inner product of another fourth order tensor and a matrix
     *  
     *  @param tensor A fourth order tensor
     *  @param matrix A Deal.II matrix
     *  @param component  The component in the fourth order tensor with which to sum
     *    (indexed from ZERO)
     * 
     *  ie. if component=0, X_{RM} T_{MNPQ} is returned  
     *  ie. if component=2, X_{RQ} T_{MNPQ} is returned
     * 
     */
    void SetAsProduct(FourthOrderTensor<DIM>& tensor, const Tensor<2,DIM>& matrix, unsigned component)
    {
        Zero();
        
        // messy repeated code but not sure how to do this neatly and efficiently..
        switch (component)
        {
            case 0:
            {
                for(unsigned M=0; M<DIM; M++)
                {
                    for(unsigned N=0; N<DIM; N++)
                    {
                        for(unsigned P=0; P<DIM; P++)
                        {
                            for(unsigned Q=0; Q<DIM; Q++)
                            {   
                                unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                                for(unsigned s=0; s<DIM; s++)
                                {
                                    mData[index] += matrix[M][s] * tensor(s,N,P,Q);
                                }
                            }
                        }
                    }
                }
                break;
            }
            case 1:
            {
                for(unsigned M=0; M<DIM; M++)
                {
                    for(unsigned N=0; N<DIM; N++)
                    {
                        for(unsigned P=0; P<DIM; P++)
                        {
                            for(unsigned Q=0; Q<DIM; Q++)
                            {   
                                unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                                for(unsigned s=0; s<DIM; s++)
                                {
                                    mData[index] += matrix[N][s] * tensor(M,s,P,Q);
                                }
                            }
                        }
                    }
                }
                break;
            }
            case 2:
            {
                for(unsigned M=0; M<DIM; M++)
                {
                    for(unsigned N=0; N<DIM; N++)
                    {
                        for(unsigned P=0; P<DIM; P++)
                        {
                            for(unsigned Q=0; Q<DIM; Q++)
                            {   
                                unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                                for(unsigned s=0; s<DIM; s++)
                                {
                                    mData[index] += matrix[P][s] * tensor(M,N,s,Q);
                                }
                            }
                        }
                    }
                }
                break;
            }
            case 3:
            {
                for(unsigned M=0; M<DIM; M++)
                {
                    for(unsigned N=0; N<DIM; N++)
                    {
                        for(unsigned P=0; P<DIM; P++)
                        {
                            for(unsigned Q=0; Q<DIM; Q++)
                            {   
                                unsigned index = M*mDimCubed + N*mDimSqd + P*DIM + Q;
                                for(unsigned s=0; s<DIM; s++)
                                {
                                    mData[index] += matrix[Q][s] * tensor(M,N,P,s);
                                }
                            }
                        }
                    }
                }
                break;
            }
            default:
            {
                EXCEPTION("Component not 0, 1, 2, or 3");
            }
        }
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
