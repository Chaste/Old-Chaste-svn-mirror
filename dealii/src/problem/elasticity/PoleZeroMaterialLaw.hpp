#ifndef POLEZEROMATERIALLAW_HPP_
#define POLEZEROMATERIALLAW_HPP_


#include "AbstractIncompressibleMaterialLaw.hpp"
#include "Exception.hpp"



/**
 *  Pole-zero material law
 *
 *  NOTE: CHANGED THE CODE AT THE MOMENT TO NOT USE THE POSITIVE PART BIT BELOW
 *  AS THEN IT IS NOT TWICE DIFFERENTIABLE
 * 
 *  W = \Sum_{M,N=1..3} k_{MN} ([ E_{MN} ]_+)^2 / (a_{MN} - E_{MN})^b_{MN}
 * 
 *  Note the positive part operator in the numerator, so that
 *  the term in W corresponding to M,N is zero if E_MN < 0. (This differs
 *  from the original pole-zero paper but seems to be what they meant..)
 * 
 *  Note that is the parameters k4,k5,k6,a4,a5,a6 etc are known, then
 *  k01=k10=0.5*k4 and similarly with k5,k6, but a01=a10=a4 etc.
 * 
 *  Not isotropic, so inherits directly from AbstractIncompressibleMaterialLaw
 */
template<unsigned DIM>
class PoleZeroMaterialLaw : public AbstractIncompressibleMaterialLaw<DIM>
{
friend class TestMaterialLaws;

private :
    std::vector<std::vector<double> > mK;
    std::vector<std::vector<double> > mA;
    std::vector<std::vector<double> > mB;
    
    Tensor<2,DIM> mIdentity;

protected :
    /** 
     *  Protected default constructor doing nothing. Just saw inherited classes
     *  can be instantiated and THEN set up the parameters 
     */  
    PoleZeroMaterialLaw()
    {
    }
    
    /** 
     *  Set k, a, and b. To be called by the constuctor or a child class 
     *  Set comments for constructor.
     */
    void SetParameters(std::vector<std::vector<double> > k,
                       std::vector<std::vector<double> > a,
                       std::vector<std::vector<double> > b)
    {
        if (DIM!=2 && DIM !=3)
        {
            EXCEPTION("Can only have a 2 or 3d incompressible pole-zero law");
        }
        
        assert(k.size()==DIM);
        assert(a.size()==DIM);
        assert(b.size()==DIM);
 
        for(unsigned i=0; i<DIM; i++)
        {
            assert(k[i].size()==DIM);
            assert(a[i].size()==DIM);
            assert(b[i].size()==DIM);

            for(unsigned j=0; j<DIM; j++)
            {
                assert( k[i][j] = k[j][i] );
                assert( a[i][j] = a[j][i] );
                assert( b[i][j] = b[j][i] );
            }
        }
        
        mK = k;
        mA = a;
        mB = b;

        for(unsigned M=0; M<DIM; M++)
        {
            for(unsigned N=0; N<DIM; N++)
            {
                mIdentity[M][N] = M==N ? 1.0 : 0.0;
            }
        }
    }
    
public :
    /**
     *  Constructor, taking in parameters k_i, a_i, b_i as matrices.
     *  These matrices must be of size DIM-by-DIM and must be symmetric
     * 
     *  Note: using the k_1..k_6 convention,  k_4 = 2*k[0][1] = 2*k[1][0], etc
     */
     PoleZeroMaterialLaw(std::vector<std::vector<double> > k,
                         std::vector<std::vector<double> > a,
                         std::vector<std::vector<double> > b)
    {
        SetParameters(k,a,b);
    }
    
    void ComputeStressAndStressDerivative(Tensor<2,DIM>&          C,
                                          Tensor<2,DIM>&          invC,
                                          double                  pressure,
                                          SymmetricTensor<2,DIM>& T,
                                          FourthOrderTensor<DIM>& dTdE,
                                          bool                    computeDTdE)
    {
        assert(fabs(C[0][1]-C[1][0]) < 1e-6);
        
        Tensor<2,DIM> E = 0.5*(C-mIdentity);
        
        for(unsigned M=0; M<DIM; M++)
        {
            for(unsigned N=0; N<DIM; N++)
            {
                double e = E[M][N];
              //  if(e > 0)
                {
                    double b = mB[M][N];
                    double a = mA[M][N];
                    double k = mK[M][N];
                    
                    //if this fails one of the strain values got too large for the law
                    assert(e < a);
                    
                    T[M][N] =   k  
                              * e 
                              * (2*(a-e) + b*e)
                              * pow(a-e,-b-1)
                              - pressure*invC[M][N];
                }
//                else
//                {
//                    T[M][N] = 0.0;
//                }
            }
        }

        if(computeDTdE)
        {
            for(unsigned M=0; M<DIM; M++)
            {
                for(unsigned N=0; N<DIM; N++)
                {
                    for(unsigned P=0; P<DIM; P++)
                    {
                        for(unsigned Q=0; Q<DIM; Q++)
                        {
                            dTdE(M,N,P,Q) = 2 * pressure * invC[M][P] * invC[Q][N];
                        }
                    }
                    
                    double e = E[M][N];
                 //   if(e > 0)
                    {
                        double b = mB[M][N];
                        double a = mA[M][N];
                        double k = mK[M][N];
    
                        dTdE(M,N,M,N) +=   k 
                                         * pow(a-e, -b-2)
                                         * (   
                                              2*(a-e)*(a-e)
                                            + 4*b*e*(a-e)
                                            + b*(b+1)*e*e
                                           );
                    }
                }
            }
        }
    }
    
    double GetZeroStrainPressure()
    {
        return 0.0;
    }
    
    /** Scale the dimensional material parameters (ie the K's) */
    void ScaleMaterialParameters(double scaleFactor)
    {
        assert(scaleFactor > 0.0);
        for(unsigned i=0; i<mK.size(); i++)
        {
            for(unsigned j=0; j<mK[i].size(); j++)
            {
                mK[i][j] /= scaleFactor;
            }
        }
    }    
};


#endif /*POLEZEROMATERIALLAW_HPP_*/
