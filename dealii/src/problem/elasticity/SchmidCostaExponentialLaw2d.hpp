#ifndef SCHMIDCOSTAEXPONENTIALLAW2D_HPP_
#define SCHMIDCOSTAEXPONENTIALLAW2D_HPP_

#include "AbstractIncompressibleMaterialLaw.hpp"

/**
 *  A 2d version of the material law in Costa, Holmes, McCulloch "Modelling Cardiac
 *  Mechanical Properties in Three Dimensions" Philo. Trans. R. Soc.
 * 
 *  W = aexp(Q-1)
 *  where Q = bff*Eff^2 + bfs*Efs^2 + bsf*Esf^2 + bss*Ess^2
 * 
 *  where the parameters are taken from the fitting in Schmid,Nash,Young,Hunter
 *  "Myocardial Material Parameter Estimation - A Comparative Study for Simple
 *  Shear" Transactions of the ASME.
 */
class SchmidCostaExponentialLaw2d : public AbstractIncompressibleMaterialLaw<2>
{

private :
    double mA;         // should be kPa as the assembler gets the active tension in kPa 
    std::vector<std::vector<double> > mB;
    Tensor<2,2> mIdentity;
    
public :
    SchmidCostaExponentialLaw2d()
    {
        mA = 0.221;    // kiloPascals, presumably, althought the paper doesn't say.
                       // gives results matching Pole-zero anyway.

        double bff = 42.5; // dimensionless
        double bfs = 11.0; // dimensionless
        double bss = 18.6; // dimensionless
        
        mB.resize(2);
        mB[0].resize(2);
        mB[1].resize(2);
        
        mB[0][0] = bff;
        mB[0][1] = bfs;
        mB[1][0] = bfs;
        mB[1][1] = bss;
        
        for(unsigned M=0; M<2; M++)
        {
            for(unsigned N=0; N<2; N++)
            {
                mIdentity[M][N] = M==N ? 1.0 : 0.0;
            }
        }
    }

    void ComputeStressAndStressDerivative(Tensor<2,2>&          C,
                                          Tensor<2,2>&          invC,
                                          double                pressure,
                                          SymmetricTensor<2,2>& T,
                                          FourthOrderTensor<2>& dTdE,
                                          bool                  computeDTdE)
    {
        Tensor<2,2> E = 0.5*(C-mIdentity);
        
        double Q = 0;
        for(unsigned M=0; M<2; M++)
        {
            for(unsigned N=0; N<2; N++)
            {
                Q += mB[M][N]*E[M][N]*E[M][N];
            }
        }
                
        double multiplier = mA*exp(Q-1);
        dTdE.Zero();

        for(unsigned M=0; M<2; M++)
        {
            for(unsigned N=0; N<2; N++)
            {
                T[M][N] = multiplier*mB[M][N]*E[M][N] - pressure*invC[M][N];
            
                if(computeDTdE)
                {
                    //dTdE(M,N,M,N) = multiplier * mB[M][N];
                    for(unsigned P=0; P<2; P++)
                    {
                        for(unsigned Q=0; Q<2; Q++)
                        {
                            dTdE(M,N,P,Q) =    multiplier * mB[M][N] * (M==P)*(N==Q)
                                            +  2*multiplier*mB[M][N]*mB[P][Q]*E[M][N]*E[P][Q]
                                            +  2*pressure*invC[M][P]*invC[Q][N];
                        }
                    }
                }
            }
        }
    }
    
    double GetA()
    {
        return mA;
    }

    std::vector<std::vector<double> > GetB()
    {
        return mB;
    }
    
    double GetZeroStrainPressure()
    {
        return 0.0;
    }
};

#endif /* SCHMIDCOSTAEXPONENTIALLAW2D_HPP_*/
