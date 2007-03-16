#ifndef TESTMATERIALLAWS_HPP_
#define TESTMATERIALLAWS_HPP_

#include <cxxtest/TestSuite.h>

#include "MooneyRivlinMaterialLaw.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include <cassert>


class TestMaterialLaws : public CxxTest::TestSuite
{
public:
    void testMooneyRivlinLaw()
    {
        TS_ASSERT_THROWS_ANYTHING(MooneyRivlinMaterialLaw<2> bad_mr_law(-3.0));
        TS_ASSERT_THROWS_ANYTHING(MooneyRivlinMaterialLaw<3> bad_mr_law2(3.0));
        TS_ASSERT_THROWS_ANYTHING(MooneyRivlinMaterialLaw<1> bad_mr_law3( 3.0, 1.0));
        
        double c1 = 2.0;
        
        MooneyRivlinMaterialLaw<2> ml_law_2d(c1);
        
        TS_ASSERT_DELTA(ml_law_2d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(ml_law_2d.Get_dW_dI1(1.0,0.0), c1, 1e-12);
        TS_ASSERT_DELTA(ml_law_2d.Get_d2W_dI1(1.0,0.0), 0.0, 1e-12);
        
        // compute the stress given C=delta_{MN} and p=zero_strain_pressure,
        // obviously it should be zero
        TS_ASSERT_DELTA(ml_law_2d.GetZeroStrainPressure(), 2*c1, 1e-12);
        Tensor<2,2> identity_strain_2d;
        for(unsigned i=0; i<2; i++)
        {
            for(unsigned j=0; j<2; j++)
            {
                identity_strain_2d[i][j]=0;
            }
            identity_strain_2d[i][i]=1;
        }

        SymmetricTensor<2,2> T_2d;
        ml_law_2d.Compute2ndPiolaKirchoffStress(identity_strain_2d,
                                                ml_law_2d.GetZeroStrainPressure(),
                                                T_2d);
        for(unsigned i=0; i<2; i++)
        {
            for(unsigned j=0; j<2; j++)
            {
                TS_ASSERT_DELTA(T_2d[i][j],0.0,1e-12);
            }
        }
                
                
        double c2 = 3.0;
        
        MooneyRivlinMaterialLaw<3> ml_law_3d(c1, c2);
        
        TS_ASSERT_DELTA(ml_law_3d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(ml_law_3d.GetC2(), c2, 1e-12);
        TS_ASSERT_DELTA(ml_law_3d.Get_dW_dI1(1.0,0.0), c1, 1e-12);
        TS_ASSERT_DELTA(ml_law_3d.Get_dW_dI2(1.0,0.0), c2, 1e-12);
        TS_ASSERT_DELTA(ml_law_3d.Get_d2W_dI1(1.0,0.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(ml_law_3d.Get_d2W_dI2(1.0,0.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(ml_law_3d.Get_d2W_dI1I2(1.0,0.0), 0.0, 1e-12);
        
        TS_ASSERT_DELTA(ml_law_3d.GetZeroStrainPressure(), 2*c1+4*c2, 1e-12);

        // compute the stress given C=delta_{MN} and p=zero_strain_pressure,
        // obviously it should be zero
        Tensor<2,3> identity_strain_3d;
        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                identity_strain_3d[i][j]=0;
            }
            identity_strain_3d[i][i]=1;
        }
        SymmetricTensor<2,3> T_3d;
        ml_law_3d.Compute2ndPiolaKirchoffStress(identity_strain_3d,
                                                ml_law_3d.GetZeroStrainPressure(),
                                                T_3d);
        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(T_3d[i][j],0.0,1e-12);
            }
        }

        // compute stress given a non-zero deformation
        Tensor<2,3> F;
        F[0][0] = 3.0;
        F[0][1] = 1.0;
        F[1][0] = -1.0;
        F[0][2] = 2.0;
        F[2][0] = 1.0;
        F[1][1] = 6.0;
        F[1][2] = -1.0;
        F[2][1] = 1.5;
        F[2][2] = 0.5;
        
        Tensor<2,3> C = transpose(F)*F;
        
        double I1 = C[0][0]+C[1][1]+C[2][2];
        
        Tensor<2,3> invC = invert(C);
        
        double pressure = 5.0;
        
        SymmetricTensor<2,3> T;
        SymmetricTensor<2,3> T2;
        Tensor<2,3> S;
        Tensor<2,3> sigma;

        double dTdE[3][3][3][3];
        
        ml_law_3d.ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, true);
        ml_law_3d.Compute1stPiolaKirchoffStress(F,pressure,S);
        ml_law_3d.Compute2ndPiolaKirchoffStress(C,pressure,T2);
        ml_law_3d.ComputeCauchyStress(F,pressure,sigma);
        
        // can't seem to do F*T where T is a SymmetricTensor<2,3>. Copy to
        // a Tensor<2,3>
        Tensor<2,3> T_as_unsym_tensor(T);
        Tensor<2,3> F_T_tranF_over_detF = (1.0/determinant(F))*F*T_as_unsym_tensor*transpose(F);
        Tensor<2,3> T_transposeF = T_as_unsym_tensor*transpose(F); 


        // check sigma is correct - sigma should be (1/detF) F * T * trans(F)
        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(sigma[i][j], F_T_tranF_over_detF[i][j], 1e-12);
            }
        }

        // check S is correct
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned i=0; i<3; i++)
            {
                TS_ASSERT_DELTA(S[M][i], T_transposeF[M][i], 1e-12);
            }
        }

        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                // check we gave a symmetric C
                assert(C[M][N]==C[N][M]);
                
                // check the stress
                TS_ASSERT_DELTA(T[M][N], (2*c1+2*c2*I1)*(M==N) - 2*c2*C[M][N] - pressure*invC[M][N], 1e-12);
                
                // check alternative computation of the stress
                TS_ASSERT_DELTA(T[M][N], T2[M][N], 1e-12);

                for(unsigned P=0;P<3;P++)
                {
                    for(unsigned Q=0;Q<3;Q++)
                    {
                        double true_val =   4*c2*((M==N)*(P==Q)-(M==P)*(N==Q))
                                          + 2*pressure*invC[M][P]*invC[Q][N];
                                          
                        TS_ASSERT_DELTA(dTdE[M][N][P][Q], true_val, 1e-12);
                    }
                }
            }
        }
    }
    
    
    void testExponentialLaw()
    {
        TS_ASSERT_THROWS_ANYTHING(ExponentialMaterialLaw<2> bad_exp_law(-3.0,1));
        TS_ASSERT_THROWS_ANYTHING(ExponentialMaterialLaw<1> bad_exp_law3(3.0,1));
        
        double a = 2.0;
        double b = 3.0;
        double I1 = 2.0;
        double I2 = 1.0;
        
        ExponentialMaterialLaw<2> exp_law_2d(a,b);
        
        TS_ASSERT_DELTA(exp_law_2d.GetA(), a, 1e-12);
        TS_ASSERT_DELTA(exp_law_2d.GetB(), b, 1e-12);
        
        TS_ASSERT_DELTA(exp_law_2d.Get_dW_dI1(I1,I2),  a*b*exp(b*(I1-2)),   1e-12);
        TS_ASSERT_DELTA(exp_law_2d.Get_d2W_dI1(I1,I2), b*exp_law_2d.Get_dW_dI1(I1,I2), 1e-12);
        
        ExponentialMaterialLaw<3> exp_law_3d(a,b);
        
        TS_ASSERT_DELTA(exp_law_3d.GetA(), a, 1e-12);
        TS_ASSERT_DELTA(exp_law_3d.GetB(), b, 1e-12);
        
        TS_ASSERT_DELTA(exp_law_3d.Get_dW_dI1(I1,I2),   a*b*exp(b*(I1-3)),   1e-12);
        TS_ASSERT_DELTA(exp_law_3d.Get_dW_dI2(I1,I2),   0.0,                 1e-12);

        TS_ASSERT_DELTA(exp_law_3d.Get_d2W_dI1(I1,I2),  b*exp_law_3d.Get_dW_dI1(I1,I2), 1e-12);
        TS_ASSERT_DELTA(exp_law_3d.Get_d2W_dI2(I1,I2),  0.0,                 1e-12);
        TS_ASSERT_DELTA(exp_law_3d.Get_d2W_dI1I2(I1,I2),0.0,                 1e-12);


        TS_ASSERT_DELTA(exp_law_3d.GetZeroStrainPressure(), 2*a*b, 1e-12);

        // compute the stress given C=delta_{MN} and p=zero_strain_pressure,
        // obviously it should be zero
        Tensor<2,3> identity_strain_3d;
        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                identity_strain_3d[i][j]=0;
            }
            identity_strain_3d[i][i]=1;
        }
        SymmetricTensor<2,3> T_3d;
        exp_law_3d.Compute2ndPiolaKirchoffStress(identity_strain_3d,
                                                exp_law_3d.GetZeroStrainPressure(),
                                                T_3d);
        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(T_3d[i][j],0.0,1e-12);
            }
        }
    }
    
    
    // The polynomial material law with N=1 is exactly a mooney-rivlin law and should
    // agree
    void testPolynomailMaterialLawAgainstMooneyRivlin()
    {
        unsigned N = 1;
        std::vector< std::vector<double> > alpha = PolynomialMaterialLaw3d::GetZeroedAlpha(N);
        
        // test GetZeroedAlpha
        TS_ASSERT_EQUALS(alpha.size(),2);
        TS_ASSERT_EQUALS(alpha[0].size(),2);
        TS_ASSERT_EQUALS(alpha[1].size(),2);

        TS_ASSERT_DELTA(alpha[0][0],0.0,1e-12);
        TS_ASSERT_DELTA(alpha[0][1],0.0,1e-12);
        TS_ASSERT_DELTA(alpha[1][0],0.0,1e-12);
        TS_ASSERT_DELTA(alpha[1][1],0.0,1e-12);
        

        double c1 = 3.0;
        double c2 = 2.0;

        alpha[1][0] = c1;
        alpha[0][1] = c2;
        
        PolynomialMaterialLaw3d poly_mr_law(N,alpha);
        MooneyRivlinMaterialLaw<3> mooney_rivlin_law(c1,c2);
        
        double I1 = 4;
        double I2 = 2.4;

        TS_ASSERT_DELTA(mooney_rivlin_law.Get_dW_dI1(I1,I2),    poly_mr_law.Get_dW_dI1(I1,I2),    1e-12);
        TS_ASSERT_DELTA(mooney_rivlin_law.Get_dW_dI2(I1,I2),    poly_mr_law.Get_dW_dI2(I1,I2),    1e-12);
        TS_ASSERT_DELTA(mooney_rivlin_law.Get_d2W_dI1(I1,I2),   poly_mr_law.Get_d2W_dI1(I1,I2),   1e-12);
        TS_ASSERT_DELTA(mooney_rivlin_law.Get_d2W_dI2(I1,I2),   poly_mr_law.Get_d2W_dI2(I1,I2),   1e-12);
        TS_ASSERT_DELTA(mooney_rivlin_law.Get_d2W_dI1I2(I1,I2), poly_mr_law.Get_d2W_dI1I2(I1,I2), 1e-12);
    }                
    
    
    
    // Test the Polynomial Material Law with a quadratic law
    //   W = c20 (I1-3)^2  +  c11 (I1-3)(I2-3)  +  c02 (I2-3)^2   -   p C^{-1}/2
    // 
    // We test ComputeStressAndStressDerivative() with this law because it uses all the
    // bits in that method. 
    void testQuadraticPolynomialLaw()
    {
        unsigned N = 2;
        std::vector< std::vector<double> > alpha = PolynomialMaterialLaw3d::GetZeroedAlpha(N);
        
        double c20 = 3.0;
        double c11 = 4.0;
        double c02 = 5.0;
        
        alpha[2][0] = c20;
        alpha[1][1] = c11;
        alpha[0][2] = c02;

        PolynomialMaterialLaw3d poly_law(N,alpha);

        Tensor<2,3> C;
        C[0][0] = 3.0;
        C[0][1] = 1.0;
        C[1][0] = 1.0;
        C[0][2] = 2.0;
        C[2][0] = 2.0;
        C[1][1] = 6.0;
        C[1][2] = -1.0;
        C[2][1] = -1.0;
        C[2][2] = 0.5;
        
        double I1 =   C[0][0]+C[1][1]+C[2][2];
        double I2 =   C[0][0]*C[1][1] + C[1][1]*C[2][2] + C[2][2]*C[0][0]
                    - C[0][1]*C[1][0] - C[1][2]*C[2][1] - C[2][0]*C[0][2];


        double true_dWdI1    = 2*c20*(I1-3) +   c11*(I2-3);
        double true_dWdI2    =   c11*(I1-3) + 2*c02*(I2-3);
        double true_d2WdI1   = 2*c20;
        double true_d2WdI1I2 =   c11;
        double true_d2WdI2   = 2*c02;

        TS_ASSERT_DELTA(true_dWdI1,    poly_law.Get_dW_dI1(I1,I2),    1e-12);
        TS_ASSERT_DELTA(true_dWdI2,    poly_law.Get_dW_dI2(I1,I2),    1e-12);
        TS_ASSERT_DELTA(true_d2WdI1,   poly_law.Get_d2W_dI1(I1,I2),   1e-12);
        TS_ASSERT_DELTA(true_d2WdI2,   poly_law.Get_d2W_dI2(I1,I2),   1e-12);
        TS_ASSERT_DELTA(true_d2WdI1I2, poly_law.Get_d2W_dI1I2(I1,I2), 1e-12);


        Tensor<2,3> invC = invert(C);
        
        double pressure = 5.0;
        
        SymmetricTensor<2,3> T;
        double dTdE[3][3][3][3];

        poly_law.ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, true);
        
        for(unsigned M=0; M<3; M++)
        {
            for(unsigned N=0; N<3; N++)
            {
                // check we gave a symmetric C
                assert(C[M][N]==C[N][M]);

                double dI1_dC_MN = (M==N); 
                double dI2_dC_MN = I1*(M==N)-C[M][N]; 

                
                double true_val =   2 * true_dWdI1 * dI1_dC_MN 
                                  + 2 * true_dWdI2 * dI2_dC_MN
                                  - pressure*invC[M][N];
                
                TS_ASSERT_DELTA(  T[M][N], true_val, 1e-12 );
                
                for(unsigned P=0;P<3;P++)
                {
                    for(unsigned Q=0;Q<3;Q++)
                    {
                        double dI1_dC_MN = (M==N); 
                        double dI1_dC_PQ = (P==Q); 

                        double d2I1_dC2  = 0;
                        
                        double dI2_dC_MN = I1*(M==N)-C[M][N]; 
                        double dI2_dC_PQ = I1*(P==Q)-C[P][Q];                          
                        
                        double d2I2_dC2  = (M==N)*(P==Q)-(M==P)*(N==Q);
                        
                        double d_invC_dC = -invC[M][P]*invC[Q][N];
                        
                        double true_val =    4 * true_d2WdI1 * dI1_dC_MN * dI1_dC_PQ
                                           + 4 * true_dWdI1  * d2I1_dC2
                                           + 4 * true_d2WdI2 * dI2_dC_MN * dI2_dC_PQ
                                           + 4 * true_dWdI2  * d2I2_dC2
                                           + 4 * true_d2WdI1I2 * (dI1_dC_MN*dI2_dC_PQ + dI1_dC_PQ*dI2_dC_MN)
                                           - 2 * pressure * d_invC_dC;
                                          
                        TS_ASSERT_DELTA(dTdE[M][N][P][Q], true_val, 1e-12);
                    }
                }
            }
        }
        
        // cover GetAlpha
        TS_ASSERT_DELTA(poly_law.GetAlpha(0,1),0.0,1e-12);
        TS_ASSERT_DELTA(poly_law.GetAlpha(0,2),c02,1e-12);
        TS_ASSERT_DELTA(poly_law.GetAlpha(1,1),c11,1e-12);
        TS_ASSERT_DELTA(poly_law.GetAlpha(1,0),0.0,1e-12);
        TS_ASSERT_DELTA(poly_law.GetAlpha(2,0),c20,1e-12);
        
        
        // check exception thrown if N=0
        TS_ASSERT_THROWS_ANYTHING(PolynomialMaterialLaw3d bad_poly1(0,alpha));

        // check exception thrown if alpha is not correctly sized
        std::vector< std::vector<double> > bad_alpha(2);
        bad_alpha[0].resize(1);
        bad_alpha[1].resize(1);
        TS_ASSERT_THROWS_ANYTHING(PolynomialMaterialLaw3d bad_poly2(2,bad_alpha));

        // check exception thrown if alpha[p][q]!=0, when p+q>N
        alpha[2][2] = 1.0;
        TS_ASSERT_THROWS_ANYTHING(PolynomialMaterialLaw3d bad_poly3(2,alpha));
     
       // compute the stress given C=delta_{MN} and p=zero_strain_pressure,
        // obviously it should be zero
        Tensor<2,3> identity_strain_3d;
        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                identity_strain_3d[i][j]=0;
            }
            identity_strain_3d[i][i]=1;
        }
        SymmetricTensor<2,3> T_3d;
        

        poly_law.Compute2ndPiolaKirchoffStress(identity_strain_3d,
                                               poly_law.GetZeroStrainPressure(),
                                               T_3d);
        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(T_3d[i][j],0.0,1e-12);
            }
        }
    }
};

#endif /*TESTMATERIALLAWS_HPP_*/
