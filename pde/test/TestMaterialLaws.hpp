/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTMATERIALLAWS_HPP_
#define TESTMATERIALLAWS_HPP_

#include <cxxtest/TestSuite.h>

#include "MooneyRivlinMaterialLaw.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "PoleZeroMaterialLaw.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "SchmidCostaExponentialLaw2d.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "ToyCompressibleMaterialLaw.hpp"
#include <cassert>


class TestMaterialLaws : public CxxTest::TestSuite
{
public:
    void TestMooneyRivlinLaw()
    {
        TS_ASSERT_THROWS_THIS(MooneyRivlinMaterialLaw<2> bad_mr_law(-3.0),"c1 must be positive in mooney-rivlin");
        TS_ASSERT_THROWS_THIS(MooneyRivlinMaterialLaw<3> bad_mr_law2(3.0),"Two parameters needed for 3d Mooney-Rivlin");

        double c1 = 2.0;

        MooneyRivlinMaterialLaw<2> law_2d(c1);

        TS_ASSERT_DELTA(law_2d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI1(1.0,0.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI1(1.0,0.0), 0.0, 1e-12);

        // compute the stress given C=delta_{MN} and p=zero_strain_pressure,
        // obviously it should be zero
        TS_ASSERT_DELTA(law_2d.GetZeroStrainPressure(), 2*c1, 1e-12);
        c_matrix<double,2,2> identity_strain_2d = identity_matrix<double>(2);

        c_matrix<double,2,2> T_2d;
        law_2d.Compute2ndPiolaKirchoffStress(identity_strain_2d,
                                                law_2d.GetZeroStrainPressure(),
                                                T_2d);
        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                TS_ASSERT_DELTA(T_2d(i,j),0.0,1e-12);
            }
        }


        double c2 = 3.0;

        MooneyRivlinMaterialLaw<3> law_3d(c1, c2);

        TS_ASSERT_DELTA(law_3d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.GetC2(), c2, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI1(1.0,0.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI2(1.0,0.0), c2, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_d2W_dI1(1.0,0.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_d2W_dI2(1.0,0.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_d2W_dI1I2(1.0,0.0), 0.0, 1e-12);

        TS_ASSERT_DELTA(law_3d.GetZeroStrainPressure(), 2*c1+4*c2, 1e-12);

        // compute the stress given C=delta_{MN} and p=zero_strain_pressure,
        // obviously it should be zero
        c_matrix<double,3,3> identity_strain_3d = identity_matrix<double>(3);
        c_matrix<double,3,3> T_3d;
        law_3d.Compute2ndPiolaKirchoffStress(identity_strain_3d,
                                                law_3d.GetZeroStrainPressure(),
                                                T_3d);
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(T_3d(i,j),0.0,1e-12);
            }
        }

        // compute stress given a non-zero deformation
        c_matrix<double,3,3> F;
        F(0,0) = 3.0;
        F(0,1) = 1.0;
        F(1,0) = -1.0;
        F(0,2) = 2.0;
        F(2,0) = 1.0;
        F(1,1) = 6.0;
        F(1,2) = -1.0;
        F(2,1) = 1.5;
        F(2,2) = 0.5;

        c_matrix<double,3,3> C = prod(trans(F),F);

        double I1 =  Trace(C);

        c_matrix<double,3,3> invC = Inverse(C);

        double pressure = 5.0;

        c_matrix<double,3,3> T;
        c_matrix<double,3,3> T2;
        c_matrix<double,3,3> S;
        c_matrix<double,3,3> sigma;

        FourthOrderTensor<3,3,3,3> dTdE;

        law_3d.ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, true);
        law_3d.Compute1stPiolaKirchoffStress(F,pressure,S);
        law_3d.Compute2ndPiolaKirchoffStress(C,pressure,T2);
        law_3d.ComputeCauchyStress(F,pressure,sigma);

        c_matrix<double,3,3> FT = prod(F,T);
        c_matrix<double,3,3> F_T_tranF_over_detF = (1.0/Determinant(F))*prod(FT,trans(F));//F*T_as_unsym_tensor*transpose(F);

        c_matrix<double,3,3> T_transposeF = prod(T,trans(F));//T_as_unsym_tensor*transpose(F);


        // check sigma is correct - sigma should be (1/detF) F * T * trans(F)
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(sigma(i,j), F_T_tranF_over_detF(i,j), 1e-12);
            }
        }

        // check S is correct
        for (unsigned M=0; M<3; M++)
        {
            for (unsigned i=0; i<3; i++)
            {
                TS_ASSERT_DELTA(S(M,i), T_transposeF(M,i), 1e-12);
            }
        }

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                // check we gave a symmetric C
                assert(C(M,N)==C(N,M));

                // check the stress
                TS_ASSERT_DELTA(T(M,N), (2*c1+2*c2*I1)*(M==N) - 2*c2*C(M,N) - pressure*invC(M,N), 1e-12);

                // check alternative computation of the stress
                TS_ASSERT_DELTA(T(M,N), T2(M,N), 1e-12);

                for (unsigned P=0;P<3;P++)
                {
                    for (unsigned Q=0;Q<3;Q++)
                    {
                        double true_val =   4*c2*((M==N)*(P==Q)-(M==P)*(N==Q))
                                            + 2*pressure*invC(M,P)*invC(Q,N);

                        TS_ASSERT_DELTA(dTdE(M,N,P,Q), true_val, 1e-12);
                    }
                }
            }
        }

        law_3d.ScaleMaterialParameters(10);
        TS_ASSERT_DELTA(law_3d.GetC1(), c1/10, 1e-12);
        TS_ASSERT_DELTA(law_3d.GetC2(), c2/10, 1e-12);
    }


    void TestExponentialLaw()
    {
        TS_ASSERT_THROWS_THIS(ExponentialMaterialLaw<2> bad_exp_law(-3.0,1),"a must be positive");

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
        c_matrix<double,3,3> identity_strain_3d = identity_matrix<double>(3);
        c_matrix<double,3,3> T_3d;
        exp_law_3d.Compute2ndPiolaKirchoffStress(identity_strain_3d,
                                                 exp_law_3d.GetZeroStrainPressure(),
                                                 T_3d);
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(T_3d(i,j),0.0,1e-12);
            }
        }
    }


    // The polynomial material law with N=1 is exactly a mooney-rivlin law and should
    // agree
    void TestPolynomailMaterialLawAgainstMooneyRivlin()
    {
        unsigned N = 1;
        std::vector< std::vector<double> > alpha = PolynomialMaterialLaw3d::GetZeroedAlpha(N);

        // test GetZeroedAlpha
        TS_ASSERT_EQUALS(alpha.size(),2u);
        TS_ASSERT_EQUALS(alpha[0].size(),2u);
        TS_ASSERT_EQUALS(alpha[1].size(),2u);

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
    void TestQuadraticPolynomialLaw()
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

        c_matrix<double,3,3> C;
        C(0,0) = 3.0;
        C(0,1) = 1.0;
        C(1,0) = 1.0;
        C(0,2) = 2.0;
        C(2,0) = 2.0;
        C(1,1) = 6.0;
        C(1,2) = -1.0;
        C(2,1) = -1.0;
        C(2,2) = 0.5;

        double I1 =  Trace(C);
        double I2 =  SecondInvariant(C);


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


        c_matrix<double,3,3> invC = Inverse(C);

        double pressure = 5.0;

        c_matrix<double,3,3> T;
        FourthOrderTensor<3,3,3,3> dTdE;

        poly_law.ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, true);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                // check we gave a symmetric C
                assert(C(M,N)==C(N,M));

                double dI1_dC_MN = (M==N);
                double dI2_dC_MN = I1*(M==N)-C(M,N);


                double true_val =   2 * true_dWdI1 * dI1_dC_MN
                                    + 2 * true_dWdI2 * dI2_dC_MN
                                    - pressure*invC(M,N);

                TS_ASSERT_DELTA(  T(M,N), true_val, 1e-12 );

                for (unsigned P=0;P<3;P++)
                {
                    for (unsigned Q=0;Q<3;Q++)
                    {
                        double dI1_dC_MN = (M==N);
                        double dI1_dC_PQ = (P==Q);

                        double d2I1_dC2  = 0;

                        double dI2_dC_MN = I1*(M==N)-C(M,N);
                        double dI2_dC_PQ = I1*(P==Q)-C(P,Q);

                        double d2I2_dC2  = (M==N)*(P==Q)-(M==P)*(N==Q);

                        double d_invC_dC = -invC(M,P)*invC(Q,N);

                        double true_val =    4 * true_d2WdI1 * dI1_dC_MN * dI1_dC_PQ
                                             + 4 * true_dWdI1  * d2I1_dC2
                                             + 4 * true_d2WdI2 * dI2_dC_MN * dI2_dC_PQ
                                             + 4 * true_dWdI2  * d2I2_dC2
                                             + 4 * true_d2WdI1I2 * (dI1_dC_MN*dI2_dC_PQ + dI1_dC_PQ*dI2_dC_MN)
                                             - 2 * pressure * d_invC_dC;

                        TS_ASSERT_DELTA(dTdE(M,N,P,Q), true_val, 1e-12);
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
        TS_ASSERT_THROWS_THIS(PolynomialMaterialLaw3d bad_poly1(0,alpha),"n must be positive");

        // check exception thrown if alpha is not correctly sized
        std::vector< std::vector<double> > bad_alpha(2);
        bad_alpha[0].resize(1);
        bad_alpha[1].resize(1);
        TS_ASSERT_THROWS_THIS(PolynomialMaterialLaw3d bad_poly2(2,bad_alpha),"alpha not big enough");

        // check exception thrown if alpha[p][q]!=0, when p+q>N
        alpha[2][2] = 1.0;
        TS_ASSERT_THROWS_THIS(PolynomialMaterialLaw3d bad_poly3(2,alpha),"alpha[2][2] should be zero, as p+q > 2");

        // compute the stress given C=delta_{MN} and p=zero_strain_pressure,
        // obviously it should be zero
        c_matrix<double,3,3> identity_strain_3d = identity_matrix<double>(3);
        c_matrix<double,3,3> T_3d;

        poly_law.Compute2ndPiolaKirchoffStress(identity_strain_3d,
                                               poly_law.GetZeroStrainPressure(),
                                               T_3d);
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(T_3d(i,j),0.0,1e-12);
            }
        }
    }

    void TestPoleZeroMaterialLaw()
    {
        std::vector<std::vector<double> > k(2),a(2),b(2);
        for (unsigned i=0; i<2; i++)
        {
            k[i].resize(2);
            a[i].resize(2);
            b[i].resize(2);
        }

        k[0][0] = 1;
        k[1][0] = k[0][1] = 2;
        k[1][1] = 3;

        a[0][0] = 4;
        a[1][0] = a[0][1] = 5;
        a[1][1] = 6;

        b[0][0] = 7;
        b[1][0] = b[0][1] = 6;
        b[1][1] = 5;

        PoleZeroMaterialLaw<2> pole_zero_law(k,a,b);

        c_matrix<double,2,2> C;
        C(0,0) = 1;
        C(0,1) = C(1,0) = 0;
        C(1,1) = 1;
        c_matrix<double,2,2> invC = Inverse(C);

        c_matrix<double,2,2> T;
        FourthOrderTensor<2,2,2,2> dTdE;
        double pressure = pole_zero_law.GetZeroStrainPressure();

        pole_zero_law.ComputeStressAndStressDerivative(C,invC,pressure,T,dTdE,true);

        TS_ASSERT_DELTA( T(0,0), 0.0, 1e-9 );
        TS_ASSERT_DELTA( T(1,0), 0.0, 1e-9 );
        TS_ASSERT_DELTA( T(0,1), 0.0, 1e-9 );
        TS_ASSERT_DELTA( T(1,1), 0.0, 1e-9 );

        // non-trivial deformation, (checking all components have such that E_MN < a_MN)
        C(0,0) = 0.5;
        C(0,1) = -0.1;
        C(1,0) = -0.1;
        C(1,1) = 0.5;
        invC = Inverse(C);

//// currently been changed on that pole-zero law DOESN'T return T=0 if E<0
//        // C such that E_MN < 0, p=0 => T=0, dTdE=0;
//        pole_zero_law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);
//
//        for (unsigned M=0; M<2; M++)
//        {
//            for (unsigned N=0; N<2; N++)
//            {
//                TS_ASSERT_DELTA(T(M,N), 0.0, 1e-9);
//                for (unsigned P=0; P<2; P++)
//                {
//                    for (unsigned Q=0; Q<2; Q++)
//                    {
//                        TS_ASSERT_DELTA(dTdE(M,N,P,Q), 0.0, 1e-9);
//                    }
//                }
//            }
//        }

        // non-trivial deformation, (checking all components have such that E_MN < a_MN)
        C(0,0) = 2;
        C(0,1) = 2;
        C(1,0) = 2;
        C(1,1) = 5;

        pole_zero_law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        // T_MN = ke(2+be/(a-e))/(a-e)^b
        double t00 = 1*0.5*3/6433.92969;
        double t10 = 2.0*1.0*(2+6.0/4.0)/4096.0;
        double t11 = 3.0*2.0*(2+5*2.0/4.0)/1024.0;

        TS_ASSERT_DELTA( T(0,0), t00, 1e-9 );
        TS_ASSERT_DELTA( T(1,0), t10, 1e-9 );
        TS_ASSERT_DELTA( T(0,1), t10, 1e-9 );
        TS_ASSERT_DELTA( T(1,1), t11, 1e-9 );

        // test dTdE
        double dtde00 = ( 2*3.5*3.5 + 14*3.5 + 0.25*7*8 )/78815.6387;
        double dtde10 = 2*( 2*4*4 + 4*1*6*4 + 42)/65536.0;
        double dtde11 = 3*( 2*4*4 + 4*2*5*4 + 120)/16384.0;

        TS_ASSERT_DELTA(dTdE(0,0,0,0), dtde00, 1e-9);
        TS_ASSERT_DELTA(dTdE(0,1,0,1), dtde10, 1e-9);
        TS_ASSERT_DELTA(dTdE(1,0,1,0), dtde10, 1e-9);
        TS_ASSERT_DELTA(dTdE(1,1,1,1), dtde11, 1e-9);

        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        if ((P!=M) || (Q!=N))
                        {
                            TS_ASSERT_DELTA(dTdE(M,N,P,Q), 0.0, 1e-9);
                        }
                    }
                }
            }
        }

        // test dTdE via a numerical derivative
        double old_T00 = T(0,0);
        double h=0.001;
        C(0,0)+=h;
        invC = Inverse(C);
        pole_zero_law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,false);

        double dTdC00_numerical = (T(0,0) - old_T00)/h;  // dC not dE
        TS_ASSERT_DELTA( dTdE(0,0,0,0), 2*dTdC00_numerical, 1e-4);


        // test the pressure terms in the stress and stress-deriv, by calling with
        // p=0 and p=1 and verifying the difference is what it should be
        c_matrix<double,2,2> T2;
        FourthOrderTensor<2,2,2,2> dTdE2;
        pole_zero_law.ComputeStressAndStressDerivative(C, invC, 0.0, T,  dTdE,  true);
        pole_zero_law.ComputeStressAndStressDerivative(C, invC, 1.0, T2, dTdE2, true);

        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                TS_ASSERT_DELTA(T(M,N) - T2(M,N), invC(M,N), 1e-6);
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(dTdE(M,N,P,Q)-dTdE2(M,N,P,Q), -2*invC(M,P)*invC(Q,N), 1e-6);
                    }
                }
            }
        }
    }

    void TestPoleZeroMaterialLaw3d()
    {
        std::vector<std::vector<double> > k(3),a(3),b(3);
        for (unsigned i=0; i<3; i++)
        {
            k[i].resize(3);
            a[i].resize(3);
            b[i].resize(3);
        }

        k[0][0] = 1;
        k[1][0] = k[0][1] = 2;
        k[0][2] = k[2][0] = 4;
        k[1][1] = 3;
        k[1][2] = k[2][1] = 5;
        k[2][2] = 3;

        a[0][0] = 4;
        a[1][0] = a[0][1] = 5;
        a[2][0] = a[0][2] = 6;
        a[1][1] = 6;
        a[2][1] = a[1][2] = 4;
        a[2][2] = 9;

        b[0][0] = 7;
        b[1][0] = b[0][1] = 6;
        b[2][0] = b[0][2] = 2;
        b[1][1] = 5;
        b[2][1] = b[1][2] = 4;
        b[2][2] = 2;

        PoleZeroMaterialLaw<3> pole_zero_law(k,a,b);

        c_matrix<double,3,3> C;
        C(0,0) = 2;
        C(0,1) = C(1,0) = 2;
        C(0,2) = C(2,0) = 3;
        C(1,1) = 5;
        C(1,2) = C(2,1) = 4;
        C(2,2) = 3;

        c_matrix<double,3,3> invC = Inverse(C);

        c_matrix<double,3,3> T;
        FourthOrderTensor<3,3,3,3> dTdE;

        pole_zero_law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        //// same as previous test, except the t22 and dtde22 bits are new
        double t00 = 1*0.5*3/6433.92969;
        double t10 = 2.0*1.0*(2+6.0/4.0)/4096.0;
        double t11 = 3.0*2.0*(2+5*2.0/4.0)/1024.0;
        double t22 = 3.0*(2+2.0/8.0)/64.0;

        TS_ASSERT_DELTA( T(0,0), t00, 1e-9 );
        TS_ASSERT_DELTA( T(1,0), t10, 1e-9 );
        TS_ASSERT_DELTA( T(0,1), t10, 1e-9 );
        TS_ASSERT_DELTA( T(1,1), t11, 1e-9 );
        TS_ASSERT_DELTA( T(2,2), t22, 1e-9 );

        // test dTdE
        double dtde00 = ( 2*3.5*3.5 + 14*3.5 + 0.25*7*8 )/78815.6387;
        double dtde10 = 2*( 2*4*4 + 4*1*6*4 + 42)/65536.0;
        double dtde11 = 3*( 2*4*4 + 4*2*5*4 + 120)/16384.0;
        double dtde22 = 3*( 2*64 + 4*1*2*8 + 2*3)/4096.0;

        TS_ASSERT_DELTA(dTdE(0,0,0,0), dtde00, 1e-9);
        TS_ASSERT_DELTA(dTdE(0,1,0,1), dtde10, 1e-9);
        TS_ASSERT_DELTA(dTdE(1,0,1,0), dtde10, 1e-9);
        TS_ASSERT_DELTA(dTdE(1,1,1,1), dtde11, 1e-9);
        TS_ASSERT_DELTA(dTdE(2,2,2,2), dtde22, 1e-9);


        pole_zero_law.ScaleMaterialParameters(10);
        TS_ASSERT_DELTA(pole_zero_law.mK[0][0], 0.1, 1e-12);
        TS_ASSERT_DELTA(pole_zero_law.mA[0][0], 4, 1e-12);
        TS_ASSERT_DELTA(pole_zero_law.mB[0][0], 7, 1e-12);
    }

//    void dontTestPoleZero3dIn1dLaw() throw(Exception)
//    {
//        PoleZero3dIn1dLaw law;
//
//        TS_ASSERT_DELTA( law.GetT(0), 0.0, 1e-12 );
//        TS_ASSERT_DELTA( law.GetT(0.1), 2.0809, 1e-3 );
//        TS_ASSERT_DELTA( law.GetT(0.2), 8.5158, 1e-3 );
//        TS_ASSERT_DELTA( law.GetT(0.3), 37.0291, 1e-3 );
//
//        TS_ASSERT_DELTA( law.GetT(-0.1), -0.5023, 1e-3 );
//        TS_ASSERT_DELTA( law.GetT(-0.2), -2.2589, 1e-3 );
//
//        law.SetUpStores();
//
//        TS_ASSERT_DELTA( law.GetT(-0.1), -0.5023, 1e-3 );
//        TS_ASSERT_DELTA( law.GetT(-0.2), -2.2589, 1e-3 );
//    }

    void TestNashHunterPoleZeroLaw3d() throw(Exception)
    {
        NashHunterPoleZeroLaw<3> law;

        c_matrix<double,3,3> C;
        c_matrix<double,3,3> invC;
        C(0,0) = 1.2;
        C(0,1) = C(1,0) = 0.1;
        C(0,2) = C(2,0) = 0.3;
        C(1,1) = 1.1;
        C(1,2) = C(2,1) = -0.1;
        C(2,2) = 1.3;
        invC = Inverse(C);

        c_matrix<double,3,3> T;
        FourthOrderTensor<3,3,3,3> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        // hard-coded test to test nothing changes (and stresses are of
        // the correct magnitude, which is dependent on whether the params
        // have been entered a Pa or KPa)
        TS_ASSERT_DELTA(T(0,0),2.0902,1e-3);
    }

    void TestNashHunterPoleZeroLaw2d() throw(Exception)
    {
        NashHunterPoleZeroLaw<2> law;

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.2;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 1.1;
        invC = Inverse(C);

        c_matrix<double,2,2> T;
        FourthOrderTensor<2,2,2,2> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        // hard-coded test to test nothing changes (and stresses are of
        // the correct magnitude, which is dependent on whether the params
        // have been entered a Pa or KPa)
        TS_ASSERT_DELTA(T(0,0),2.0902,1e-3);

        C(0,0) = 10;
        TS_ASSERT_THROWS_CONTAINS(law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true), "strain unacceptably large");
    }

    void TestNashHunterPoleZeroChangeOfBasis() throw(Exception)
    {
        NashHunterPoleZeroLaw<2> law;

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.2;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 1.1;
        invC = Inverse(C);

        c_matrix<double,2,2> T_Xfibres;
        c_matrix<double,2,2> T_Yfibres;
        FourthOrderTensor<2,2,2,2> dTdE_Xfibres;
        FourthOrderTensor<2,2,2,2> dTdE_Yfibres;

        double p = 1.0;
        law.ComputeStressAndStressDerivative(C,invC,p,T_Xfibres,dTdE_Xfibres,true); // no change of basis no fibres in X-dir

        // now assume fibres in Y-dir. first set up equivalent C
        C(0,0) = 1.1;
        C(1,1) = 1.2;
        invC = Inverse(C);

        // change of basis matrix
        c_matrix<double,2,2> P;
        P(0,0) = P(1,1) = 0.0;
        P(1,0) = P(0,1) = 1.0;

        law.SetChangeOfBasisMatrix(P);
        law.ComputeStressAndStressDerivative(C,invC,p,T_Yfibres,dTdE_Yfibres,true);

        TS_ASSERT_DELTA(T_Xfibres(0,0), T_Yfibres(1,1), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(1,1), T_Yfibres(0,0), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(0,1), T_Yfibres(1,0), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(1,0), T_Yfibres(0,1), 1e-8);

        // dTdE_Xfibres(0,1,1,0) should be able to dTdE_Yfibres(1,0,0,1), etc
        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(dTdE_Xfibres(M,N,P,Q), dTdE_Yfibres((M+1)%2,(N+1)%2,(P+1)%2,(Q+1)%2), 1e-8);
                    }
                }
            }
        }
    }

    void TestDerivateInPoleZeroLaw2d() throw(Exception)
    {
        NashHunterPoleZeroLaw<2> law;

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.1;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 0.9;
        invC = Inverse(C);

        c_matrix<double,2,2> T_base;
        FourthOrderTensor<2,2,2,2> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T_base,dTdE,false);

        double h=0.0001;

        for (unsigned M=0; M<2; M++)
        {
            c_matrix<double,2,2> C;
            c_matrix<double,2,2> invC;

            C(0,0) = 1.1;
            C(0,1) = C(1,0) = 0.1;
            C(1,1) = 0.9;

            C(M,M) += h;     // just change C00 and C11. Can't see how to compute numerical
                              // derivative of wrt C01,C10, given the C has to be passed in symmetric
            invC = Inverse(C);

            c_matrix<double,2,2> T;

            law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

            for (unsigned P=0; P<2; P++)
            {
                for (unsigned Q=0; Q<2; Q++)
                {
                    double dtdc = (T(P,Q)-T_base(P,Q))/h;
                    //std::cout << P << Q << M << M << " " << dTdE(P,Q,M,M) << "\n";
                    TS_ASSERT_DELTA(2*dtdc, dTdE(P,Q,M,M), dTdE(P,Q,M,M)*1e-3);
                }
            }
            C(M,M) -= h;
        }
    }


    void TestSchmidCostaExponentialLaw()
    {
        SchmidCostaExponentialLaw2d law;

        double a = law.GetA();
        assert(a>0);
        double bff = law.GetB()[0][0];
        double bfs = law.GetB()[0][1];
        double bsf = law.GetB()[1][0];
        double bss = law.GetB()[1][1];

        TS_ASSERT(bsf == bfs);

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.1;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 0.9;
        invC = Inverse(C);

        c_matrix<double,2,2> T_base;
        FourthOrderTensor<2,2,2,2> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T_base,dTdE,false);

        double e00 = 0.5*(C(0,0)-1);
        double e01 = 0.5*C(0,1);
        double e11 = 0.5*(C(1,1)-1);
        double Q = bff*e00*e00 + 2*bfs*e01*e01 + bss*e11*e11;
        TS_ASSERT_DELTA(T_base(0,0), a*exp(Q-1)*bff*e00, 1e-9);
        TS_ASSERT_DELTA(T_base(0,1), a*exp(Q-1)*bfs*e01, 1e-9);
        TS_ASSERT_DELTA(T_base(1,0), a*exp(Q-1)*bsf*e01, 1e-9);
        TS_ASSERT_DELTA(T_base(1,1), a*exp(Q-1)*bss*e11, 1e-9);

        double h=0.00001;

        for (unsigned M=0; M<2; M++)
        {
            c_matrix<double,2,2> C;
            c_matrix<double,2,2> invC;

            C(0,0) = 1.1;
            C(0,1) = C(1,0) = 0.1;
            C(1,1) = 0.9;

            C(M,M) += h;     // just change C00 and C11. Can't see how to compute numerical
                              // derivative of wrt C01,C10, given the C is assumed symmetric
            invC = Inverse(C);

            c_matrix<double,2,2> T;

            law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

            for (unsigned P=0; P<2; P++)
            {
                for (unsigned Q=0; Q<2; Q++)
                {
                    double dtdc = (T(P,Q)-T_base(P,Q))/h;
                    //std::cout << P << Q << M << M << " " << dTdE(P,Q,M,M) << "\n";
                    TS_ASSERT_DELTA(2*dtdc, dTdE(P,Q,M,M), fabs(dTdE(P,Q,M,M)*1e-3));
                }
            }

            C(M,M) -= h;
        }

        // test get zero strain pressure
        C(0,0) = 1.0;
        C(0,1) = C(1,0) = 0.0;
        C(1,1) = 1.0;
        invC = Inverse(C);
        double p = law.GetZeroStrainPressure();
        law.ComputeStressAndStressDerivative(C,invC,p,T_base,dTdE,false);
        TS_ASSERT_DELTA(T_base(0,0), 0.0, 1e-9);
        TS_ASSERT_DELTA(T_base(0,1), 0.0, 1e-9);
        TS_ASSERT_DELTA(T_base(0,1), 0.0, 1e-9);
        TS_ASSERT_DELTA(T_base(1,1), 0.0, 1e-9);
    }

    void TestSchmidCostaChangeOfBasis() throw(Exception)
    {
        SchmidCostaExponentialLaw2d law;

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.2;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 1.1;
        invC = Inverse(C);

        c_matrix<double,2,2> T_Xfibres;
        c_matrix<double,2,2> T_Yfibres;
        FourthOrderTensor<2,2,2,2> dTdE_Xfibres;
        FourthOrderTensor<2,2,2,2> dTdE_Yfibres;

        double p = 1.0;
        law.ComputeStressAndStressDerivative(C,invC,p,T_Xfibres,dTdE_Xfibres,true); // no change of basis no fibres in X-dir

        // now assume fibres in Y-dir. first set up equivalent C
        C(0,0) = 1.1;
        C(1,1) = 1.2;
        invC = Inverse(C);

        // change of basis matrix
        c_matrix<double,2,2> P;
        P(0,0) = P(1,1) = 0.0;
        P(1,0) = P(0,1) = 1.0;

        law.SetChangeOfBasisMatrix(P);
        law.ComputeStressAndStressDerivative(C,invC,p,T_Yfibres,dTdE_Yfibres,true);

        TS_ASSERT_DELTA(T_Xfibres(0,0), T_Yfibres(1,1), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(1,1), T_Yfibres(0,0), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(0,1), T_Yfibres(1,0), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(1,0), T_Yfibres(0,1), 1e-8);

        // dTdE_Xfibres(0,1,1,0) should be able to dTdE_Yfibres(1,0,0,1), etc
        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(dTdE_Xfibres(M,N,P,Q), dTdE_Yfibres((M+1)%2,(N+1)%2,(P+1)%2,(Q+1)%2), 1e-8);
                    }
                }
            }
        }
    }


    void TestCompressibleLawsUsingToyCompressibleMaterialLaw()
    {
        double c1 = 2.0;

        ToyCompressibleMaterialLaw<2> law_2d(c1, 0.0, -c1);

        TS_ASSERT_DELTA(law_2d.Get_dW_dI1(2.0,1.0,1.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI3(2.0,1.0,1.0), -c1, 1e-12);

        // compute the stress given C=delta_{MN}, obviously it should be zero
        c_matrix<double,2,2> identity_strain_2d = identity_matrix<double>(2);

        c_matrix<double,2,2> T_2d;
        law_2d.Compute2ndPiolaKirchoffStress(identity_strain_2d,
                                                0.0,
                                                T_2d);
        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                TS_ASSERT_DELTA(T_2d(i,j),0.0,1e-12);
            }
        }


        double c2 = 3.0;

        ToyCompressibleMaterialLaw<3> law_3d(c1, c2, -c1-c2);

        TS_ASSERT_DELTA(law_3d.Get_dW_dI1(2.0,2.0,1.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI2(2.0,2.0,1.0), c2, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI3(2.0,2.0,1.0),  -c1-c2, 1e-12);


        // compute stress given a non-zero deformation
        c_matrix<double,3,3> F;
        F(0,0) = 3.0;
        F(0,1) = 1.0;
        F(1,0) = -1.0;
        F(0,2) = 2.0;
        F(2,0) = 1.0;
        F(1,1) = 6.0;
        F(1,2) = -1.0;
        F(2,1) = 1.5;
        F(2,2) = 0.5;

        c_matrix<double,3,3> C = prod(trans(F),F);

        double I1 =  Trace(C);

        c_matrix<double,3,3> invC = Inverse(C);

        c_matrix<double,3,3> T;
        c_matrix<double,3,3> T2;
        c_matrix<double,3,3> S;
        c_matrix<double,3,3> sigma;

        FourthOrderTensor<3,3,3,3> dTdE;

        law_3d.ComputeStressAndStressDerivative(C, invC, 0.0, T, dTdE, true);
        law_3d.Compute1stPiolaKirchoffStress(F,0.0,S);
        law_3d.Compute2ndPiolaKirchoffStress(C,0.0,T2);
        law_3d.ComputeCauchyStress(F,0.0,sigma);

        c_matrix<double,3,3> FT = prod(F,T);
        c_matrix<double,3,3> F_T_tranF_over_detF = (1.0/Determinant(F))*prod(FT,trans(F));//F*T_as_unsym_tensor*transpose(F);

        c_matrix<double,3,3> T_transposeF = prod(T,trans(F));//T_as_unsym_tensor*transpose(F);


        // check sigma is correct - sigma should be (1/detF) F * T * trans(F)
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(sigma(i,j), F_T_tranF_over_detF(i,j), 1e-12);
            }
        }

        // check S is correct
        for (unsigned M=0; M<3; M++)
        {
            for (unsigned i=0; i<3; i++)
            {
                TS_ASSERT_DELTA(S(M,i), T_transposeF(M,i), 1e-12);
            }
        }

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                // check we gave a symmetric C
                assert(C(M,N)==C(N,M));

                // check the stress
                TS_ASSERT_DELTA(T(M,N), (2*c1+2*c2*I1)*(M==N) - 2*c2*C(M,N) - 2*(c1+c2)*invC(M,N), 1e-12);

                // check alternative computation of the stress
                TS_ASSERT_DELTA(T(M,N), T2(M,N), 1e-12);

                for (unsigned P=0;P<3;P++)
                {
                    for (unsigned Q=0;Q<3;Q++)
                    {
                        double true_val =   4*c2*((M==N)*(P==Q)-(M==P)*(N==Q))
                                            + 4*(c1+c2)*invC(M,P)*invC(Q,N);

                        TS_ASSERT_DELTA(dTdE(M,N,P,Q), true_val, 1e-12);
                    }
                }
            }
        }
    }


    void TestCompressibleMooneyRivlinLaw()
    {
        double c1 = 3.0;
        double c3 = 2.0;

        CompressibleMooneyRivlinMaterialLaw<2> law_2d(c1, c3);

        TS_ASSERT_DELTA(law_2d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.GetC3(), c3, 1e-12);

        TS_ASSERT_DELTA(law_2d.Get_dW_dI1(2.0,2.0,1.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI2(2.0,2.0,1.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI3(2.0,2.0,1.0), -c1, 1e-12);

        double i1 = 2.4;
        double i3 = 0.6;
        TS_ASSERT_DELTA(law_2d.Get_dW_dI1(i1,2.0,i3), c1/sqrt(i3), 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI2(i1,2.0,i3), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI3(i1,2.0,i3), -c1*i1*pow(i3,-1.5)/2 + c3 - c3*pow(i3,-0.5), 1e-12);

        TS_ASSERT_DELTA(law_2d.Get_d2W_dI1I3(i1,2.0,i3), -0.5*c1*pow(i3,-1.5), 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI3(i1,2.0,i3),0.75*c1*i1*pow(i3,-2.5) + 0.5*c3*pow(i3,-1.5), 1e-12);

        TS_ASSERT_DELTA(law_2d.Get_d2W_dI1(i1,2.0,i3), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI2(i1,2.0,i3), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI1I2(i1,2.0,i3), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI2I3(i1,2.0,i3), 0.0, 1e-12);


        law_2d.ScaleMaterialParameters(10);
        TS_ASSERT_DELTA(law_2d.GetC1(), c1/10, 1e-12);
        TS_ASSERT_DELTA(law_2d.GetC3(), c3/10, 1e-12);


        CompressibleMooneyRivlinMaterialLaw<3> law_3d(c1, c3);

        TS_ASSERT_DELTA(law_3d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.GetC3(), c3, 1e-12);

        TS_ASSERT_DELTA(law_3d.Get_dW_dI1(3.0,3.0,1.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI2(3.0,3.0,1.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI3(3.0,3.0,1.0), -c1, 1e-12);
    }
};

#endif /*TESTMATERIALLAWS_HPP_*/
