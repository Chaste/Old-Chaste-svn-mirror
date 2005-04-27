#ifndef _TESTNONNECROTICTUMOURNONLINEARELLIPTICEQUATION_HPP_
#define _TESTNONNECROTICTUMOURNONLINEARELLIPTICEQUATION_HPP_

#include "EquationForNonnecroticTumour.hpp"
#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"

/**
 *  Test for the non-necrotic Tumour equation which is a linear elliptic equation.
 */
class TestNonnecroticTumourNonlinearEllipticEquation : public CxxTest::TestSuite
{
public:
    void testNonnecroticTumourC()
    {
       
        Point<1> zero1(0);
        Point<2> zero2(0,0);
        Point<3> zero3(0,0,0);
        double u = 2.0;
        
        EquationForNonnecroticTumour<1> tumour_equ1;
        EquationForNonnecroticTumour<2> tumour_equ2;
        EquationForNonnecroticTumour<3> tumour_equ3;
        
        TS_ASSERT_DELTA(tumour_equ1.ComputeLinearSourceTerm(zero1), 0.0, 1e-12);
        TS_ASSERT_DELTA(tumour_equ2.ComputeLinearSourceTerm(zero2), 0.0, 1e-12);
        TS_ASSERT_DELTA(tumour_equ3.ComputeLinearSourceTerm(zero3), 0.0, 1e-12);
        
        TS_ASSERT_DELTA(tumour_equ1.ComputeNonlinearSourceTerm(zero1,u),(double)2/3,1e-12);
        TS_ASSERT_DELTA(tumour_equ2.ComputeNonlinearSourceTerm(zero2,u),(double)2/3,1e-12);
        TS_ASSERT_DELTA(tumour_equ3.ComputeNonlinearSourceTerm(zero3,u),(double)2/3,1e-12);        
        
        MatrixDouble diff1 = tumour_equ1.ComputeDiffusionTerm(zero1,u);
        MatrixDouble diff2 = tumour_equ2.ComputeDiffusionTerm(zero2,u);
        MatrixDouble diff3 = tumour_equ3.ComputeDiffusionTerm(zero3,u);
        
        TS_ASSERT_DELTA(diff1(0,0),1,1e-12);

        TS_ASSERT_DELTA(diff2(0,0),1,1e-12);
        TS_ASSERT_DELTA(diff2(1,1),1,1e-12);
        TS_ASSERT_DELTA(diff2(0,1),0,1e-12);

        TS_ASSERT_DELTA(diff3(0,0),1,1e-12);
        TS_ASSERT_DELTA(diff3(1,1),1,1e-12);
        TS_ASSERT_DELTA(diff3(2,2),1,1e-12);
        TS_ASSERT_DELTA(diff3(0,1),0,1e-12);
        TS_ASSERT_DELTA(diff3(0,2),0,1e-12);
        TS_ASSERT_DELTA(diff3(1,2),0,1e-12);            

    }
};


#endif //_TESTNONNECROTICTUMOURNONLINEARELLIPTICEQUATION_HPP_
