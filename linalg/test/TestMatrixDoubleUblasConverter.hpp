#ifndef TESTMATRIXDOUBLEUBLASCONVERTER_HPP_
#define TESTMATRIXDOUBLEUBLASCONVERTER_HPP_

#include "MatrixDouble.hpp"
#include "MatrixDoubleUblasConverter.hpp"

class TestMatrixDoubleUblasConverter : public CxxTest::TestSuite
{
    public:
    
    void TestConvertor()
    {
        MatrixDoubleUblasConverter<2> converter2;
        MatrixDouble a(2,2);
        a(0,0)=3; a(0,1)=4; a(1,0)=2; a(1,1)=5; 
        c_matrix<double, 2,2> a_ublas_matrix = converter2.ConvertToUblas(a);
        TS_ASSERT_EQUALS(a(0,0), a_ublas_matrix(0,0));
        TS_ASSERT_EQUALS(a(0,1), a_ublas_matrix(0,1));
        TS_ASSERT_EQUALS(a(1,0), a_ublas_matrix(1,0));
        TS_ASSERT_EQUALS(a(1,1), a_ublas_matrix(1,1));
        
        MatrixDoubleUblasConverter<3> converter3;
        MatrixDouble b(3,3);
        b(0,0)=3; b(2,1)=4; b(1,2)=2; 
        c_matrix<double, 3,3> b_ublas_matrix = converter3.ConvertToUblas(b);
        TS_ASSERT_EQUALS(b(0,0), b_ublas_matrix(0,0));
        TS_ASSERT_EQUALS(b(2,1), b_ublas_matrix(2,1));
        TS_ASSERT_EQUALS(b(1,2), b_ublas_matrix(1,2));
    }

};        
        
#endif /*TESTMATRIXDOUBLEUBLASCONVERTER_HPP_*/
