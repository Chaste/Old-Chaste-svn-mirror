#ifndef TESTVECTORDOUBLEUBLASCONVERTER_
#define TESTVECTORDOUBLEUBLASCONVERTER_


#include "VectorDouble.hpp"
#include "VectorDoubleUblasConverter.hpp"

class TestVectorDoubleUblasConverter : public CxxTest::TestSuite
{
public:
	
	void TestConvertor()
    {
        VectorDoubleUblasConverter<2> converter2;
        VectorDouble a(2);
        a(0)=3; a(1)=4;

        c_vector<double, 2>& a_ublas_vector = converter2.rConvertToUblas(a);
        TS_ASSERT_EQUALS(a(0), a_ublas_vector(0));
        TS_ASSERT_EQUALS(a(1), a_ublas_vector(1));
        
        VectorDoubleUblasConverter<3> converter3;
        VectorDouble b(3);
        b(0)=3; b(1)=4; b(2)=10;

        c_vector<double, 3>& b_ublas_vector = converter3.rConvertToUblas(b);
        TS_ASSERT_EQUALS(b(0), b_ublas_vector(0));
        TS_ASSERT_EQUALS(b(1), b_ublas_vector(1));
        TS_ASSERT_EQUALS(b(2), b_ublas_vector(2));   
    }

    void TestSetWithConvertor()
    {
        VectorDoubleUblasConverter<1> converter1;
        VectorDouble a(1);
        a(0)=3; 
 
        c_vector<double, 1>& a_ublas_vector = converter1.rConvertToUblas(a);
            
        // set via ublas handle
        a_ublas_vector(0)=5;
        
        // test initial vector has changed
        TS_ASSERT_EQUALS( a(0), 5 );
   }
};
#endif /*TESTVECTORDOUBLEUBLASCONVERTER_*/
