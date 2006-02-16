#ifndef TESTUBLAS_HPP_
#define TESTUBLAS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
class TestUblas : public CxxTest::TestSuite 
{
public:
    
    void TestUblas1( void )
    {
        using namespace boost::numeric::ublas;
        matrix<double> m (3, 3);
        for (unsigned i = 0; i < m.size1 (); ++ i)
            for (unsigned j = 0; j < m.size2 (); ++ j)
                m (i, j) = 3 * i + j;
        // std::cout << m << std::endl;
        
        
    }
    
    void TestUblas2( void )
    {
        using namespace boost::numeric::ublas;
        c_matrix<double,3,3> m (3, 3);
        c_vector<double,3> v (3);
        c_vector<double,3> w (3);
        for (unsigned i = 0; i < std::min (m.size1 (), v.size ()); ++ i) {
            for (unsigned j = 0; j < m.size2 (); ++ j)
                m (i, j) = 3 * i + j;
            v (i) = i;
        }

        TS_ASSERT_EQUALS(v(2),2);

        // run for 5E7 times for a sensible benchmark
        for (int k = 0; k < 5E1; k++) {        
            w = prod (m, v);
            w = prod (v, m);
        }
        
        // std::cout << w << std::endl;
        


    }    
};

#endif /*TESTUBLAS_HPP_*/
