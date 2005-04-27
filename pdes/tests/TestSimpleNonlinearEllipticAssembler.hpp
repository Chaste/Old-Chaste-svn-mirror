#ifndef _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_

#include "SimpleNonlinearEllipticAssembler.hpp"
#include <cxxtest/TestSuite.h>
//#include "petscvec.h"
//#include "petscmat.h"
  
class TestSimpleNonlinearEllipticAssembler : public CxxTest::TestSuite 
{
	
public:
    void testASimpleNonlinearEllipticAssembler( void )
    {
        //create a new SimpleNonlinearEllipticAssembler
        SimpleNonlinearEllipticAssembler<1,1> assembler;
    }
        
};

#endif //_TESTSIMPLENONLINEARELLIPTICASSEMBLER_HPP_
