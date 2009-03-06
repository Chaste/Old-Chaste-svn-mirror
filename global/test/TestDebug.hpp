#ifndef TESTDEBUG_HPP_
#define TESTDEBUG_HPP_

#include <cxxtest/TestSuite.h>
#include "Debug.hpp"

class TestDebug : public CxxTest::TestSuite
{
public:
    // Can't really test these other than that they compile and visually looking at output.
    void TestDebugMacros()
    {
        TRACE("Some trace");
        
        unsigned my_var = 3141;
        PRINT_VARIABLE(my_var);
        
        double another_var = 2.81;
        PRINT_VARIABLES(my_var, another_var);
        
        double cancer_curing_constant = 0.053450242435;
        PRINT_3_VARIABLES(my_var, another_var, cancer_curing_constant);

        double heart_disease_ending_constant = -3e-141;
        PRINT_4_VARIABLES(my_var, another_var, cancer_curing_constant, heart_disease_ending_constant);
        
        for(unsigned i=0; i<10; i++)
        {
            HOW_MANY_TIMES_HERE("inside for loop");
            
            for(unsigned j=0; j<2; j++)
            {
                HOW_MANY_TIMES_HERE("nested loop");
            }
        }
        
        for(unsigned j=0; j<10 /*change to 11 and it should quit*/; j++)
        {
            QUIT_AFTER_N_VISITS(11);
        }
    }
};

#endif /*TESTDEBUG_HPP_*/
