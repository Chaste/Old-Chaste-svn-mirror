#ifndef _BASISFUNCTIONSCHECKERS_HPP_
#define _BASISFUNCTIONSCHECKERS_HPP_
#include <cxxtest/TestSuite.h>
#include "LinearBasisFunction.hpp"
#include <vector>

template <int SPACE_DIM>
class BasisFunctionsCheckers
{
public:
    void checkBasisFunctions(std::vector<Point<SPACE_DIM>*> evaluationPoints)
    {
        int size = evaluationPoints.size();		// number of evalutation points and basis functions too
        std::vector<double> basis_function_vector; // store results of evalutation
        
        double expected_evaluation;
        
        for (int point_index=0; point_index<size; point_index++)
        {
            c_vector<double, SPACE_DIM+1> basis_function_vector
                = LinearBasisFunction<SPACE_DIM>::ComputeBasisFunctions(*(evaluationPoints[point_index]));
            for (int func_index=0; func_index<size; func_index ++)
            {
                if (func_index==point_index)
                {
                    expected_evaluation=1.0;
                }
                else
                {
                    expected_evaluation = 0.0;
                }
                TS_ASSERT_DELTA(basis_function_vector(func_index),
                                expected_evaluation,
                                1e-12);
                                
                                
                TS_ASSERT_DELTA(LinearBasisFunction<SPACE_DIM>::ComputeBasisFunction(*(evaluationPoints[point_index]),func_index),
                                expected_evaluation,
                                1e-12);
            }
        }
    }
};
#endif //_BASISFUNCTIONSCHECKERS_HPP_
