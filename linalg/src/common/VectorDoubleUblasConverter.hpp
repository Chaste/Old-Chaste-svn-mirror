#ifndef VECTORDOUBLEUBLASCONVERTER_HPP_
#define VECTORDOUBLEUBLASCONVERTER_HPP_


#include "VectorDouble.hpp"

/**
 * This class converts between VectorDouble and c_vector<double, DIM>
 * where DIM is the size of the vector. i.e. a uBlas vector
 * The key aspect of this class is that it has a template paramter
 * for the size of the vector. It can thus be used to avoid
 * switch statements in the caller as compared to calling VectorDouble
 * directly
 **/

template<int DIM> 
class VectorDoubleUblasConverter
{
public:
    /***
     * Convert a VectorDouble to c_vector<double,DIM>
     * where 1<=DIM<=4.
     ***/
    c_vector<double, DIM> ConvertToUblas(VectorDouble &vectorDouble)
    {
    switch (DIM)
        {
        case 1:
            return vectorDouble.GetUblasHandle1();
            break;
        case 2:
            return vectorDouble.GetUblasHandle2();
            break;
        case 3:
            return vectorDouble.GetUblasHandle3();
            break;      
        case 4:
            return vectorDouble.GetUblasHandle4();
            break;
        default:
            assert(0);
        }          
    }
};

#endif /*VECTORDOUBLEUBLASCONVERTER_HPP_*/
