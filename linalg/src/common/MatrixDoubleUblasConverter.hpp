#ifndef MATRIXDOUBLEUBLASCONVERTER_HPP_
#define MATRIXDOUBLEUBLASCONVERTER_HPP_

#include "MatrixDouble.hpp"

/**
 * This class converts between MatrixDouble and c_matrix<double, DIM, DIM>
 * where DIM is the size of the square matrix. i.e. a uBlas matrix
 * The key aspect of this class is that it has a template paramter
 * for the size of the matrix. It can thus be used to avoid
 * switch statements in the caller as compared to calling MatrixDouble
 * directly
 **/

template<int DIM> 
class MatrixDoubleUblasConverter
{
public:
    /***
     * Convert a MatrixDouble to c_matrix<double,DIM,DIM>
     * where 1<=DIM<=4.
     ***/
    c_matrix<double, DIM, DIM> ConvertToUblas(MatrixDouble &matrixDouble)
    {
    switch (DIM)
        {
        case 1:
            return matrixDouble.GetUblasHandle1();
            break;
        case 2:
            return matrixDouble.GetUblasHandle2();
            break;
        case 3:
            return matrixDouble.GetUblasHandle3();
            break;      
        case 4:
            return matrixDouble.GetUblasHandle4();
            break;
        default:
            assert(0);
        }          
    }
};


#endif /*MATRIXDOUBLEUBLASCONVERTER_HPP_*/
