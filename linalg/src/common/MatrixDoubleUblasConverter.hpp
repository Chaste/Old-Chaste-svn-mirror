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

template <int DIM>
class MatrixDoubleUblasConverter
{
    // This default case has no methods, since it should never be used.
};

template <>
class MatrixDoubleUblasConverter<1>
{
public:
    /***
     * Convert a MatrixDouble to c_matrix<double,DIM,DIM>
     * where DIM==1.
     ***/
    c_matrix<double, 1, 1>* ConvertToUblas(MatrixDouble &matrixDouble)
    {
        return matrixDouble.GetUblasHandle1();
    }
    
    c_matrix<double, 1, 1>& rConvertToUblas(MatrixDouble &matrixDouble)
    {
        return matrixDouble.rGetUblasHandle1();
    }
};

template <>
class MatrixDoubleUblasConverter<2>
{
public:
    /***
     * Convert a MatrixDouble to c_matrix<double,DIM,DIM>
     * where DIM==2.
     ***/
    c_matrix<double, 2, 2>* ConvertToUblas(MatrixDouble &matrixDouble)
    {
        return matrixDouble.GetUblasHandle2();
    }
    
    c_matrix<double, 2, 2>& rConvertToUblas(MatrixDouble &matrixDouble)
    {
        return matrixDouble.rGetUblasHandle2();
    }
};

template <>
class MatrixDoubleUblasConverter<3>
{
public:
    /***
     * Convert a MatrixDouble to c_matrix<double,DIM,DIM>
     * where DIM==3.
     ***/
    c_matrix<double, 3, 3>* ConvertToUblas(MatrixDouble &matrixDouble)
    {
        return matrixDouble.GetUblasHandle3();
    }
    
    c_matrix<double, 3, 3>& rConvertToUblas(MatrixDouble &matrixDouble)
    {
        return matrixDouble.rGetUblasHandle3();
    }
};

template <>
class MatrixDoubleUblasConverter<4>
{
public:
    /***
     * Convert a MatrixDouble to c_matrix<double,DIM,DIM>
     * where DIM==4.
     ***/
    c_matrix<double, 4, 4>* ConvertToUblas(MatrixDouble &matrixDouble)
    {
        return matrixDouble.GetUblasHandle4();
    }
    
    c_matrix<double, 4, 4>& rConvertToUblas(MatrixDouble &matrixDouble)
    {
        return matrixDouble.rGetUblasHandle4();
    }
};
#endif /*MATRIXDOUBLEUBLASCONVERTER_HPP_*/
