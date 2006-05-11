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
};

template<> 
class VectorDoubleUblasConverter<1>
{
public:
    /**
     * Convert a VectorDouble to c_vector<double,1>
     **/
    c_vector<double, 1>* ConvertToUblas(VectorDouble &vectorDouble)
    {
        return vectorDouble.GetUblasHandle1();
    }

    c_vector<double, 1>& rConvertToUblas(VectorDouble &vectorDouble)
    {
        return vectorDouble.rGetUblasHandle1();
    }

};

template<> 
class VectorDoubleUblasConverter<2>
{
public:
    /**
     * Convert a VectorDouble to c_vector<double,2>
     **/
    c_vector<double, 2>* ConvertToUblas(VectorDouble &vectorDouble)
    {
        return vectorDouble.GetUblasHandle2();
    }

    c_vector<double, 2>& rConvertToUblas(VectorDouble &vectorDouble)
    {
        return vectorDouble.rGetUblasHandle2();
    }

};

template<> 
class VectorDoubleUblasConverter<3>
{
public:
    /**
     * Convert a VectorDouble to c_vector<double,3>
     **/
    c_vector<double, 3>* ConvertToUblas(VectorDouble &vectorDouble)
    {
        return vectorDouble.GetUblasHandle3();
    }

    c_vector<double, 3>& rConvertToUblas(VectorDouble &vectorDouble)
    {
        return vectorDouble.rGetUblasHandle3();
    }

};

template<> 
class VectorDoubleUblasConverter<4>
{
public:
    /**
     * Convert a VectorDouble to c_vector<double,4>
     **/
    c_vector<double, 4>* ConvertToUblas(VectorDouble &vectorDouble)
    {
        return vectorDouble.GetUblasHandle4();
    }

    c_vector<double, 4>& rConvertToUblas(VectorDouble &vectorDouble)
    {
        return vectorDouble.rGetUblasHandle4();
    }

};
#endif /*VECTORDOUBLEUBLASCONVERTER_HPP_*/
