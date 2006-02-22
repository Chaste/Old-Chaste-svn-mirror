#include "VectorUblas.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cassert>
#include <iostream>
#include <math.h>

//VectorUblas::VectorUblas()
//{
//    mSize = -1;
//}

VectorUblas::VectorUblas(int Size)
{
    assert(Size > 0);
    using namespace boost::numeric::ublas;
    
    mSize = Size;
    
    switch (Size)
    {
        case 1:
            mpVectorOf1 = new c_vector<double,1>(1); 
            break;
        case 2:
            mpVectorOf2 = new c_vector<double,2> (2); 
            break; 
        case 3:
            mpVectorOf3 = new c_vector<double,3> (3); 
            break; 
        
        default:
        // Vector biger than  3 Throw Exception
        throw("Vector size larger than 3");
            break;
    }
}

VectorUblas::~VectorUblas()
{
    switch (mSize)
    {
        case 1:
            delete (mpVectorOf1);
            break;
        case 2:
            delete (mpVectorOf2);
            break;
        case 3:
            delete (mpVectorOf3);
            break;
        default:
        // Vector biger than  3 Throw Exception
        throw("Vector size larger than 3");
            break;
    }        
}    
    
double VectorUblas::operator()(int entry) const
{
    assert(entry > -1);
    assert(entry < mSize);
    
    switch(mSize) {
        case 1:
            return (*mpVectorOf1)(entry);
        case 2:
            return (*mpVectorOf2)(entry); 
        case 3:
            return (*mpVectorOf3)(entry);   
        default:
            // Vector biger than  3 Throw Exception
            throw("Vector size larger than 3");
    }

}

    
           
   

