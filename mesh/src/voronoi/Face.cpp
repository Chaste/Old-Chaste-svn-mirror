#include "Face.hpp"

template <unsigned DIM>
void Face<DIM>::Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                     Face<DIM>& rFace) const
{
    rIterator++;
    if (rIterator==rFace.mVertices.end() )
    {
        rIterator=rFace.mVertices.begin();
    }
};

template <unsigned DIM>
bool Face<DIM>::operator==(Face<DIM>& otherFace)
{
    typename std::vector< c_vector<double, DIM>* >::iterator this_iterator=mVertices.begin();
    typename std::vector< c_vector<double, DIM>* >::iterator other_iterator=otherFace.mVertices.begin();
    // find first vertex
    while ( this_iterator!=mVertices.end() && 
            other_iterator!=otherFace.mVertices.end() && 
            norm_2(**this_iterator - **other_iterator) >1e-10 )
    {
        this_iterator++;
    }
    if (this_iterator==mVertices.end() || other_iterator==otherFace.mVertices.end())
    {
        // could not find first vertex
        // faces are distinct unless they are empty
        return this_iterator==mVertices.end() && other_iterator==otherFace.mVertices.end();
    }
    
    typename std::vector< c_vector<double, DIM>* >::iterator this_start=this_iterator;
    Increment(this_iterator, *this);
    Increment(other_iterator, otherFace);
    
    // check remanining vertices are equal
    while (this_iterator!=this_start)
    {
        if (norm_2(**this_iterator - **other_iterator) >1e-10)
        {
            return false;
        }
        else
        {
            Increment(this_iterator, *this);
            Increment(other_iterator, otherFace);
        }
    }
    return other_iterator==otherFace.mVertices.begin();
};

#define COVERAGE_IGNORE //Spuriously not covered
template <unsigned DIM>
bool Face<DIM>::operator!=(Face& otherFace)
{
   return !(*this==otherFace);
};
#undef COVERAGE_IGNORE

template <unsigned DIM>
Face<DIM> Face<DIM>::operator-()
{
   Face<DIM> reversed_face;
   typename std::vector< c_vector<double, DIM>* >::iterator this_iterator=mVertices.end();
   while (this_iterator !=mVertices.begin())
   {
       this_iterator--;
       reversed_face.mVertices.push_back(*this_iterator);
   }
   return reversed_face;
};

template <unsigned DIM>
double Face<DIM>::GetPerimeter() const
{
    double perimeter_return = 0;
    for(unsigned i=0; i<mVertices.size(); i++)
    {
        perimeter_return += norm_2(*mVertices[i]-*mVertices[(i+1)%mVertices.size()]);
    }
    return perimeter_return;
};

template <unsigned DIM>
double Face<DIM>::GetArea() const
{
    assert(DIM==2);
    double area_return = 0;
    for(unsigned i=0; i<mVertices.size(); i++)
    {
        //  Area = sum ( x_i * y_i+1 - y_i * x_i+1 )/2.0 over all vertices, 
        //      assuming vertices are ordered anti-clockwise
        area_return +=   ( (*mVertices[i])(0) * (*mVertices[(i+1)%mVertices.size()])(1) 
                          -(*mVertices[i])(1) * (*mVertices[(i+1)%mVertices.size()])(0) ) / 2.0 ;
    }
    return area_return;
};


