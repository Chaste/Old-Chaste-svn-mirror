#include "Face.hpp"

const void Face::Increment(std::vector< c_vector<double, 3>* >::iterator& rIterator,
                           Face& rFace)
{
    rIterator++;
    if (rIterator==rFace.mVertices.end() )
    {
        rIterator=rFace.mVertices.begin();
    }
};

bool Face::operator==(Face& otherFace)
{
    std::vector< c_vector<double, 3>* >::iterator this_iterator=mVertices.begin();
    std::vector< c_vector<double, 3>* >::iterator other_iterator=otherFace.mVertices.begin();
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
    
    std::vector< c_vector<double, 3>* >::iterator this_start=this_iterator;
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

bool Face::operator!=(Face& otherFace)
{
   return !(*this==otherFace);
};
 
Face Face::operator-()
{
   Face reversed_face;
   std::vector< c_vector<double, 3>* >::iterator this_iterator=mVertices.end();
   while (this_iterator !=mVertices.begin())
   {
       this_iterator--;
       reversed_face.mVertices.push_back(*this_iterator);
   }
   return reversed_face;
};
