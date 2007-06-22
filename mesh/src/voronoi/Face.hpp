#ifndef FACE_HPP_
#define FACE_HPP_

#include "UblasCustomFunctions.hpp"
#include <vector>

class Face
{
public:
    std::vector< c_vector<double, 3>* > mVertices;

private:    
    const void Increment(std::vector< c_vector<double, 3>* >::iterator& rIterator,
                   Face& rFace)
    {
        rIterator++;
        if (rIterator==rFace.mVertices.end() )
        {
            rIterator=rFace.mVertices.begin();
        }
    }

public:    
    bool operator==(Face& otherFace)
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
    }
    
     bool operator!=(Face& otherFace)
     {
        return !(*this==otherFace);
     }
     
     Face operator-()
     {
        Face reversed_face;
        std::vector< c_vector<double, 3>* >::iterator this_iterator=mVertices.end();
        while (this_iterator !=mVertices.begin())
        {
            this_iterator--;
            reversed_face.mVertices.push_back(*this_iterator);
        }
        return reversed_face;
     }
    
};


#endif /*FACE_HPP_*/
