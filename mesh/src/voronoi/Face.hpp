#ifndef FACE_HPP_
#define FACE_HPP_

#include "UblasCustomFunctions.hpp"
#include <vector>

class Face
{
public:
    /**
     * The vertices of the face, in clockwise order. Each vertex must be distinct.
     */
    std::vector< c_vector<double, 3>* > mVertices;

private:    
    const void Increment(std::vector< c_vector<double, 3>* >::iterator& rIterator,
                         Face& rFace);

public:
    /**
     * Compare two faces for equality.
     * Two faces are the same if their vertices differ only by cyclic permutation.
     */
    bool operator==(Face& otherFace);
    
    /**
     * Compare two faces for inequality
     */
    bool operator!=(Face& otherFace);
    
    /**
     * Return a new face in which the order of the vertices is reversed.
     */
    Face operator-();
};


#endif /*FACE_HPP_*/
