#ifndef FACE_HPP_
#define FACE_HPP_

#include "UblasCustomFunctions.hpp"
#include <vector>

template <unsigned DIM>
class Face
{
public:
    /**
     * The vertices of the face, in clockwise order. Each vertex must be distinct.
     */
    std::vector< c_vector<double, DIM>* > mVertices;

private:    
    void Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                   Face<DIM>& rFace) const;

public:
    /**
     * Compare two faces for equality.
     * Two faces are the same if their vertices differ only by cyclic permutation.
     */
    bool operator==(Face<DIM>& otherFace);
    
    /**
     * Compare two faces for inequality
     */
    bool operator!=(Face<DIM>& otherFace);
    
    /**
     * Return a new face in which the order of the vertices is reversed.
     */
    Face<DIM> operator-();
    
    /**
     * Gets the sum of the length of all edges
     */
    double GetPerimeter();
    
    /**
     * Gets the area of a voronoi face for 2d space only
     */
    double GetArea();
    
};


#endif /*FACE_HPP_*/
