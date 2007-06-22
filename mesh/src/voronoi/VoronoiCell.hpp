#ifndef VORONOICELL_HPP_
#define VORONOICELL_HPP_

#include "UblasCustomFunctions.hpp"
#include "Face.hpp"
#include <cxxtest/TestSuite.h>

#include <cmath>
#include <vector>

class VoronoiCell
{
public:

    /***
     * Faces of the cell, which should be distinct.
     */
    std::vector< Face* > mFaces;
    
    /***
     * How each face is oriented.
     * From the perspective of the centre of the cell, the vertices of each face should be ordered clockwise.
     * If and only if this is false, the order of vertices in the corresponding face should be reversed.
     * 
     * N.B. Most faces belong to two cells, but with opposite orientations. This allows us to reuse the face data
     * across the two cells.
     */
    std::vector< bool > mOrientations;

private:
    bool EqualFaces(Face& face1, bool orientation1, Face& face2, bool orientation2);

public:

    /***
     * Test whether two cells are equal.
     * 
     * Two cells are equal if their set of faces are equal (including whether the faces have the same orientations).
     */
    bool operator==(VoronoiCell& otherCell);
};

#endif /*VORONOICELL_HPP_*/
