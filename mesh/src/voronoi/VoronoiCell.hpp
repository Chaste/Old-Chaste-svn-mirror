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
    std::vector< Face* > mFaces; // all faces should be distinct
    std::vector< bool > mOrientations;
    
public:
    bool EqualFaces(Face& face1, bool orientation1, Face& face2, bool orientation2)
    {
        if ( orientation1 == orientation2)
        {
            return face1==face2;
        }
        else
        {
            Face face3=-face2;
            return face1==face3;
        }
    }

    bool operator==(VoronoiCell& otherCell)
    {
        if ( mFaces.size() != otherCell.mFaces.size() )
        {
            return false;
        }
        
        std::vector< bool > other_faces_matched;
        
        std::vector< Face* >::iterator this_face_iterator=mFaces.begin();
        std::vector< bool >::iterator this_orientation_iterator=mOrientations.begin();
        
        while (this_face_iterator!=mFaces.end())
        {
            std::vector< Face* >::iterator other_face_iterator=otherCell.mFaces.begin();
            std::vector< bool >::iterator other_orientation_iterator=otherCell.mOrientations.begin();
            while ( other_face_iterator != otherCell.mFaces.end()
                    && !EqualFaces(**this_face_iterator, *this_orientation_iterator,
                                     **other_face_iterator, *other_orientation_iterator) )
            {
                other_face_iterator++;
                other_orientation_iterator++;
            }
            if (other_face_iterator == otherCell.mFaces.end())
            {
                return false;
            }
            this_face_iterator++;
            this_orientation_iterator++;
        }
        return true;
    }
    
};

#endif /*VORONOICELL_HPP_*/
