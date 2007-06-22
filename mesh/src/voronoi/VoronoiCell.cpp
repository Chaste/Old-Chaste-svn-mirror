#include "VoronoiCell.hpp"

bool VoronoiCell::EqualFaces(Face& face1, bool orientation1, Face& face2, bool orientation2)
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
};

bool VoronoiCell::operator==(VoronoiCell& otherCell)
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
};
