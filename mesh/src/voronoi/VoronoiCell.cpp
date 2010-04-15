/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#include "VoronoiCell.hpp"


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


bool VoronoiCell::EqualFaces(Face<3>& face1, bool orientation1, Face<3>& face2, bool orientation2)
{
    if (orientation1 == orientation2)
    {
        return (face1 == face2);
    }
    else
    {
        Face<3> face3 =- face2;
        return (face1 == face3);
    }
}

bool VoronoiCell::operator==(VoronoiCell& rOtherCell)
{
    if (mFaces.size() != rOtherCell.GetNumFaces())
    {
        return false;
    }

    std::vector<bool> other_faces_matched;

    std::vector< Face<3>* >::iterator this_face_iterator = mFaces.begin();
    std::vector<bool>::iterator this_orientation_iterator = mOrientations.begin();

    while (this_face_iterator != mFaces.end())
    {
        std::vector< Face<3>* >::iterator other_face_iterator = rOtherCell.mFaces.begin();
        std::vector<bool>::iterator other_orientation_iterator = rOtherCell.mOrientations.begin();

        while ( other_face_iterator != rOtherCell.mFaces.end()
                && !EqualFaces(**this_face_iterator, *this_orientation_iterator,
                                 **other_face_iterator, *other_orientation_iterator) )
        {
            ++other_face_iterator;
            ++other_orientation_iterator;
        }
        if (other_face_iterator == rOtherCell.mFaces.end())
        {
            return false;
        }
        ++this_face_iterator;
        ++this_orientation_iterator;
    }
    return true;
}

c_vector<double, 3>& VoronoiCell::rGetVoronoiCellCentre()
{
    return mCellCentre;
}

unsigned VoronoiCell::GetNumFaces() const
{
    return mFaces.size();
}

const Face<3>& VoronoiCell::rGetFace(unsigned index) const
{
    return *(mFaces[index]);
}

bool VoronoiCell::FaceIsOrientatedClockwise(unsigned index) const
{
    return mOrientations[index];
}

void VoronoiCell::AddFace(Face<3>* pFace)
{
    mFaces.push_back(pFace);
}

void VoronoiCell::AddOrientation(bool isOrientedClockwise)
{
    mOrientations.push_back(isOrientedClockwise);
}

void VoronoiCell::SetCellCentre(c_vector<double, 3> cellCentre)
{
    mCellCentre = cellCentre;
}
