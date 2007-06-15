#include "VoronoiCell.hpp"

VoronoiCell::VoronoiCell(c_vector<double, 3> cell_centre, 
                         std::set< Node<3>* > vertices, 
                         std::vector< std::vector < unsigned > > faces, 
                         unsigned colour)
    : mCellCentre(cell_centre),
      mpVertices(vertices),
      mFaces(faces),
      mColour(colour)
{   
}

std::set< Node<3>* >& VoronoiCell::GetVertices()
{
    return mpVertices;
}
