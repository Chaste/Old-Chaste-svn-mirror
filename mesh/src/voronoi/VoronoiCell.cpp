#include "VoronoiCell.hpp"

VoronoiCell::VoronoiCell(c_vector<double, 3> cell_centre, 
                         std::vector< c_vector<double, 3> > vertices, 
                         std::vector< std::vector < unsigned > > faces, 
                         unsigned colour)
    : mCellCentre(cell_centre),
      mVertices(vertices),
      mFaces(faces),
      mColour(colour)
{   
}
