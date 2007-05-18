#ifndef VORONOICELL_HPP_
#define VORONOICELL_HPP_

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>

#include <cmath>
#include <vector>

class VoronoiCell
{
private:
    c_vector<double, 3> mCellCentre;
    std::vector< c_vector<double, 3> > mVertices;
    std::vector< std::vector < unsigned > > mFaces;
    unsigned mColour;
    
public:
    VoronoiCell(c_vector<double, 3> cell_centre, 
                std::vector< c_vector<double, 3> > vertices, 
                std::vector< std::vector < unsigned > > faces, 
                unsigned colour = 1);
                
    
    
};

#endif /*VORONOICELL_HPP_*/
