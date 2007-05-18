#ifndef TESTVORONOICELL_HPP_
#define TESTVORONOICELL_HPP_

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"

#include <cmath>
#include <vector>

class TestVoronoiCell : public CxxTest::TestSuite
{
public:
    void TestCreateCell()
    {
        //Test for a 1x1x1 cube centered at (0.5,0.5,0.5)  
        
        c_vector<double, 3> cell_centre;
        cell_centre(0)=0.5;
        cell_centre(1)=0.5;
        cell_centre(2)=0.5;
        
        std::vector< c_vector<double, 3> > vertices;
        for (unsigned i=0; i<=1; i++)
        {
            for (unsigned j=0; i<=1; i++)
            {
                for (unsigned k=0; i<=1; i++)
                {
                    c_vector<double, 3> vertex;
                    vertex(0)=(double) i;
                    vertex(1)=(double) j;
                    vertex(2)=(double) k;
                    vertices.push_back(vertex);
                }
            }
        }
        
        std::vector< std::vector < unsigned > > faces;
        
        std::vector<unsigned> empty;
        faces.push_back(empty);
        faces[0].push_back(0);
        faces[0].push_back(1);
        faces[0].push_back(5);
        faces[0].push_back(4);
        
        faces.push_back(empty);
        faces[1].push_back(4);
        faces[1].push_back(5);
        faces[1].push_back(7);
        faces[1].push_back(6);
        
        faces.push_back(empty);
        faces[2].push_back(6);
        faces[2].push_back(7);
        faces[2].push_back(3);
        faces[2].push_back(2);
        
        faces.push_back(empty);
        faces[3].push_back(2);
        faces[3].push_back(3);
        faces[3].push_back(1);
        faces[3].push_back(0);
        
        faces.push_back(empty);
        faces[4].push_back(1);
        faces[4].push_back(5);
        faces[4].push_back(7);
        faces[4].push_back(3);
        
        faces.push_back(empty);
        faces[5].push_back(0);
        faces[5].push_back(4);
        faces[5].push_back(6);
        faces[5].push_back(2);
        
        unsigned colour = 1u;
        
        VoronoiCell cell(cell_centre, vertices, faces, colour);
        
//        c_vector<double, 3> return_cell_centre = cell.GetCellCentre();
//        std::vector< c_vector<double, 3> > return_vertices = cell.GetVertices();
//        std::vector< std::vector < unsigned > > return_faces = cell.GetFaces();
//        unsigned return_colour = cell.GetColour();
    }
    


};


#endif /*TESTVORONOICELL_HPP_*/
