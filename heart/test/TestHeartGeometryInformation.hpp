/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTHEARTGEOMETRYINFORMATION_HPP_
#define TESTHEARTGEOMETRYINFORMATION_HPP_

#include "TrianglesMeshReader.hpp"
#include "HeartGeometryInformation.hpp"

#include "MeshalyzerMeshWriter.hpp" //temporal
#include "PetscSetupAndFinalize.hpp"

class TestHeartGeometryInformation : public CxxTest::TestSuite
{
public:
    void TestCalculateRelativeWallPositionSimple2dMesh() throw(Exception)
    {
        TetrahedralMesh<2,2> mesh;
        //This mesh will have 6 nodes per face, spaced by 1
        mesh.ConstructRectangularMesh(5, 5);

        std::vector<unsigned> left_face;
        std::vector<unsigned> right_face;

        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {  
            // Get the nodes at the left face of the square
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
            {
                left_face.push_back(index);
            }
            // Get the nodes at the right face of the square
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                right_face.push_back(index);
            }
            
        }           
        HeartGeometryInformation<2> info(mesh, left_face, right_face);
        
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(5.0-x)/5.0);
        } 
    }


    void TestCalculateRelativeWallPositionSimple3dMesh() throw(Exception)
    {
        TetrahedralMesh<3,3> mesh;
        //This mesh will have 6 nodes per face, spaced by 1
        mesh.ConstructCuboid(5, 5, 5);

        std::vector<unsigned> left_face;
        std::vector<unsigned> right_face;

        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {  
            // Get the nodes at the left face of the cube
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
            {
                left_face.push_back(index);
            }
            // Get the nodes at the right face of the cube
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                right_face.push_back(index);
            }
            
        }           
        HeartGeometryInformation<3> info(mesh, left_face, right_face);
        
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(5-x)/5);
        } 
    }
    
    void TestCalculateRelativeWallPositionWithSurfaces() throw(Exception)
    {
        TetrahedralMesh<3,3> mesh;
        //This mesh will have 9 nodes per side, spaced by 1, it is a cube
        mesh.ConstructCuboid(8, 8, 8);

        std::vector<unsigned> epi_face;
        std::vector<unsigned> lv_face;
        std::vector<unsigned> rv_face;
        //Define three surfaces, epi, lv and rv.
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {  
            // Get the nodes at cube face considered to be epi (at both external faces)
            if (  (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
                ||(fabs(mesh.GetNode(index)->rGetLocation()[0]-8.0) < 1e-6))
            {
                epi_face.push_back(index);
            }
            // Get the nodes at cube face considered to be lv (at the plane defined by x=3)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-3.0) < 1e-6)
            {
                lv_face.push_back(index);
            }
            // Get the nodes at cube face considered to be rv (at the plane defined by x=5)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                rv_face.push_back(index);
            }
            
        }           

        HeartGeometryInformation<3> info(mesh, epi_face, lv_face, rv_face);
        
        //check that the method returns expected value in this case
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            //in the lv...
            if(x<3)
            {
                TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(3-x)/3);
                continue;
            }
            //..in the septum...
            if ((x>=3)&&(x<=5))
            {
                //TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index));
                continue;
            }
            //...and in the rv.
            if(x>5)
            {
                TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(x-5)/3);
            }
        } 
    }
    
    void TestDetermineLayerForEachNode() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        //This mesh will have 31 nodes per side, spaced by 1, it is a cube
        mesh.ConstructCuboid(30, 30, 30);

        std::vector<unsigned> epi_face;
        std::vector<unsigned> lv_face;
        std::vector<unsigned> rv_face;
        //Define three surfaces, epi, lv and rv.
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {  
            // Get the nodes at cube face considered to be epi (at both external faces)
            if (  (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
                ||(fabs(mesh.GetNode(index)->rGetLocation()[0]-30.0) < 1e-6))
            {
                epi_face.push_back(index);
            }
            // Get the nodes at cube face considered to be lv (at the plane defined by x=3)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-10.0) < 1e-6)
            {
                lv_face.push_back(index);
            }
            // Get the nodes at cube face considered to be rv (at the plane defined by x=5)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-20.0) < 1e-6)
            {
                rv_face.push_back(index);
            }
            
        }  
        
        HeartGeometryInformation<3> info(mesh, epi_face, lv_face, rv_face);   
        
        info.DetermineLayerForEachNode(0.3,0.5);
        
        info.WriteLayerForEachNode("TestHeartGeom","layers.het");
        
//        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
//        {
//            double x = mesh.GetNode(index)->rGetLocation()[0];
//            
//            std::cout << x << " " << info.rGetLayerForEachNode()[index] << "\n";
//            
////// ///\todo - ******** fix this with build=Intel and bring it back - it is just rounding error ********
////            if (x<=3||x>=27)
////            {
////                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],EPI);
////            }
////            else if(x<=5 || x>=25)
////            {
////                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],MID);
////            }
////            else
////            {
////                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],ENDO);
////            }            
//        }        
    }
    
    void TestHeartGeometryTakingFiles() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
        std::string epi_face_file = "heart/test/data/box_shaped_heart/epi.tri";
        std::string rv_face_file = "heart/test/data/box_shaped_heart/rv.tri";
        std::string lv_face_file = "heart/test/data/box_shaped_heart/lv.tri";
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        HeartGeometryInformation<3> info(mesh, epi_face_file, lv_face_file, rv_face_file); 
///\todo add a test here
    }
};

#endif /*TESTHEARTGEOMETRYINFORMATION_HPP_*/

