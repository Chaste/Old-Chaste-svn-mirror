/*

Copyright (C) University of Oxford, 2008

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


#ifndef TESTVORONOITESSELLATION_HPP_
#define TESTVORONOITESSELLATION_HPP_


#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"
#include "VoronoiTessellation.hpp"
#include "RefinableMesh.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include "TrianglesMeshWriter.hpp"

#include <cmath>
#include <vector>

class TestVoronoiTessellation : public CxxTest::TestSuite
{
public:
    void TestReturnPolarAngle() throw (Exception)
    {
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0, true,  1.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(1, true, -1.0, -1.0,  1.0));
        nodes.push_back(new Node<3>(2, true, -1.0,  1.0, -1.0));
        nodes.push_back(new Node<3>(3, true,  1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(4, false, 0.0,0.0,0.0));

        RefinableMesh<3,3> mesh(nodes);

        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);

        // Four cases to test:
        // x> 0, y>0
        double angle = tessellation.ReturnPolarAngle(1.0,sqrt(3.0));
        TS_ASSERT_DELTA(angle, M_PI/3.0, 1e-7);

        // x> 0, y<0
        angle = tessellation.ReturnPolarAngle(1.0,-sqrt(3.0));
        TS_ASSERT_DELTA(angle, -M_PI/3.0, 1e-7);

        // x< 0, y>0
        angle = tessellation.ReturnPolarAngle(-1.0,sqrt(3.0));
        TS_ASSERT_DELTA(angle, 2.0*M_PI/3.0, 1e-7);

        // x< 0, y<0
        angle = tessellation.ReturnPolarAngle(-1.0,-sqrt(3.0));
        TS_ASSERT_DELTA(angle, -2.0*M_PI/3.0, 1e-7);

        // check boundary cases
        angle = tessellation.ReturnPolarAngle(1.0,0.0);
        TS_ASSERT_DELTA(angle, 0.0, 1e-7);

        angle = tessellation.ReturnPolarAngle(0.0,1.0);
        TS_ASSERT_DELTA(angle, M_PI/2.0, 1e-7);

        angle = tessellation.ReturnPolarAngle(-1.0,0.0);
        TS_ASSERT_DELTA(angle, M_PI, 1e-7);

        angle = tessellation.ReturnPolarAngle(0.0,-1.0);
        TS_ASSERT_DELTA(angle, -M_PI/2.0, 1e-7);

        TS_ASSERT_THROWS_ANYTHING(tessellation.ReturnPolarAngle(0.0, 0.0));
    }

    void TestGenerateVerticesFromElementCircumcentres() throw (Exception)
    {
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0, true,  1.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(1, true, -1.0, -1.0,  1.0));
        nodes.push_back(new Node<3>(2, true, -1.0,  1.0, -1.0));
        nodes.push_back(new Node<3>(3, true,  1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(4, false, 0.0,0.0,0.0));

        RefinableMesh<3,3> mesh(nodes);

        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);

        tessellation.GenerateVerticesFromElementCircumcentres();

        c_vector<double,3> this_vertex = *(tessellation.mVertices[0]);

        TS_ASSERT_DELTA(this_vertex[0], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);

        this_vertex = *(tessellation.mVertices[1]);

        TS_ASSERT_DELTA(this_vertex[0], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], 1.5, 1e-7);

        this_vertex = *(tessellation.mVertices[2]);

        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], 1.5, 1e-7);

        this_vertex = *(tessellation.mVertices[3]);

        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);
    }

    void TestSimpleTessellation() throw (Exception)
    {
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;

        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,0.5,0.5));

        RefinableMesh<3,3> mesh(nodes);
        //TrianglesMeshWriter<3,3> mesh_writer("","Simple_Delaunay_Tet");
        //mesh_writer.WriteFilesUsingMesh(mesh);
        TS_ASSERT(mesh.CheckVoronoi());

        // create expected cell
        c_vector<double, 3> vertex1;
        vertex1(0)= -0.2500;
        vertex1(1)=-0.2500;
        vertex1(2)=1.2500;
        c_vector<double, 3> vertex2;
        vertex2(0)=1.2500;
        vertex2(1)=-0.2500;
        vertex2(2)=-0.2500;
        c_vector<double, 3> vertex3;
        vertex3(0)= -0.2500;
        vertex3(1)=1.2500;
        vertex3(2)=-0.2500;
        c_vector<double, 3> vertex4;
        vertex4(0)= 1.2500;
        vertex4(1)=1.2500;
        vertex4(2)=1.2500;

        Face<3> face1;
        face1.mVertices.push_back(&vertex2);
        face1.mVertices.push_back(&vertex3);
        face1.mVertices.push_back(&vertex4);
        Face<3> face2;
        face2.mVertices.push_back(&vertex1);
        face2.mVertices.push_back(&vertex4);
        face2.mVertices.push_back(&vertex3);
        Face<3> face3;
        face3.mVertices.push_back(&vertex1);
        face3.mVertices.push_back(&vertex2);
        face3.mVertices.push_back(&vertex4);
        Face<3> face4;
        face4.mVertices.push_back(&vertex1);
        face4.mVertices.push_back(&vertex3);
        face4.mVertices.push_back(&vertex2);

        VoronoiCell cell;
        cell.mFaces.push_back(&face1);
        cell.mOrientations.push_back(true);
        cell.mFaces.push_back(&face2);
        cell.mOrientations.push_back(true);
        cell.mFaces.push_back(&face3);
        cell.mOrientations.push_back(true);
        cell.mFaces.push_back(&face4);
        cell.mOrientations.push_back(true);

        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);

        // check tesellation is correct
        for (unsigned cell_index=0; cell_index<mesh.GetNumNodes(); cell_index++)
        {
            if ( mesh.GetNode(cell_index)->IsBoundaryNode() )
            {
                // for a boundary node we get a null cell
                TS_ASSERT(tessellation.rGetCell(cell_index).mFaces.size()==0);
            }
            else
            {
                // only node 4 is a non-boundary node
                TS_ASSERT_EQUALS(cell_index, 4u);
                // this gives the expected cell
                TS_ASSERT_EQUALS(tessellation.rGetCell(cell_index), cell);

                VoronoiCell cell = tessellation.rGetCell(cell_index);
                // there will only be one non-boundary cell so we know what the position is.
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[0], 0.5,1e-5);
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[1], 0.5,1e-5);
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[2], 0.5,1e-5);
            }
        }
    }


    void TestTessellation2d() throw (Exception)
    {
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true,  0,  0));
        nodes.push_back(new Node<2>(0, true,  0,  1));
        nodes.push_back(new Node<2>(0, true,  1,  1));
        nodes.push_back(new Node<2>(0, true,  1,  0));
        nodes.push_back(new Node<2>(0, true,  0.5,0.5));

        RefinableMesh<2,2> mesh(nodes);

        TS_ASSERT(mesh.CheckVoronoi());

        // Create Voronoi Tesselation
        VoronoiTessellation<2> tessellation(mesh);

        Face<2> expected_face;
        c_vector<double, 2> vertex1;
        vertex1(0)= 0.5;
        vertex1(1)= 0.0;
        expected_face.mVertices.push_back(&vertex1);
        c_vector<double, 2> vertex2;
        vertex2(0)= 1.0;
        vertex2(1)= 0.5;
        expected_face.mVertices.push_back(&vertex2);
        c_vector<double, 2> vertex3;
        vertex3(0)= 0.5;
        vertex3(1)= 1.0;
        expected_face.mVertices.push_back(&vertex3);
        c_vector<double, 2> vertex4;
        vertex4(0)= 0.0;
        vertex4(1)= 0.5;
        expected_face.mVertices.push_back(&vertex4);

        TS_ASSERT_EQUALS(*(tessellation.GetFace(4)), expected_face);
        TS_ASSERT_EQUALS( tessellation.GetFace(4)->GetNumVertices(), 4u);

        std::vector< c_vector<double, 2>*> vertices_of_face4 = tessellation.GetFace(4)->GetVertices();
        c_vector<double, 2> first_vertex_of_face4 = *(vertices_of_face4[0]);
        TS_ASSERT_DELTA( first_vertex_of_face4(0), 0.5,1e-4);
        //  Calculate length of voronoi edge between nodes 4 and 2

        TS_ASSERT_DELTA(tessellation.GetEdgeLength(4u, 2u), pow(2.0, -0.5), 1e-7);
    }

    void TestTessellation2dComplex() throw (Exception)
    {
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true,  0,  0));
        nodes.push_back(new Node<2>(0, true,  0,  1));
        nodes.push_back(new Node<2>(0, true,  -1,  0));
        nodes.push_back(new Node<2>(0, true,  1,  0));
        nodes.push_back(new Node<2>(0, true,  0.5,-pow(3,0.5)/2.0));
        nodes.push_back(new Node<2>(0, true,  -0.5,-pow(3,0.5)/2.0));

        RefinableMesh<2,2> mesh(nodes);

        TS_ASSERT(mesh.CheckVoronoi());

        // Create Voronoi Tesselation
        VoronoiTessellation<2> tessellation(mesh);

        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 1u), 1.0, 1e-6);
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 2u), 0.5 + pow(3,-0.5)/2.0, 1e-6);
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 3u), 0.5 + pow(3,-0.5)/2.0, 1e-6);
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 4u), pow(3,-0.5), 1e-6);
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 5u), pow(3,-0.5), 1e-6);

        TS_ASSERT_DELTA(tessellation.GetFace(0)->GetArea(), pow(3, 0.5)/4.0+0.5, 1e-6);
        TS_ASSERT_DELTA(tessellation.GetFace(0)->GetPerimeter(), 2.0 + pow(3, 0.5) , 1e-6);

    }


    void TestOrderVerticesAntiClockwise()
    {

        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true,  0,  0));
        nodes.push_back(new Node<2>(0, true,  0,  1));
        nodes.push_back(new Node<2>(0, true,  1,  1));
        nodes.push_back(new Node<2>(0, true,  1,  0));
        nodes.push_back(new Node<2>(0, true,  0.5,0.5));

        RefinableMesh<2,2> mesh(nodes);

        TS_ASSERT(mesh.CheckVoronoi());

        // Create Voronoi Tesselation
        VoronoiTessellation<2> tessellation(mesh);

        Face<2> original_face = *(tessellation.GetFace(0u));

        Face<2>* p_reordered_face = tessellation.mFaces[0];


        c_vector<double,2>* p_temp = p_reordered_face->mVertices[0];

        //std::cout<< "\n Before " << *(p_reordered_face->mVertices[0]);
        //std::cout<< "\n"<< *(original_face.mVertices[0]);
        //std::cout<< "\n Before " << *(p_reordered_face->mVertices[1]);
        //std::cout<< "\n"<< *(original_face.mVertices[1]);

        p_reordered_face->mVertices[0] = p_reordered_face->mVertices[1];
        p_reordered_face->mVertices[1] = p_temp;

        //std::cout<< "\n After " << *(p_reordered_face->mVertices[0]);
        //std::cout<< "\n"<< *(original_face.mVertices[0]);
        //std::cout<< "\n After " << *(p_reordered_face->mVertices[1]);
        //std::cout<< "\n"<< *(original_face.mVertices[1]);

        //tessellation.OrderFaceVerticesAntiClockwise(0u);

        TS_ASSERT_EQUALS(*p_reordered_face,original_face);

        //TS_ASSERT_DIFFERS(original_face,original_face);


    }
//    void TestWhetherMutationsSpread() throw (Exception)
//    {
//        SimulationTime* p_simulation_time = SimulationTime::Instance();
//        p_simulation_time->SetStartTime(0.0);
//
//        /*
//         * We load the steady state that the profiled test uses
//         */
//        std::string test_to_profile = "NiceCryptSim";
//        double load_time = 350;   // this is the folder and time that the stored results were archived (needed to know foldernames)
//        double time_of_each_run = 10; // run for 10 hours.
//        double end_of_simulation = 1000;
//
//
//        // Call a function to label a cell
//        TissueSimulation<2>* p_simulator = TissueSimulation<2>::Load(test_to_profile,load_time);
//        unsigned label_this = Label();
//        p_simulator->rGetTissue().rGetCellAtNodeIndex(label_this).SetMutationState(LABELLED);
//        p_simulator->Save();
//
//        // write out to file which cell it was
//        OutputFileHandler results_handler("NiceCryptSim",false);
//        out_stream file=results_handler.OpenOutputFile("overall_results.dat");
//        std::vector<double> position = p_simulator->GetNodeLocation(label_this);
//        (*file) << "Node = " << label_this << " at x = " << position[0] << "\ty = " << position[1] << "\n" << std::flush;
//
//        delete p_simulator;
//
//
//        SimulationTime::Destroy();
//        RandomNumberGenerator::Destroy();
//    }
};

#endif /*TESTVORONOITESSELLATION_HPP_*/
