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

#ifndef TESTFINECOARSEMESHPAIR_HPP_
#define TESTFINECOARSEMESHPAIR_HPP_

#include <cxxtest/TestSuite.h>
#include "FineCoarseMeshPair.hpp"
#include "QuadraturePointsGroup.hpp"


class TestFineCoarseMeshPair : public CxxTest::TestSuite
{
public:
    // simple test where the whole of the coarse mesh is in one fine element
    void TestComputeFineElemsAndWeightsForQuadPointsSimple() throw(Exception)
    {
        TetrahedralMesh<2,2> fine_mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        fine_mesh.ConstructFromMeshReader(mesh_reader);

        QuadraticMesh<2> coarse_mesh(0.1, 0.1, 1, 1);
        coarse_mesh.Translate(0.5,0.0); // whole of the coarse mesh in now in fine element with index 1

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);

        mesh_pair.SetUpBoxesOnFineMesh();
        GaussianQuadratureRule<2> quad_rule(3);
        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);

        // all coarse quadrature points should have been found in the fine mesh
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 0u);

        // check the elements and weights have been set up correctly
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 18u);

        for(unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[i].ElementNum, 1u);
        }

        for(unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());

            // All the weights should be between 0 and 1 as no coarse nodes are 
            // Note weights = (1-psi_x-psi_y, psi_x, psi_y), where psi is the position of the
            // point in that element when transformed to the canonical element
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_LESS_THAN(-1e14, mesh_pair.rGetElementsAndWeights()[i].Weights(j));
                TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].Weights(j), 1.0+1e-14);
            }
        }

        TS_ASSERT_EQUALS(mesh_pair.mCounters[0], 18u);
        TS_ASSERT_EQUALS(mesh_pair.mCounters[1], 0u);
        TS_ASSERT_EQUALS(mesh_pair.mCounters[2], 0u);
    }

    void TestWithCoarseContainedInFine() throw(Exception)
    {
        // fine mesh is has h=0.1, on unit cube (so 6000 elements)
        TetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructCuboid(10,10,10);
        fine_mesh.Scale(0.1, 0.1, 0.1);

        // coarse mesh is has h=1 on unit cube (so 6 elements)
        QuadraticMesh<3> coarse_mesh(1.0, 1.0, 1.0, 1, 1, 1);

        FineCoarseMeshPair<3> mesh_pair(fine_mesh,coarse_mesh);

        mesh_pair.SetUpBoxesOnFineMesh(0.3);

        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 4*4*4u);

        // For each node, find containing box. That box should contain any element that node is in.
        for(unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            unsigned box_index = mesh_pair.mpFineMeshBoxCollection->CalculateContainingBox(fine_mesh.GetNode(i));

            assert(fine_mesh.GetNode(i)->rGetContainingElementIndices().size() > 0);

            for(std::set<unsigned>::iterator iter = fine_mesh.GetNode(i)->rGetContainingElementIndices().begin();
                iter != fine_mesh.GetNode(i)->rGetContainingElementIndices().end();
                ++iter)
            {
                Element<3,3>* p_element = fine_mesh.GetElement(*iter);
                TS_ASSERT_DIFFERS( mesh_pair.mpFineMeshBoxCollection->rGetBox(box_index).rGetElementsContained().find(p_element), mesh_pair.mpFineMeshBoxCollection->rGetBox(box_index).rGetElementsContained().end() )
            }
        }

        GaussianQuadratureRule<3> quad_rule(3);
        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);

        // all coarse quadrature points should have been found in the fine mesh
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 0u);

        // check the elements and weights have been set up correctly
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 6*3*3*3u);

        // some hardcoded values, just to check element_nums not all zero
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[0].ElementNum, 4846u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[10].ElementNum, 3149u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[20].ElementNum, 1209u);

        for(unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());

            // as all the quadrature points should have been found in the fine mesh, all the weights
            // should be between 0 and 1.
            // Note weights = (1-psi_x-psi_y-psi_z, psi_x, psi_y, psi_z), where psi is the position of the
            // point in that element when transformed to the canonical element
            for(unsigned j=0; j<4; j++)
            {
                TS_ASSERT_LESS_THAN(-1e14, mesh_pair.rGetElementsAndWeights()[i].Weights(j));
                TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].Weights(j), 1.0+1e-14);
            }
        }

        TS_ASSERT_EQUALS(mesh_pair.mCounters[0], 6*3*3*3u);
        TS_ASSERT_EQUALS(mesh_pair.mCounters[1], 0u);
        TS_ASSERT_EQUALS(mesh_pair.mCounters[2], 0u);
        mesh_pair.PrintStatistics();

        mesh_pair.DeleteFineBoxCollection();
        TS_ASSERT(mesh_pair.mpFineMeshBoxCollection==NULL);
    }

    void TestWithCoarseSlightlyOutsideFine() throw(Exception)
    {
        // fine mesh is has h=0.1, on unit cube (so 6000 elements)
        TetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructCuboid(10,10,10);
        fine_mesh.Scale(0.1, 0.1, 0.1);

        // coarse mesh is slightly bigger than in previous test
        QuadraticMesh<3> coarse_mesh(1.03, 1.0, 1.0, 1, 1, 1); // xmax > 1.0

        FineCoarseMeshPair<3> mesh_pair(fine_mesh,coarse_mesh);

        //TS_ASSERT_EQUALS(mesh_pair.mIdenticalMeshes, false);

        GaussianQuadratureRule<3> quad_rule(3);
        // need to call SetUpBoxesOnFineMesh first
        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true), "Call");

        mesh_pair.SetUpBoxesOnFineMesh(0.3);
        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 4*4*4u);

        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);

        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 16u); // hardcoded
        TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 16u);

        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 6*3*3*3u);

        for(unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());

            // comment out this test as now some of the weights are negative/greater than one
            //for(unsigned j=0; j<4; j++)
            //{
            //    TS_ASSERT_LESS_THAN(-1e14, mesh_pair.rGetElementsAndWeights()[i].Weights(j));
            //    TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].Weights(j), 1.0+1e-14);
            //}
        }

        // for each quadrature point that was not found in the fine mesh, check that it's x-value is greater
        // than one - this is the only way it could be outside the fine mesh
        QuadraturePointsGroup<3> quad_point_posns(coarse_mesh, quad_rule);
        for(unsigned i=0; i<mesh_pair.mNotInMesh.size(); i++)
        {
            double x = quad_point_posns.Get(mesh_pair.mNotInMesh[i])(0);
            TS_ASSERT_LESS_THAN(1.0, x);
        }


        mesh_pair.PrintStatistics();
    }

////bring back this functionality if needed    
//    void dontTestWithIdenticalMeshes() throw(Exception)
//    {
//        TrianglesMeshReader<1,1> reader1("mesh/test/data/1D_0_to_1_10_elements");
//        TetrahedralMesh<1,1> fine_mesh;
//        fine_mesh.ConstructFromMeshReader(reader1);
//        
//        TrianglesMeshReader<1,1> reader2("mesh/test/data/1D_0_to_1_10_elements_quadratic",2);
//        QuadraticMesh<1> coarse_mesh;
//        coarse_mesh.ConstructFromMeshReader(reader2);
//        
//        FineCoarseMeshPair<1> mesh_pair(fine_mesh,coarse_mesh);
//        TS_ASSERT_EQUALS(mesh_pair.mIdenticalMeshes, true);
//
//        GaussianQuadratureRule<1> quad_rule(1);
//        mesh_pair.SetUpBoxesOnFineMesh(0.3);
//
//        // Covers the mIdenticalMeshes=true part of this method. Would throw exception if can't find
//        // quad point in first choice of element.
//        TS_ASSERT_THROWS_NOTHING(mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true));
//    }

    void TestWithDefaultBoxWidth() throw(Exception)
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(10,10);
        fine_mesh.Scale(0.1, 0.1);

        QuadraticMesh<2> coarse_mesh(1.0, 1.0, 1, 1);

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);

        mesh_pair.SetUpBoxesOnFineMesh();

        // With this mesh the proposed width - the width that would correspond to 20 boxes
        // in the x-direction is:     
        //   Proposed width = 0.0552632
        // but
        //   max_edge_length = 0.141421
        // and we want width > max_edge_length, so end up with 
        //   box width = 0.155563
        // (1.1 times max edge length)
        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 64u);
        
        // now use a mesh with a smaller edge length
        TetrahedralMesh<2,2> fine_mesh2;
        fine_mesh2.ConstructRectangularMesh(100,100);
        fine_mesh2.Scale(0.01, 0.01);
        
        // Can use smaller boxes
        //  Proposed width = 0.0552632
        //  max_edge_length = 0.0141421
        //  box width = 0.0552632
        FineCoarseMeshPair<2> mesh_pair2(fine_mesh2,coarse_mesh);
        mesh_pair2.SetUpBoxesOnFineMesh();
        TS_ASSERT_EQUALS(mesh_pair2.mpFineMeshBoxCollection->GetNumBoxes(), 20*20u);
    }

    // Test when calling ComputeFineElementsAndWeightsForCoarseQuadPoints() in non-safe mode,
    // but using the default value of box width.
    // It is difficult to get the class to run incorrectly (ie fail without an assertion failing)
    // in non-safe mode (ie we can't just specify boxes that are too small), so we just test we
    // get the same results as in safe mode.
    void TestNonSafeMode() throw(Exception)
    {
        // fine mesh is has h=0.1, on unit cube (so 6000 elements)
        TetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructCuboid(10,10,10);
        fine_mesh.Scale(0.1, 0.1, 0.1);

        QuadraticMesh<3> coarse_mesh(1.03, 1.0, 1.0, 1, 1, 1); // xmax > 1.0

        FineCoarseMeshPair<3> mesh_pair(fine_mesh,coarse_mesh);
        GaussianQuadratureRule<3> quad_rule(3);

        // call SetUpBoxesOnFineMesh() without providing a width
        mesh_pair.SetUpBoxesOnFineMesh();

        // whereas before 4 by 4 by 4 boxes were explicitly chosen, here it has been determined
        // that 6 by 6 by 6 are needed
        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 6*6*6u);

        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, false /* non-safe mode*/);

        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 16u); // hardcoded
        TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 16u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 6*3*3*3u);

        for(unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());
        }

        // for each quadrature point that was not found in the fine mesh, check that it's x-value is greater
        // than one - this is the only way it could be outside the fine mesh
        QuadraturePointsGroup<3> quad_point_posns(coarse_mesh, quad_rule);
        for(unsigned i=0; i<mesh_pair.mNotInMesh.size(); i++)
        {
            double x = quad_point_posns.Get(mesh_pair.mNotInMesh[i])(0);
            TS_ASSERT_LESS_THAN(1.0, x);
        }

        mesh_pair.PrintStatistics();
    }
    
    // covers some bits that aren't covered in the tests above, 
    void TestOtherCoverage() throw(Exception)
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(10,10);
        fine_mesh.Scale(0.1, 0.1);

        QuadraticMesh<2> coarse_mesh(1.0, 1.0, 1, 1);
    
        // rotate the mesh by 45 degrees, makes it possible (since boxes no longer lined up with elements)
        // for the containing element of a quad point to be in a *local* box, ie not an element 
        // contained in the box containing this point 
        c_matrix<double,2,2> rotation_mat;
        rotation_mat(0,0) = 1.0/sqrt(2);
        rotation_mat(1,0) = -1.0/sqrt(2);
        rotation_mat(0,1) = 1.0/sqrt(2);
        rotation_mat(1,1) = 1.0/sqrt(2);
        
        fine_mesh.Rotate(rotation_mat);
        coarse_mesh.Rotate(rotation_mat);

        GaussianQuadratureRule<2> quad_rule(3);

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
        mesh_pair.SetUpBoxesOnFineMesh();
        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);

        // repeat again with smaller boxes, covers the bit requiring the whole mesh to be searched to
        // find an element for a particular quad point
        FineCoarseMeshPair<2> mesh_pair2(fine_mesh,coarse_mesh);
        mesh_pair2.SetUpBoxesOnFineMesh(0.01);
        mesh_pair2.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);
    }
    
    
    void TestComputeCoarseElementsForFineNodes() throw(Exception)
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(5,5);
        fine_mesh.Scale(0.2, 0.2);

        QuadraticMesh<2> coarse_mesh(1.0, 1.0, 1, 1); // 2 triangular elements
    
        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeCoarseElementsForFineNodes(true),"Call SetUpBoxesOnCoarseMesh()");
        
        mesh_pair.SetUpBoxesOnCoarseMesh();
        mesh_pair.ComputeCoarseElementsForFineNodes(true);
        
        for(unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            double x = fine_mesh.GetNode(i)->rGetLocation()[0];
            double y = fine_mesh.GetNode(i)->rGetLocation()[1];
            
            
            if( x+y < 1.0 - 1e-5 )  // x+y < 1
            {
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], 0u);
            }
            else if ( x+y > 1.0 + 1e-5 )  // x+y > 1
            {
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], 1u);
            }
            else // x=1-y, so in both elements, result could be either. However, it should find 0 first
            {
                //TS_ASSERT_LESS_THAN(mesh_pair.rGetCoarseElementsForFineNodes()[i], 2u);
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], 0u);
            }
        }


        // translate the fine mesh in the (-1, -1) direction --> all
        // fine nodes nearest to (not contained in) element 0. We have to 
        // make the fine mesh tiny and then translate a small amount so 
        // that it is still in the box collection for the coarse (normally the 
        // two meshes should overlap) 
        fine_mesh.Scale(1e-2, 1e-2);
        fine_mesh.Translate(-1.1e-2, -1.1e-2);
        mesh_pair.ComputeCoarseElementsForFineNodes(true);
        for(unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], 0u);
        }
        
        
        // call again with safeMode=false this time (same results, faster)
        mesh_pair.rGetCoarseElementsForFineNodes()[0] = 189342958;
        mesh_pair.ComputeCoarseElementsForFineNodes(false);
        for(unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], 0u);
        }
    }

    
    void TestComputeCoarseElementsForFineElementCentroids() throw(Exception)
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(5,5);
        fine_mesh.Scale(0.2, 0.2);

        QuadraticMesh<2> coarse_mesh(1.0, 1.0, 1, 1); // 2 triangular elements

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
        
        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeCoarseElementsForFineElementCentroids(true),"Call SetUpBoxesOnCoarseMesh()");
        
        mesh_pair.SetUpBoxesOnCoarseMesh();
        mesh_pair.ComputeCoarseElementsForFineElementCentroids(true);
        
        TS_ASSERT_EQUALS( mesh_pair.rGetCoarseElementsForFineElementCentroids().size(), fine_mesh.GetNumElements());
        for(unsigned i=0; i<fine_mesh.GetNumElements(); i++)
        {
            double x = fine_mesh.GetElement(i)->CalculateCentroid()(0);
            double y = fine_mesh.GetElement(i)->CalculateCentroid()(1);
            if(x+y < 1.0) 
            {
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineElementCentroids()[i], 0u);
            }
            else
            {
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineElementCentroids()[i], 1u);
            }
        }

        // translate the fine mesh in the (-1, -1) direction --> all
        // fine elements nearest to (not contained in) element 0. We have to 
        // make the fine mesh tiny and then translate a small amount so 
        // that it is still in the box collection for the coarse (normally the 
        // two meshes should overlap) 
        fine_mesh.Scale(1e-2, 1e-2);
        fine_mesh.Translate(-1.1e-2, -1.1e-2);
        mesh_pair.ComputeCoarseElementsForFineElementCentroids(true);
        TS_ASSERT_EQUALS( mesh_pair.rGetCoarseElementsForFineElementCentroids().size(), fine_mesh.GetNumElements());
        for(unsigned i=0; i<fine_mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS( mesh_pair.rGetCoarseElementsForFineElementCentroids()[i], 0u);
        }
    }

    void TestComputeFineElemsAndWeightsForCoarseNodes() throw(Exception)
    {
        TetrahedralMesh<2,2> fine_mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        fine_mesh.ConstructFromMeshReader(mesh_reader);

        QuadraticMesh<2> coarse_mesh(0.5, 0.5, 1, 1);
        coarse_mesh.Translate(0.2,0.1);

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);

         // need to call SetUpBoxesOnFineMesh first
        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeFineElementsAndWeightsForCoarseNodes(true), "Call");

        mesh_pair.SetUpBoxesOnFineMesh();
        mesh_pair.ComputeFineElementsAndWeightsForCoarseNodes(true);

        // all coarse quadrature points should have been found in the fine mesh
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 0u);

        // check the elements and weights have been set up correctly
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 9u);

        // check the first four nodes against what they should be
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[0].ElementNum, 1u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[1].ElementNum, 1u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[2].ElementNum, 0u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[3].ElementNum, 2u);

        for(unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());

            // All the weights should be between 0 and 1 as no coarse nodes are 
            // Note weights = (1-psi_x-psi_y-psi_z, psi_x, psi_y, psi_z), where psi is the position of the
            // point in that element when transformed to the canonical element
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_LESS_THAN(-1e14, mesh_pair.rGetElementsAndWeights()[i].Weights(j));
                TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].Weights(j), 1.0+1e-14);
            }
        }

        TS_ASSERT_EQUALS(mesh_pair.mCounters[0], 9u);
        TS_ASSERT_EQUALS(mesh_pair.mCounters[1], 0u);
        TS_ASSERT_EQUALS(mesh_pair.mCounters[2], 0u);
    }


//
//    // Uncomment this test and stick in your own two meshes to see if they will be suitable for
//    // an electromechanics simulation. There shouldn't be two many quadrature points reported to
//    // be outside the fine mesh, and for those that are, the weights (corresponding to the nearest
//    // elements) shouldn't be too far away from the interval [0,1]
//    void verylongTestExperiment() throw(Exception)
//    {
//        TetrahedralMesh<3,3> electrics_mesh;
//        QuadraticMesh<3> mechanics_mesh;
//
//        {
//            TrianglesMeshReader<3,3> reader1("/home/chaste/Desktop/heartmeshes/pras_rat_mesh/heart_d1_16_i_triangles");
//            electrics_mesh.ConstructFromMeshReader(reader1);
//
//            TrianglesMeshReader<3,3> reader2("projects/pras/test/data/Rat/rat_d1_quadratic",2,2);
//            mechanics_mesh.ConstructFromMeshReader(reader2);
//        }
//
//        FineCoarseMeshPair<3> mesh_pair(electrics_mesh,mechanics_mesh);
//
//        std::cout << "Min & max x = " << mesh_pair.mMinValuesFine(0) << " " << mesh_pair.mMaxValuesFine(0) << "\n";
//
//        mesh_pair.SetUpBoxesOnFineMesh();
//
//        std::cout << "Num boxes = " << mesh_pair.mpFineMeshBoxCollection->GetNumBoxes() << "\n";
//
//        GaussianQuadratureRule<3> quad_rule(3);
//        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, false);
//
//        mesh_pair.PrintStatistics();
//    }

};

#endif /*TESTFINECOARSEMESHPAIR_HPP_*/
