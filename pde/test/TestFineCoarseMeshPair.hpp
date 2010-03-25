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
    void TestWithCoarseContainedInFine() throw(Exception)
    {
        // fine mesh is has h=0.1, on unit cube (so 6000 elements)
        TetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructCuboid(10,10,10);
        fine_mesh.Scale(0.1, 0.1, 0.1);

        // coarse mesh is has h=1 on unit cube (so 6 elements)
        QuadraticMesh<3> coarse_mesh(1.0, 1.0, 1.0, 1, 1, 1);

        FineCoarseMeshPair<3> mesh_pair(fine_mesh,coarse_mesh);

        TS_ASSERT_EQUALS(mesh_pair.mIdenticalMeshes, false);

        // check min values on fine mesh have been computed correctly
        TS_ASSERT_DELTA(mesh_pair.mMinValuesFine(0), 0.0, 1e-8);
        TS_ASSERT_DELTA(mesh_pair.mMinValuesFine(1), 0.0, 1e-8);
        TS_ASSERT_DELTA(mesh_pair.mMinValuesFine(2), 0.0, 1e-8);
        TS_ASSERT_DELTA(mesh_pair.mMaxValuesFine(0), 1.0, 1e-8);
        TS_ASSERT_DELTA(mesh_pair.mMaxValuesFine(1), 1.0, 1e-8);
        TS_ASSERT_DELTA(mesh_pair.mMaxValuesFine(2), 1.0, 1e-8);

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
        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule);

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

        mesh_pair.DeleteBoxCollection();
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

        TS_ASSERT_EQUALS(mesh_pair.mIdenticalMeshes, false);

        GaussianQuadratureRule<3> quad_rule(3);
        // need to call SetUpBoxesOnFineMesh first
        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule), "Call");

        mesh_pair.SetUpBoxesOnFineMesh(0.3);
        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 4*4*4u);

        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule);

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
    
    void TestWithIdenticalMeshes() throw(Exception)
    {
        TrianglesMeshReader<1,1> reader1("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> fine_mesh;
        fine_mesh.ConstructFromMeshReader(reader1);
        
        TrianglesMeshReader<1,1> reader2("mesh/test/data/1D_0_to_1_10_elements_quadratic",2);
        QuadraticMesh<1> coarse_mesh;
        coarse_mesh.ConstructFromMeshReader(reader2);
        
        FineCoarseMeshPair<1> mesh_pair(fine_mesh,coarse_mesh);
        TS_ASSERT_EQUALS(mesh_pair.mIdenticalMeshes, true);

        GaussianQuadratureRule<1> quad_rule(1);
        mesh_pair.SetUpBoxesOnFineMesh(0.3);

        // Covers the mIdenticalMeshes=true part of this method. Would throw exception if can't find
        // quad point in first choice of element.
        TS_ASSERT_THROWS_NOTHING(mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule));
    }

    void TestWithDefaultBoxWidth() throw(Exception)
    {
        // fine mesh is has h=0.1, on unit cube (so 6000 elements)
        TetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructCuboid(10,10,10);
        fine_mesh.Scale(0.1, 0.1, 0.1);

        // coarse mesh is has h=1 on unit cube (so 6 elements)
        QuadraticMesh<3> coarse_mesh(1.0, 1.0, 1.0, 1, 1, 1);

        FineCoarseMeshPair<3> mesh_pair(fine_mesh,coarse_mesh);

        mesh_pair.SetUpBoxesOnFineMesh();

        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 15*15*15u);
    }


////    // Uncomment this test and stick in your own two meshes to see if they will be suitable for
////    // an electromechanics simulation. There shouldn't be two many quadrature points reported to
////    // be outside the fine mesh, and for those that are, the weights (corresponding to the nearest
////    // elements) shouldn't be too far away from the interval [0,1]
////    void longrunningTestExperiment() throw(Exception)
////    {
////        TetrahedralMesh<3,3> electrics_mesh;
////        QuadraticMesh<3> mechanics_mesh;
////
////        {
////            TrianglesMeshReader<3,3> reader1("projects/pras/test/data/Rat/rat_d1");
////            electrics_mesh.ConstructFromMeshReader(reader1);
////
////            TrianglesMeshReader<3,3> reader2("/home/chaste/Desktop/rat_8125nodes_quadratic",2,2);
////            mechanics_mesh.ConstructFromMeshReader(reader2);
////        }
////
////        FineCoarseMeshPair<3> mesh_pair(electrics_mesh,mechanics_mesh);
////
////        std::cout << "Min & max x = " << mesh_pair.mMinValuesFine(0) << " " << mesh_pair.mMaxValuesFine(0) << "\n";
////
////        mesh_pair.SetUpBoxesOnFineMesh();
////
////        std::cout << "Num boxes = " << mesh_pair.mpFineMeshBoxCollection->GetNumBoxes() << "\n";
////
////        GaussianQuadratureRule<3> quad_rule(3);
////        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule);
////
////        mesh_pair.PrintStatistics();
////    }

};

#endif /*TESTFINECOARSEMESHPAIR_HPP_*/
