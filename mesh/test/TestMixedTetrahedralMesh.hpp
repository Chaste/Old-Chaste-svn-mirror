#ifndef TESTMIXEDTETRAHEDRALMESH_HPP_
#define TESTMIXEDTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "MixedTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "DistributedVector.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestMixedTetrahedralMesh : public CxxTest::TestSuite
{
public:

//   void TestMeshConstructionFromMeshReader(void)
// is not yet implemented
    void TestConstructionFromCuboidMeshes3D()
    {
        // create fine mesh as CTM
        
        ConformingTetrahedralMesh<3,3> fine_mesh;
        
        fine_mesh.ConstructCuboid(6, 6, 6);
        double sixth=1.0L/6.0L;
        fine_mesh.Scale(sixth, sixth, sixth);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<3,3> coarse_mesh;
        
        coarse_mesh.ConstructCuboid(3, 3, 3);
        double third=1.0L/3.0L;
        coarse_mesh.Scale(third, third, third);
        
        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        TS_ASSERT_EQUALS(coarse_mesh.GetFineMesh(), &fine_mesh);
        
        const NodeMap &node_map = coarse_mesh.rGetCoarseFineNodeMap();
        
        TS_ASSERT_EQUALS(node_map.GetNewIndex(0), 0u);
        TS_ASSERT_EQUALS(node_map.GetNewIndex(1), 2u);
        TS_ASSERT_EQUALS(node_map.GetNewIndex(4), 14u);
        TS_ASSERT_EQUALS(node_map.GetNewIndex(63), 342u);
        //Top node is 4^3-1 and 7^3-1 respectively

        TS_ASSERT_EQUALS(coarse_mesh.GetFineNodeIndexForCoarseNode(0), 0u);
        TS_ASSERT_EQUALS(coarse_mesh.GetFineNodeIndexForCoarseNode(1), 2u);
        TS_ASSERT_EQUALS(coarse_mesh.GetFineNodeIndexForCoarseNode(4), 14u);
        TS_ASSERT_EQUALS(coarse_mesh.GetFineNodeIndexForCoarseNode(63), 342u);
        
        // We're not allowed to call SetFineMesh twice
        TS_ASSERT_THROWS_ANYTHING(coarse_mesh.SetFineMesh(&fine_mesh));
    }
    
    void TestCoarseFineElementsMap2D(void)
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(2, 2, false);
        double half=1.0L/2.0L;
        fine_mesh.Scale(half, half, 0.0);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(1, 1, false);
        
        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        std::set<Element<2,2>*> expected_elements;
        expected_elements.insert(fine_mesh.GetElement(6));
        expected_elements.insert(fine_mesh.GetElement(0));
        expected_elements.insert(fine_mesh.GetElement(3));
        expected_elements.insert(fine_mesh.GetElement(2));
        // Elements 0,2,3 and 6 form the upper right half of the space
        // (added in a funny order since set equality ought to cope with this)
        
        TS_ASSERT(expected_elements ==
                  coarse_mesh.GetFineElementsForCoarseElementIndex(0));
    }
    
    void TestTransferFlags()
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(4, 4, false);
        double half=1.0L/2.0L;
        fine_mesh.Scale(half, half, 0.0);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(2, 2, false);
        
        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        bool any_elements_flagged = coarse_mesh.TransferFlags();
        TS_ASSERT_EQUALS(any_elements_flagged, false);
        
        // Flag the right half of the coarse mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            ChastePoint<2> centroid = ChastePoint<2>(element.CalculateCentroid());
            if (centroid[0] > 1.0)
            {
                element.Flag();
            }
            else
            {
                element.Unflag();
            }
        }
        
        // Flag the corresponding region of the fine mesh
        any_elements_flagged = coarse_mesh.TransferFlags();
        TS_ASSERT_EQUALS(any_elements_flagged, true);
        
        // Check
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_fine_element;
        for (i_fine_element = fine_mesh.GetElementIteratorBegin();
             i_fine_element != fine_mesh.GetElementIteratorEnd();
             i_fine_element++)
        {
            Element<2,2> &element = **i_fine_element;
            ChastePoint<2> centroid = ChastePoint<2>(element.CalculateCentroid());
            if (centroid[0] > 1.0)
            {
                TS_ASSERT(element.IsFlagged());
            }
            else
            {
                TS_ASSERT(!element.IsFlagged());
            }
        }
    }
    
    void TestFineNodesCoarseElementsMap2D(void)
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(2, 2, false);
        double half=1.0L/2.0L;
        fine_mesh.Scale(half, half, 0.0);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(1, 1, false);
        
        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        //node 1 is on the top edge of the fine mesh
        TS_ASSERT(coarse_mesh.GetElement(0) ==
                  coarse_mesh.GetACoarseElementForFineNodeIndex(1));
        //node 3 is on the left edge of the fine mesh
        TS_ASSERT(coarse_mesh.GetElement(1) ==
                  coarse_mesh.GetACoarseElementForFineNodeIndex(3));
    }
    
    void TestFineMeshIncorrect3D(void)
    {
        ConformingTetrahedralMesh<3,3> fine_mesh;
        
        fine_mesh.ConstructCuboid(5, 5, 5);
        double fifth=1.0L/5.0L;
        fine_mesh.Scale(fifth, fifth, fifth);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<3,3> coarse_mesh;
        
        coarse_mesh.ConstructCuboid(3, 3, 3);
        double third=1.0L/3.0L;
        coarse_mesh.Scale(third, third, third);
        
        // give fine mesh to coarse mesh and calculate node map
        // should throw because not every coarse node has a coincident fine node
        TS_ASSERT_THROWS_ANYTHING(coarse_mesh.SetFineMesh(&fine_mesh));
    }
    
    void TestFineAndCoarseDisc(void)
    {
        TrianglesMeshReader<2,2> fine_mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructFromMeshReader(fine_mesh_reader);
        
        TrianglesMeshReader<2,2> coarse_mesh_reader("mesh/test/data/DecimatedDisk");
        MixedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructFromMeshReader(coarse_mesh_reader);
        
        TS_ASSERT_THROWS_ANYTHING(coarse_mesh.TransferFlags());
        
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        // Flag the right semicircle of the coarse mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            ChastePoint<2> centroid = ChastePoint<2>(element.CalculateCentroid());
            if (centroid[0] > 0)
            {
                element.Flag();
            }
            else
            {
                element.Unflag();
            }
        }
        
        // Flag the corresponding region of the fine mesh
        coarse_mesh.TransferFlags();
        
        // Check - flagged fine elements must have at least one node in a flagged coarse element
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_fine_element;
        for (i_fine_element = fine_mesh.GetElementIteratorBegin();
             i_fine_element != fine_mesh.GetElementIteratorEnd();
             i_fine_element++)
        {
            Element<2,2> &fine_element = **i_fine_element;
            
            unsigned count = 0;
            for (unsigned i=0; i<fine_element.GetNumNodes(); i++)
            {
                unsigned node_index = fine_element.GetNodeGlobalIndex(i);
                const Element<2,2> *p_coarse_element = coarse_mesh.GetACoarseElementForFineNodeIndex(node_index);
                if (p_coarse_element->IsFlagged())
                {
                    count++;
                }
            }
            if (count == 0)
            {
                TS_ASSERT(!fine_element.IsFlagged());
            }
            else if (count > 1)
            {
                if (fine_element.GetIndex() != 421)
                {
                    TS_ASSERT(fine_element.IsFlagged());
                }
            }
            // Some fine elements just have 1 node touching the edge of the flagged
            // region in the coarse mesh, so are not flagged.
            // Element 421 is a special case - one edge lies on the border of the
            // flagged region.
        }
    }

    void TestInterpolateOnUnflaggedRegion()
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        
        unsigned num_elem = 16;
        fine_mesh.ConstructRectangularMesh(num_elem,num_elem);
        fine_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<2,2> coarse_mesh;
        
        num_elem = 4;
        coarse_mesh.ConstructRectangularMesh(num_elem, num_elem);
        coarse_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // give fine mesh to coarse mesh
        coarse_mesh.SetFineMesh(&fine_mesh);

        // Set linear initial condition on coarse mesh 
        DistributedVector::SetProblemSize(coarse_mesh.GetNumNodes());
        Vec coarse_ic_petsc = DistributedVector::CreateVec();
        DistributedVector coarse_ic(coarse_ic_petsc);
        
        c_vector<double,2> linear_combination;
        linear_combination(0)=3;
        linear_combination(1)=-23;
        
        for (DistributedVector::Iterator index=DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            coarse_ic[index] = inner_prod(linear_combination, 
                                          coarse_mesh.GetNode(index.Global)->rGetLocation());
        }
        coarse_ic.Restore();
        
        // Set zero initial condition on the fine mesh, so we know if interpolation changed things
        DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
        Vec fine_ic_petsc = DistributedVector::CreateVec();
        #if (PETSC_VERSION_MINOR == 2) //Old API
        PetscScalar zero = 0;
        VecSet(&zero, fine_ic_petsc);
        #else
        VecSet(fine_ic_petsc, 0);
        #endif  
        VecAssemblyBegin(fine_ic_petsc);
        VecAssemblyEnd(fine_ic_petsc);
        
        
        // Flag the right half of the coarse mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            ChastePoint<2> centroid = ChastePoint<2>(element.CalculateCentroid());
            if (centroid[0] > 0.5)
            {
                element.Flag();
            }
            else
            {
                element.Unflag();
            }
        }
        
        // Flag the corresponding region of the fine mesh
        coarse_mesh.TransferFlags();
        
        // Interpolate coarse initial condition onto UNflagged region of fine mesh
        coarse_mesh.InterpolateOnUnflaggedRegion(coarse_ic_petsc, fine_ic_petsc);

        // Check that unflagged nodes were interpolated, but flagged ones were not
        DistributedVector fine_ic(fine_ic_petsc);
        bool found_unflagged_elt = false;
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_fine_element;
        for (i_fine_element = fine_mesh.GetElementIteratorBegin();
             i_fine_element != fine_mesh.GetElementIteratorEnd();
             i_fine_element++)
        {
            Element<2,2> &element = **i_fine_element;
            for (unsigned element_node_index=0; element_node_index<3; element_node_index++)
            {
                const c_vector<double, 2>& r_node_loc = element.GetNode(element_node_index)->rGetLocation();

                double expected_value;
                if (r_node_loc[0] < 0.5)
                {
                    // Node lies only on unflagged elements
                    expected_value = inner_prod(linear_combination, r_node_loc);
                    found_unflagged_elt = true;
                }
                else
                {
                    // Node is in a flagged element
                    expected_value = 0;
                }

                unsigned element_node_global_index = element.GetNodeGlobalIndex(element_node_index);
                if (DistributedVector::IsGlobalIndexLocal(element_node_global_index))
                {
                    TS_ASSERT_DELTA(fine_ic[element_node_global_index], expected_value, 1e-12);
                }
            }
        }
        // Check the test really did something useful!
        TS_ASSERT(found_unflagged_elt);
        
        VecDestroy(coarse_ic_petsc);
        VecDestroy(fine_ic_petsc);
    }


    void TestUpdateCoarseSolutionOnflaggedRegion()
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        
        unsigned num_elem = 16;
        fine_mesh.ConstructRectangularMesh(num_elem,num_elem);
        fine_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<2,2> coarse_mesh;
        
        num_elem = 4;
        coarse_mesh.ConstructRectangularMesh(num_elem, num_elem);
        coarse_mesh.Scale(1.0/num_elem, 1.0/num_elem);
        
        // give fine mesh to coarse mesh
        coarse_mesh.SetFineMesh(&fine_mesh);

        // Set linear 'solution' vector on fine mesh 
        DistributedVector::SetProblemSize(fine_mesh.GetNumNodes());
        Vec fine_soln_petsc = DistributedVector::CreateVec();
        DistributedVector fine_soln(fine_soln_petsc);
        
        c_vector<double,2> linear_combination;
        linear_combination(0)=3;
        linear_combination(1)=-23;
        
        for (DistributedVector::Iterator index=DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            fine_soln[index] = inner_prod(linear_combination, 
                                          fine_mesh.GetNode(index.Global)->rGetLocation());
        }
        fine_soln.Restore();
        
        // Set zero solution on the coarse mesh, so we know if updating changed things
        DistributedVector::SetProblemSize(coarse_mesh.GetNumNodes());
        Vec coarse_soln_petsc = DistributedVector::CreateVec();
#if (PETSC_VERSION_MINOR == 2) //Old API
        PetscScalar zero = 0;
        VecSet(&zero, coarse_soln_petsc);
#else
        VecSet(coarse_soln_petsc, 0);
#endif  
        VecAssemblyBegin(coarse_soln_petsc);
        VecAssemblyEnd(coarse_soln_petsc);
        
        
        // Flag the right half of the coarse mesh
        ConformingTetrahedralMesh<2, 2>::ElementIterator i_coarse_element;
        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            ChastePoint<2> centroid = ChastePoint<2>(element.CalculateCentroid());
            if (centroid[0] > 0.5)
            {
                element.Flag();
            }
            else
            {
                element.Unflag();
            }
        }        

        coarse_mesh.UpdateCoarseSolutionOnFlaggedRegion(coarse_soln_petsc,fine_soln_petsc);
        
        // Check that flagged nodes were copied but unflagged ones were not
        DistributedVector coarse_soln(coarse_soln_petsc);
        bool found_flagged_elt = false;

        for (i_coarse_element = coarse_mesh.GetElementIteratorBegin();
             i_coarse_element != coarse_mesh.GetElementIteratorEnd();
             i_coarse_element++)
        {
            Element<2,2> &element = **i_coarse_element;
            for (unsigned element_node_index=0; element_node_index<3; element_node_index++)
            {
                const c_vector<double, 2>& r_node_loc = element.GetNode(element_node_index)->rGetLocation();

                double expected_value;
                if (r_node_loc[0] < 0.5)
                {
                    // Node lies only on unflagged elements
                    expected_value = 0;
                }
                else
                {
                    // Node is in a flagged element
                    expected_value = inner_prod(linear_combination, r_node_loc);
                    found_flagged_elt = true;
                }

                unsigned element_node_global_index = element.GetNodeGlobalIndex(element_node_index);
                if (DistributedVector::IsGlobalIndexLocal(element_node_global_index))
                {
                    TS_ASSERT_DELTA(coarse_soln[element_node_global_index], expected_value, 1e-12);
                }
            }
        }
        
        // Check the test really did something useful!
        TS_ASSERT(found_flagged_elt);
        
        VecDestroy(coarse_soln_petsc);
        VecDestroy(fine_soln_petsc);
    }
};





#endif /*TESTMIXEDTETRAHEDRALMESH_HPP_*/
