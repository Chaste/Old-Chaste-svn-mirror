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

#ifndef TESTPARALLELTETRAHEDRALMESH_HPP_
#define TESTPARALLELTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

class TestParallelTetrahedralMesh : public CxxTest::TestSuite
{
    
public:

    // ticket #922: 1D parallel meshes not supported (!). Since culling internal faces is mandatory here and there's not a parallel implementation of that yet
    void dontTestConstructFromMeshReader1D()
    {        
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");
        
        ParallelTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 11U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10U);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1U);
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            try
            {
                unsigned region = mesh.GetElement(i)->GetRegion();
                TS_ASSERT_EQUALS(region, i%5+1);
            }
            catch(Exception& e)
            {
                // I don't own this element do I?               
            }
        }
    }

    void TestConstructFromMeshReader2DWithoutReordering()
    {
        /*
         * In this test we don't use reordering since we want to check that a TetrahedralMesh and
         * a ParallelTetrahedralMesh create the same geometry from the same file.
         */
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ParallelTetrahedralMesh<2,2> mesh(ParallelTetrahedralMesh<2,2>::DUMB); // No reordering
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100U);
        
        // For coverage purposes
        mesh.SetElementOwnerships(0,1); // see comment in ParallelTetrahedralMesh
        
        for(ParallelTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
            iter != mesh.GetElementIteratorEnd();
            ++iter)
        {
            TS_ASSERT((*iter)->GetOwnership());
        }            
        
        // Check the inverse Jacobian
        c_matrix<double, 2, 2> jacobian;
        double jacobian_determinant;
        c_matrix<double, 2, 2> inverse_jacobian;
        
        c_matrix<double, 2, 2> element_jacobian;
        double element_jacobian_determinant;
        c_matrix<double, 2, 2> element_inverse_jacobian;
        
        mesh.GetInverseJacobianForElement(0, jacobian, jacobian_determinant, inverse_jacobian);
        
        mesh.GetElement(0)->CalculateInverseJacobian(element_jacobian, element_jacobian_determinant, element_inverse_jacobian);
        
        TS_ASSERT_EQUALS(element_jacobian_determinant, jacobian_determinant);
        
        for (unsigned row=0; row<2; row++)
        {
            for (unsigned col=0; col<2; col++)
            {
                TS_ASSERT_EQUALS(element_inverse_jacobian(row,col), inverse_jacobian(row,col));                
            }            
        }
        
        c_vector<double, 2> direction;
        c_vector<double, 2> element_direction;
        
        mesh.GetWeightedDirectionForBoundaryElement(0, direction, jacobian_determinant);
        mesh.GetBoundaryElement(0)->CalculateWeightedDirection(element_direction, element_jacobian_determinant);
        
        TS_ASSERT_EQUALS(element_jacobian_determinant, jacobian_determinant);
        
        for (unsigned row=0; row<2; row++)
        {
            TS_ASSERT_EQUALS(element_direction(row), direction(row));
        }
    
        mesh_reader.Reset();        
        TetrahedralMesh<2,2> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);
        
        for (ParallelTetrahedralMesh<2,2>::ElementIterator it=mesh.GetElementIteratorBegin(); 
             it!=mesh.GetElementIteratorEnd(); 
             ++it)
        {
            Element<2,2>* p_para_element = *it;
            unsigned element_index = p_para_element->GetIndex();
            
            Element<2,2>* p_sequ_element = seq_mesh.GetElement(element_index);            
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());
            
            for (unsigned node_local_index=0; node_local_index < p_para_element->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_EQUALS(p_para_element->GetNodeGlobalIndex(node_local_index), 
                                 p_sequ_element->GetNodeGlobalIndex(node_local_index));                                 

                TS_ASSERT_EQUALS(p_para_element->GetNode(node_local_index)->GetPoint()[0], 
                                 p_sequ_element->GetNode(node_local_index)->GetPoint()[0]);                                 
            }
        }

        for (ParallelTetrahedralMesh<2,2>::BoundaryElementIterator it=mesh.GetBoundaryElementIteratorBegin(); 
             it!=mesh.GetBoundaryElementIteratorEnd(); 
             ++it)
        {
            BoundaryElement<1,2>* p_para_boundary_element = *it;
            unsigned boundary_element_index = p_para_boundary_element->GetIndex();
            
            BoundaryElement<1,2>* p_sequ_boundary_element = seq_mesh.GetBoundaryElement(boundary_element_index);            
            TS_ASSERT_EQUALS(boundary_element_index, p_sequ_boundary_element->GetIndex());
            
            for (unsigned node_local_index=0; node_local_index < p_para_boundary_element->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_EQUALS(p_para_boundary_element->GetNodeGlobalIndex(node_local_index), 
                                 p_sequ_boundary_element->GetNodeGlobalIndex(node_local_index));                                 

                TS_ASSERT_EQUALS(p_para_boundary_element->GetNode(node_local_index)->GetPoint()[0], 
                                 p_sequ_boundary_element->GetNode(node_local_index)->GetPoint()[0]);                                 
            }
            
            
        }
         
    }

    void TestConstructFromMeshReader3D()
    {
        /*
         * In this test we let METIS reorder the ParallelTetrahedralMesh. We want to check that although
         * the indices of the nodes have changed, the location of the nodes is consistent with a 
         * TetrahedralMesh representation of the same mesh.
         */
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        
        ParallelTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 51U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 96U);
    
        mesh_reader.Reset();        
        TetrahedralMesh<3,3> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);
        
        for (ParallelTetrahedralMesh<3,3>::ElementIterator it=mesh.GetElementIteratorBegin(); 
             it!=mesh.GetElementIteratorEnd(); 
             ++it)
        {
            Element<3,3>* p_para_element = *it;
            unsigned element_index = p_para_element->GetIndex();
            
            Element<3,3>* p_sequ_element = seq_mesh.GetElement(element_index);            

            // The elements have the same index and the nodes are located in the same position.
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());            
            for (unsigned node_local_index=0; node_local_index < p_para_element->GetNumNodes(); node_local_index++)
            {
                for (unsigned dim=0; dim<3; dim++)
                {
                    TS_ASSERT_EQUALS(p_para_element->GetNode(node_local_index)->GetPoint()[dim], 
                                     p_sequ_element->GetNode(node_local_index)->GetPoint()[dim]);
                }                                 
            }
        }

        for (ParallelTetrahedralMesh<3,3>::BoundaryElementIterator it=mesh.GetBoundaryElementIteratorBegin(); 
             it!=mesh.GetBoundaryElementIteratorEnd(); 
             ++it)
        {
            BoundaryElement<2,3>* p_para_boundary_element = *it;
            unsigned boundary_element_index = p_para_boundary_element->GetIndex();
            
            BoundaryElement<2,3>* p_sequ_boundary_element = seq_mesh.GetBoundaryElement(boundary_element_index);            

            // The boundary elements have the same index and the nodes are located in the same position.
            TS_ASSERT_EQUALS(boundary_element_index, p_sequ_boundary_element->GetIndex());            
            for (unsigned node_local_index=0; node_local_index < p_para_boundary_element->GetNumNodes(); node_local_index++)
            {
                for (unsigned dim=0; dim<3; dim++)
                {
                    TS_ASSERT_EQUALS(p_para_boundary_element->GetNode(node_local_index)->GetPoint()[dim], 
                                     p_sequ_boundary_element->GetNode(node_local_index)->GetPoint()[dim]);
                }                                 
            }
            
            
        }
         
    }

    void TestEverythingIsAssignedMetisLibrary()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");

        ParallelTetrahedralMesh<3,3> mesh(ParallelTetrahedralMesh<3,3>::METIS_LIBRARY);        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());
        
        /*
         * All the nodes have been assigned
         */
        {
            unsigned num_global_nodes = mesh.GetNumNodes();
            unsigned nodes_owned[num_global_nodes];
    
            // Create a local map of the nodes this processor owns
            for (unsigned node_id=0; node_id<num_global_nodes; node_id++)
            {
                try
                {
                    unsigned node_index = mesh.GetNode(node_id)->GetIndex();
                    TS_ASSERT_EQUALS(node_id, node_index);
    
                    nodes_owned[node_index] = 1;
                }
                catch(Exception& e)
                {
                    nodes_owned[node_id] = 0;               
                }
            }
                    
            // Combine all the local maps by adding them up in the master process                
            unsigned nodes_reduction[num_global_nodes];        
            MPI_Reduce(&nodes_owned, &nodes_reduction, num_global_nodes, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);        
    
            //Make sure every node is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned node_id=0; node_id<num_global_nodes; node_id++)
                {
                    TS_ASSERT(nodes_reduction[node_id] > 0);
                }
            }
        }        
        
        /*
         * All elements have been assigned
         */
        {
            unsigned num_global_elements = mesh.GetNumElements();
            unsigned elements_owned[num_global_elements];
    
            // Create a local map of the elements this processor owns
            for (unsigned element_id=0; element_id<num_global_elements; element_id++)
            {
                try
                {
                    unsigned element_index = mesh.GetElement(element_id)->GetIndex();
                    TS_ASSERT_EQUALS(element_id, element_index);
    
                    elements_owned[element_index] = 1;
                }
                catch(Exception& e)
                {
                    elements_owned[element_id] = 0;               
                }
            }
                    
            // Combine all the local maps by adding them up in the master process                
            unsigned elements_reduction[num_global_elements];        
            MPI_Reduce(&elements_owned, &elements_reduction, num_global_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);        
    
            //Make sure every element is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned element_id=0; element_id<num_global_elements; element_id++)
                {
                    TS_ASSERT(elements_reduction[element_id] > 0);
                }
            }
        }
                    
        /*
         * All boundary elements have been assigned
         */
        {
            unsigned num_global_b_elements = mesh.GetNumBoundaryElements();
            unsigned b_elements_owned[num_global_b_elements];
    
            // Create a local map of the boundary elements this processor owns
            for (unsigned b_element_id=0; b_element_id<num_global_b_elements; b_element_id++)
            {
                try
                {
                    unsigned b_element_index = mesh.GetElement(b_element_id)->GetIndex();
                    TS_ASSERT_EQUALS(b_element_id, b_element_index);
    
                    b_elements_owned[b_element_index] = 1;
                }
                catch(Exception& e)
                {
                    b_elements_owned[b_element_id] = 0;               
                }
            }
                    
            // Combine all the local maps by adding them up in the master process                
            unsigned b_elements_reduction[num_global_b_elements];        
            MPI_Reduce(&b_elements_owned, &b_elements_reduction, num_global_b_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);        
    
            //Make sure every boundary element is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned b_element_id=0; b_element_id<num_global_b_elements; b_element_id++)
                {
                    TS_ASSERT(b_elements_reduction[b_element_id] > 0);
                }
            }
        }        
    }    
    
    void TestEverythingIsAssignedMetisBinary()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");

        ParallelTetrahedralMesh<3,3> mesh(ParallelTetrahedralMesh<3,3>::METIS_BINARY);        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());
        
        /*
         * All the nodes have been assigned
         */
        {
            unsigned num_global_nodes = mesh.GetNumNodes();
            unsigned nodes_owned[num_global_nodes];
    
            // Create a local map of the nodes this processor owns
            for (unsigned node_id=0; node_id<num_global_nodes; node_id++)
            {
                try
                {
                    unsigned node_index = mesh.GetNode(node_id)->GetIndex();
                    TS_ASSERT_EQUALS(node_id, node_index);
    
                    nodes_owned[node_index] = 1;
                }
                catch(Exception& e)
                {
                    nodes_owned[node_id] = 0;               
                }
            }
                    
            // Combine all the local maps by adding them up in the master process                
            unsigned nodes_reduction[num_global_nodes];        
            MPI_Reduce(&nodes_owned, &nodes_reduction, num_global_nodes, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);        
    
            //Make sure every node is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned node_id=0; node_id<num_global_nodes; node_id++)
                {
                    TS_ASSERT(nodes_reduction[node_id] > 0);
                }
            }
        }        
        
        /*
         * All elements have been assigned
         */
        {
            unsigned num_global_elements = mesh.GetNumElements();
            unsigned elements_owned[num_global_elements];
    
            // Create a local map of the elements this processor owns
            for (unsigned element_id=0; element_id<num_global_elements; element_id++)
            {
                try
                {
                    unsigned element_index = mesh.GetElement(element_id)->GetIndex();
                    TS_ASSERT_EQUALS(element_id, element_index);
    
                    elements_owned[element_index] = 1;
                }
                catch(Exception& e)
                {
                    elements_owned[element_id] = 0;               
                }
            }
                    
            // Combine all the local maps by adding them up in the master process                
            unsigned elements_reduction[num_global_elements];        
            MPI_Reduce(&elements_owned, &elements_reduction, num_global_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);        
    
            //Make sure every element is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned element_id=0; element_id<num_global_elements; element_id++)
                {
                    TS_ASSERT(elements_reduction[element_id] > 0);
                }
            }
        }
                    
        /*
         * All boundary elements have been assigned
         */
        {
            unsigned num_global_b_elements = mesh.GetNumBoundaryElements();
            unsigned b_elements_owned[num_global_b_elements];
    
            // Create a local map of the boundary elements this processor owns
            for (unsigned b_element_id=0; b_element_id<num_global_b_elements; b_element_id++)
            {
                try
                {
                    unsigned b_element_index = mesh.GetElement(b_element_id)->GetIndex();
                    TS_ASSERT_EQUALS(b_element_id, b_element_index);
    
                    b_elements_owned[b_element_index] = 1;
                }
                catch(Exception& e)
                {
                    b_elements_owned[b_element_id] = 0;               
                }
            }
                    
            // Combine all the local maps by adding them up in the master process                
            unsigned b_elements_reduction[num_global_b_elements];        
            MPI_Reduce(&b_elements_owned, &b_elements_reduction, num_global_b_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);        
    
            //Make sure every boundary element is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned b_element_id=0; b_element_id<num_global_b_elements; b_element_id++)
                {
                    TS_ASSERT(b_elements_reduction[b_element_id] > 0);
                }
            }
        }        
    }

    void TestConstruct3DWithRegions() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart_positive_flags");
        ParallelTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1U);

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            try
            {
                unsigned region = mesh.GetElement(i)->GetRegion();
                TS_ASSERT_EQUALS(region, (i+1)%3+1);
            }
            catch(Exception& e)
            {
                // I don't own this element do I?               
            }
        }

       TS_ASSERT_EQUALS(mesh_reader.GetNumFaceAttributes(), 1U);

        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            
            try
            {
                unsigned region = mesh.GetBoundaryElement(i)->GetRegion();
                TS_ASSERT_LESS_THAN(0u, region);
                TS_ASSERT_LESS_THAN(region, 5u);
            }
            catch(Exception& e)
            {
                // I don't own this element do I?               
            }
        }
    }
    
    void TestMetisPartitioning()
    {
        EXIT_IF_SEQUENTIAL;
        
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ParallelTetrahedralMesh<3,3> mesh(ParallelTetrahedralMesh<3,3>::METIS_LIBRARY);        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Check that each processor owns the number of nodes corresponding to its METIS partition
        std::vector<unsigned> nodes_per_processor = mesh.rGetNodesPerProcessor();
        TS_ASSERT_EQUALS(nodes_per_processor[PetscTools::GetMyRank()], mesh.GetNumLocalNodes());
    }
    
};
#endif /*TESTPARALLELTETRAHEDRALMESH_HPP_*/
