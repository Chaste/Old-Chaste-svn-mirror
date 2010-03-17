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

#ifndef TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_
#define TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PetscTools.hpp"
#include "ArchiveOpener.hpp"

class TestDistributedTetrahedralMesh : public CxxTest::TestSuite
{
private:

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    void CheckEverythingIsAssigned(DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
    {
        /*
         * All the nodes have been assigned
         */
        unsigned total_nodes_this_process = 0;
        {
            unsigned num_global_nodes = rMesh.GetNumNodes();
            unsigned nodes_owned[num_global_nodes];

            for (unsigned node_id=0; node_id<num_global_nodes;  node_id++)
            {
                
                try
                {
                     unsigned node_index = rMesh.GetNode(node_id)->GetIndex();
                     TS_ASSERT_EQUALS(node_id, node_index);
                     nodes_owned[node_index] = 1;
                     total_nodes_this_process++;
                }
                catch (Exception &e)
                {
                    nodes_owned[node_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalNodes(), total_nodes_this_process);
            
            // Combine all the local maps by adding them up in the master process
            unsigned nodes_reduction[num_global_nodes];
            MPI_Reduce(&nodes_owned, &nodes_reduction, num_global_nodes, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

            // Make sure every node is owned at least by one processor
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
        unsigned total_elements_this_process = 0;
        {
            unsigned num_global_elements = rMesh.GetNumElements();
            unsigned elements_owned[num_global_elements];

            // Create a local map of the elements this processor owns
            for (unsigned element_id=0; element_id<num_global_elements; element_id++)
            {
                try
                {
                    unsigned element_index = rMesh.GetElement(element_id)->GetIndex();
                    TS_ASSERT_EQUALS(element_id, element_index);

                    elements_owned[element_index] = 1;

                    total_elements_this_process++;
                }
                catch(Exception& e)
                {
                    elements_owned[element_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalElements(), total_elements_this_process);
            
            // Combine all the local maps by adding them up in the master process
            unsigned elements_reduction[num_global_elements];
            MPI_Reduce(&elements_owned, &elements_reduction, num_global_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

            // Make sure every element is owned at least by one processor
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
        unsigned total_b_elements_this_process = 0;
        {
            unsigned num_global_b_elements = rMesh.GetNumBoundaryElements();
            unsigned b_elements_owned[num_global_b_elements];

            // Create a local map of the boundary elements this processor owns
            for (unsigned b_element_id=0; b_element_id<num_global_b_elements; b_element_id++)
            {
                try
                {
                    unsigned b_element_index = rMesh.GetBoundaryElement(b_element_id)->GetIndex();
                    TS_ASSERT_EQUALS(b_element_id, b_element_index);

                    b_elements_owned[b_element_index] = 1;
                    
                    total_b_elements_this_process++;
                }
                catch(Exception& e)
                {
                    b_elements_owned[b_element_id] = 0;
                }
            }

            TS_ASSERT_EQUALS(rMesh.GetNumLocalBoundaryElements(), total_b_elements_this_process);
            
            // Combine all the local maps by adding them up in the master process
            unsigned b_elements_reduction[num_global_b_elements];
            MPI_Reduce(&b_elements_owned, &b_elements_reduction, num_global_b_elements, MPI_UNSIGNED, MPI_SUM, PetscTools::MASTER_RANK, PETSC_COMM_WORLD);

            // Make sure every boundary element is owned at least by one processor
            if (PetscTools::AmMaster())
            {
                for (unsigned b_element_id=0; b_element_id<num_global_b_elements; b_element_id++)
                {
                    TS_ASSERT(b_elements_reduction[b_element_id] > 0);
                }
            }
        }
        if (total_nodes_this_process != 0)
        {
            TS_ASSERT( 0u != total_nodes_this_process );
            TS_ASSERT( 0u != total_elements_this_process );
            TS_ASSERT( 0u != total_b_elements_this_process );
        }
        else
        {
            //Metis may allocate no nodes to a partition if the mesh is small and there are many processes
            // Look out for "You just increased the maxndoms"
            TS_ASSERT( 0u == total_nodes_this_process );
            TS_ASSERT( 0u == total_elements_this_process );
            TS_ASSERT( 0u == total_b_elements_this_process );
        }
            


    }

public:

    void TestConstructFromMeshReader1D()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");

        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 11u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 11u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);
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
         * a DistributedTetrahedralMesh create the same geometry from the same file.
         */
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMesh<2,2>::DUMB); // No reordering
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543u);
        TS_ASSERT_EQUALS(mesh.GetDistributedVectorFactory()->GetProblemSize(), 543u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100u);

        // For coverage purposes
        mesh.SetElementOwnerships(0,1); // see comment in DistributedTetrahedralMesh

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT(iter->GetOwnership());
        }

        // Check the inverse Jacobian
        c_matrix<double, 2, 2> jacobian;
        double jacobian_determinant;
        c_matrix<double, 2, 2> inverse_jacobian;

        c_matrix<double, 2, 2> element_jacobian;
        double element_jacobian_determinant;
        c_matrix<double, 2, 2> element_inverse_jacobian;

        try
        {
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
        }
        catch(Exception& e)
        {
            // I don't own this element do I?
        }

        c_vector<double, 2> direction;
        c_vector<double, 2> element_direction;

        try
        {
            mesh.GetWeightedDirectionForBoundaryElement(0, direction, jacobian_determinant);
            mesh.GetBoundaryElement(0)->CalculateWeightedDirection(element_direction, element_jacobian_determinant);

            TS_ASSERT_EQUALS(element_jacobian_determinant, jacobian_determinant);

            for (unsigned row=0; row<2; row++)
            {
                TS_ASSERT_EQUALS(element_direction(row), direction(row));
            }
        }
        catch(Exception& e)
        {
            // I don't own this boundary element do I?
        }

        TetrahedralMesh<2,2> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);

        for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();

            Element<2,2>* p_sequ_element = seq_mesh.GetElement(element_index);
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());

            for (unsigned node_local_index=0; node_local_index < iter->GetNumNodes(); node_local_index++)
            {
                TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(node_local_index),
                                 p_sequ_element->GetNodeGlobalIndex(node_local_index));

                TS_ASSERT_EQUALS(iter->GetNode(node_local_index)->GetPoint()[0],
                                 p_sequ_element->GetNode(node_local_index)->GetPoint()[0]);
            }
        }

        for (DistributedTetrahedralMesh<2,2>::BoundaryElementIterator it=mesh.GetBoundaryElementIteratorBegin();
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

    // See #1199
    void TestConstructFromMeshReader2DWithUnevenDistribution()
    {
        unsigned local_nodes = 1u;
        unsigned local_nodes_wrong = 1u;
        unsigned total_nodes = 543u;
        unsigned total_nodes_wrong = 100u;
        if (PetscTools::AmTopMost())
        {
            local_nodes = total_nodes - (PetscTools::GetNumProcs()-1);
            local_nodes_wrong = total_nodes_wrong - (PetscTools::GetNumProcs()-1);
        }
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMesh<2,2>::DUMB); // No reordering

        // Exceptions
        DistributedVectorFactory* p_wrong_factory1 = new DistributedVectorFactory(PetscTools::GetMyRank(), PetscTools::GetMyRank()+1,
                                                                                  PetscTools::GetNumProcs(), PetscTools::GetNumProcs()+1);
        TS_ASSERT_THROWS_THIS(mesh.SetDistributedVectorFactory(p_wrong_factory1),
                              "The distributed vector factory provided to the mesh is for the wrong number of processes.");
        delete p_wrong_factory1;
        
        DistributedVectorFactory* p_wrong_factory2 = new DistributedVectorFactory(total_nodes_wrong, local_nodes_wrong);
        mesh.SetDistributedVectorFactory(p_wrong_factory2);
        TS_ASSERT_THROWS_THIS(mesh.ConstructFromMeshReader(mesh_reader),
                              "The distributed vector factory size in the mesh doesn't match the total number of nodes.");
        delete p_wrong_factory2;

        // OK call
        DistributedVectorFactory* p_uneven_factory = new DistributedVectorFactory(total_nodes, local_nodes);
        mesh.SetDistributedVectorFactory(p_uneven_factory);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check the mesh is using the supplied factory
        TS_ASSERT(mesh.GetDistributedVectorFactory() == p_uneven_factory);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), total_nodes);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), total_nodes);
        TS_ASSERT_EQUALS(mesh.GetDistributedVectorFactory()->GetProblemSize(), total_nodes);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100u);
        TS_ASSERT_EQUALS(mesh.GetNumLocalNodes(), local_nodes);
        
        // Another exception
        TS_ASSERT_THROWS_THIS(mesh.SetDistributedVectorFactory(p_uneven_factory),
                              "Cannot change the mesh's distributed vector factory once it has been set.");
    }
    
    void TestConstructFromMeshReader3D()
    {
        /*
         * In this test we let METIS reorder the DistributedTetrahedralMesh. We want to check that although
         * the indices of the nodes have changed, the location of the nodes is consistent with a
         * TetrahedralMesh representation of the same mesh.
         */
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 51u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 96u);

        TetrahedralMesh<3,3> seq_mesh;
        seq_mesh.ConstructFromMeshReader(mesh_reader);

        for (AbstractTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();

            Element<3,3>* p_sequ_element = seq_mesh.GetElement(element_index);

            // The elements have the same index and the nodes are located in the same position.
            TS_ASSERT_EQUALS(element_index, p_sequ_element->GetIndex());
            for (unsigned node_local_index=0; node_local_index < iter->GetNumNodes(); node_local_index++)
            {
                for (unsigned dim=0; dim<3; dim++)
                {
                    TS_ASSERT_EQUALS(iter->GetNode(node_local_index)->GetPoint()[dim],
                                     p_sequ_element->GetNode(node_local_index)->GetPoint()[dim]);
                }
            }
        }

        for (DistributedTetrahedralMesh<3,3>::BoundaryElementIterator it=mesh.GetBoundaryElementIteratorBegin();
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

        //Scale it (for coverage)
        mesh.Scale(2.0);

    }

    void TestConstructFromMeshReaderWithBinaryFiles()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements_binary");
        TrianglesMeshReader<3,3> mesh_reader_ascii("mesh/test/data/cube_136_elements");

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<3,3> mesh_from_ascii;
        mesh_from_ascii.ConstructFromMeshReader(mesh_reader_ascii);

        // Check that we have the right number of nodes and elements
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_from_ascii.GetNumBoundaryElements());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_from_ascii.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_from_ascii.GetNumNodes());

        // Check that the nodes and elements of each mesh are identical
        for (AbstractTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();

            Element<3,3>* p_ascii_element = mesh_from_ascii.GetElement(element_index);

            // The elements have the same index and the nodes are located in the same position.
            TS_ASSERT_EQUALS(element_index, p_ascii_element->GetIndex());
            for (unsigned node_local_index=0; node_local_index < iter->GetNumNodes(); node_local_index++)
            {
//                for (unsigned dim=0; dim<3; dim++)
//                {
//                    TS_ASSERT_EQUALS(iter->GetNode(node_local_index)->GetPoint()[dim],
//                                     p_ascii_element->GetNode(node_local_index)->GetPoint()[dim]);
//                }
                TS_ASSERT_DELTA( norm_2( iter->GetNode(node_local_index)->rGetLocation() - 
                                     p_ascii_element->GetNode(node_local_index)->rGetLocation() ), 0.0, 1e-10 );
            }
        }

    }

    void TestEverythingIsAssignedMetisLibrary()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMesh<3,3>::METIS_LIBRARY);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh_reader.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh_reader.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh_reader.GetNumFaces());

        CheckEverythingIsAssigned<3,3>(mesh);
    }


    void TestConstruct3DWithRegions() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart_nonnegative_flags");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

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

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaceAttributes(), 1u);

        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            try
            {
                unsigned region = mesh.GetBoundaryElement(i)->GetRegion();
                TS_ASSERT_LESS_THAN(region, 4u);
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

        {
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
            DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMesh<3,3>::METIS_LIBRARY);
            mesh.ConstructFromMeshReader(mesh_reader);

            // Check that each processor owns the number of nodes corresponding to its METIS partition
            unsigned local_nodes = mesh.GetDistributedVectorFactory()->GetLocalOwnership();
            TS_ASSERT_EQUALS(local_nodes, mesh.GetNumLocalNodes());

            typedef DistributedTetrahedralMesh<3,3> MESH_TYPE; // To stop TS_ASSERT mistaking the comma for an argument
            TS_ASSERT_EQUALS(mesh.GetPartitionType(),MESH_TYPE::METIS_LIBRARY);
        }

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
            DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMesh<2,2>::METIS_LIBRARY);
            mesh.ConstructFromMeshReader(mesh_reader);

            // Check that each processor owns the number of nodes corresponding to its METIS partition
            unsigned local_nodes = mesh.GetDistributedVectorFactory()->GetLocalOwnership();
            TS_ASSERT_EQUALS(local_nodes, mesh.GetNumLocalNodes());
        }
    }
    
    void TestPartitioningOfEmbeddedDimensionMesh()
    {
        //Shouldn't ever use a partition other than DUMB because it's 1-D 
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/branched_1d_in_3d_mesh");
        DistributedTetrahedralMesh<1,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 31u);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 30u);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 3u);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 3u);
    }


    void TestArchiving() throw(Exception)
    {
        std::string archive_dir = "archive";
        std::string archive_file = "distributed_tetrahedral_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("distributed_tetrahedral_mesh");

        DistributedTetrahedralMesh<2,2>* p_mesh = new DistributedTetrahedralMesh<2,2>(DistributedTetrahedralMesh<2,2>::DUMB);
        //std::vector<unsigned> halo_node_indices;
        std::vector<Node<2>*> halo_nodes;
        unsigned num_nodes;
        unsigned local_num_nodes;
        unsigned num_elements;
        // archive
        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

            p_mesh->ConstructFromMeshReader(mesh_reader);
            num_nodes = p_mesh->GetNumNodes();
            local_num_nodes = p_mesh->GetNumLocalNodes();
            num_elements = p_mesh->GetNumElements();

            halo_nodes = p_mesh->mHaloNodes;

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractTetrahedralMesh<2,2>* const p_mesh_abstract = static_cast<AbstractTetrahedralMesh<2,2>* >(p_mesh);
            (*p_arch) << p_mesh_abstract;
        }

        // restore
        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<2,2>* p_mesh_abstract2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_mesh_abstract2;
            // Check we have the right number of nodes & elements
            DistributedTetrahedralMesh<2,2>* p_mesh2 = static_cast<DistributedTetrahedralMesh<2,2>*>(p_mesh_abstract2);

            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(), num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumLocalNodes(), local_num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), num_elements);

            // Check some node co-ordinates
            try
            {
                Node<2>* p_node1 = p_mesh->GetNode(0);
                Node<2>* p_node2 = p_mesh2->GetNode(0);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], p_node2->GetPoint()[1], 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Node<2>* p_node1 = p_mesh->GetNode(500);
                Node<2>* p_node2 = p_mesh2->GetNode(500);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check first element has the right nodes
            try
            {
                Element<2,2>* p_element = p_mesh->GetElement(0);
                Element<2,2>* p_element2 = p_mesh2->GetElement(0);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Element<2,2>* p_element = p_mesh->GetElement(500);
                Element<2,2>* p_element2 = p_mesh2->GetElement(500);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check the halo nodes are right
            std::vector<Node<2>*> halo_nodes2 = p_mesh2->mHaloNodes;
            TS_ASSERT_EQUALS(halo_nodes2.size(), halo_nodes.size());
            delete p_mesh2;
        }

        // restore from a single processor archive
        {
            if ( PetscTools::IsSequential() )
            {
                ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(
                        "mesh/test/data/distributed_mesh_archive/", "distributed_tetrahedral_mesh.arch", false);
                boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
                AbstractTetrahedralMesh<2,2>* p_mesh3 = NULL;
                (*p_arch) >> p_mesh3;
                delete p_mesh3;
            }
            else
            {
                typedef ArchiveOpener<boost::archive::text_iarchive, std::ifstream> InputArchiveOpener;
                if (PetscTools::GetMyRank()>0)
                {
                    /// Should not read this archive because none exists here.
                    TS_ASSERT_THROWS_CONTAINS(InputArchiveOpener arch_opener("mesh/test/data/distributed_mesh_archive/", "distributed_tetrahedral_mesh.arch", false),
                                "Cannot load secondary archive file:");
                }
                else
                {
                    /// Should not read this archive because there are two or more processes and
                    // this archive was written on one process.
                    InputArchiveOpener arch_opener("mesh/test/data/distributed_mesh_archive/", "distributed_tetrahedral_mesh.arch", false);
                    boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
                    AbstractTetrahedralMesh<2,2>* p_mesh3 = NULL;
                    TS_ASSERT_THROWS_THIS((*p_arch) >> p_mesh3,
                                          "This archive was written for a different number of processors");

                }
            }
        }

        delete p_mesh;
    }

private:

    template <unsigned DIM>
    void CompareParallelMeshOwnership(DistributedTetrahedralMesh<DIM,DIM> &readMesh, DistributedTetrahedralMesh<DIM,DIM> &constructedMesh)
    {
        typedef DistributedTetrahedralMesh<DIM,DIM> MESH_TYPE; // To stop TS_ASSERT mistaking the comma for an argument
        //The read mesh has a dumb partition in the test
        TS_ASSERT_EQUALS(readMesh.GetPartitionType(), MESH_TYPE::DUMB);
        //All constructed meshes have dumb partitioning -- so that they are invariant under archiving
        TS_ASSERT_EQUALS(constructedMesh.GetPartitionType(), MESH_TYPE::DUMB);
        TS_ASSERT_EQUALS(constructedMesh.GetDistributedVectorFactory()->GetLocalOwnership(),
                         readMesh.GetDistributedVectorFactory()->GetLocalOwnership());
        TS_ASSERT_EQUALS(constructedMesh.GetNumNodes(), readMesh.GetNumNodes());
        TS_ASSERT_EQUALS(constructedMesh.GetNumLocalNodes(), readMesh.GetNumLocalNodes());
        TS_ASSERT_EQUALS(constructedMesh.GetNumBoundaryNodes(), readMesh.GetNumBoundaryNodes());
        TS_ASSERT_EQUALS(constructedMesh.GetNumBoundaryElements(),  readMesh.GetNumBoundaryElements());
        TS_ASSERT_EQUALS(constructedMesh.GetNumElements(), readMesh.GetNumElements());
        TS_ASSERT_EQUALS(constructedMesh.GetNumLocalElements(), readMesh.GetNumLocalElements());

        for (unsigned i=0; i<readMesh.GetNumNodes(); i++)
        {
            try
            {
                unsigned index=constructedMesh.SolveNodeMapping(i);
                //Read mesh didn't throw so owns the node
                TS_ASSERT_THROWS_NOTHING(constructedMesh.GetNode(i));
                TS_ASSERT_EQUALS(index, readMesh.SolveNodeMapping(i));
             }
            catch(Exception& e)
            {
                //Read mesh threw so does not own node
                TS_ASSERT_THROWS_CONTAINS(constructedMesh.GetNode(i), "does not belong to processor");
            }
        }

        for (unsigned i=0; i<readMesh.GetNumElements(); i++)
        {
            try
            {
                unsigned index=constructedMesh.SolveElementMapping(i);
                //Read mesh didn't throw so owns the element
                TS_ASSERT_THROWS_NOTHING(constructedMesh.GetElement(i));
                TS_ASSERT_EQUALS(index, readMesh.SolveElementMapping(i));
             }
            catch(Exception& e)
            {
                //Read mesh threw so does not own element
                TS_ASSERT_THROWS_CONTAINS(constructedMesh.GetElement(i), "does not belong to processor");
            }
        }

        for (unsigned i=0; i<readMesh.GetNumBoundaryElements(); i++)
        {
            try
            {

                unsigned index=constructedMesh.SolveBoundaryElementMapping(i);
                //Read mesh didn't throw so owns the element
                TS_ASSERT_THROWS_NOTHING(constructedMesh.GetBoundaryElement(i));
                TS_ASSERT_EQUALS(index, readMesh.SolveBoundaryElementMapping(i));
             }
            catch(Exception& e)
            {
                //Read mesh threw so does not own element
                TS_ASSERT_THROWS_CONTAINS(constructedMesh.GetBoundaryElement(i), "does not belong to processor");
            }
        }

    }
public:
    void TestConstructLinearMesh()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");
        DistributedTetrahedralMesh<1,1> read_mesh;
        read_mesh.ConstructFromMeshReader(mesh_reader);
        DistributedTetrahedralMesh<1,1> constructed_mesh;
        constructed_mesh.ConstructLinearMesh(10u);

        CompareParallelMeshOwnership(read_mesh, constructed_mesh);

        unsigned owned=constructed_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        unsigned owned_in_read=read_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        TS_ASSERT_EQUALS(owned_in_read, owned);
        TS_ASSERT_EQUALS(constructed_mesh.GetNumBoundaryElements(), 2u);
        TS_ASSERT_EQUALS(constructed_mesh.GetNumLocalNodes(), owned);
        //Sequential: Process owns one fewer element than the number of nodes
        //Parallel: End processes own the same as the number of node (since one node is paired with a halo node)
        //Parallel: Middle processes own one more than the number of nodes (since two nodes are paired with a halo nodes)
        unsigned expected_elements=owned+1;
        if (PetscTools::AmMaster())
        {
            expected_elements--;
        }
         if (PetscTools::AmTopMost())
        {
            expected_elements--;
        }
        TS_ASSERT_EQUALS(constructed_mesh.GetNumLocalElements(), expected_elements);

        //Note that boundary nodes are local to the process
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(constructed_mesh.GetNumBoundaryNodes(), 2u);
        }
        else
        {
            TS_ASSERT_LESS_THAN(constructed_mesh.GetNumBoundaryNodes(), 2u);
        }

    }


    void TestConstructLinearMeshVerySmall()
    {
        DistributedTetrahedralMesh<1,1> small_mesh;
        //Coverage hack
        TS_ASSERT_THROWS_THIS(small_mesh.ConstructLinearMesh(0), "There aren't enough nodes to make parallelisation worthwhile");
        //Works with up to 3 processes
        small_mesh.ConstructLinearMesh(2);
        //Scale it
        small_mesh.Scale(10.0);
        unsigned owned=small_mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        TS_ASSERT_EQUALS(small_mesh.GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(small_mesh.GetNumLocalNodes(), owned);
        TS_ASSERT_EQUALS(small_mesh.GetNumBoundaryElements(),  2u);
        TS_ASSERT_EQUALS(small_mesh.GetNumElements(), 2u);
        //See logic in earlier test
        unsigned expected_elements=owned+1;
        
        std::vector<unsigned> halo_indices;
        small_mesh.GetHaloNodeIndices(halo_indices);
        /**
         * 1 Proc:
         * p0:  0 Ow 1 Ow 2 Ow
         * 2 Proc:
         * p0:  0 Ow 1 Ow 2 Ha
         * p1:       1 Ha 2 Ow
         * 3 Proc:
         * p0:  0 Ow 1 Ha
         * p1:  0 Ha 1 Ow 2 Ha
         * p2:       2 Ha 3 Ow
         */ 
        if (PetscTools::AmMaster())
        {
            expected_elements--;
            //Left processor always owns left node
            TS_ASSERT_EQUALS(small_mesh.GetNode(0),small_mesh.GetAnyNode(0));
            TS_ASSERT_DELTA(small_mesh.GetAnyNode(0)->rGetLocation()[0], 0.0, 1e-5);
            TS_ASSERT_DELTA(small_mesh.GetAnyNode(1)->rGetLocation()[0], 10.0, 1e-5);
            if (PetscTools::IsSequential())
            {
                TS_ASSERT_EQUALS(halo_indices.size(), 0u);
            }
            else
            {
                TS_ASSERT_EQUALS(halo_indices.size(), 1u);
                TS_ASSERT_DELTA(halo_indices[0], 1u, 1u);//Halo is at index 1 (3 procs) or index 2 (2 procs)
            }
        }
        if (PetscTools::AmTopMost() && PetscTools::GetNumProcs() <= 3)
        {
            expected_elements--;
            if (PetscTools::IsSequential())
            {
                TS_ASSERT_EQUALS(halo_indices.size(), 0u);
            }
            else
            {
                TS_ASSERT_EQUALS(halo_indices.size(), 1u);
                TS_ASSERT_EQUALS(halo_indices[0], 1u); //Halo is at index 1 (2 or 3 procs)
                TS_ASSERT_THROWS_CONTAINS(small_mesh.GetAnyNode(0), "Requested node/halo");
                //Right processor has  node 1 as halo
                TS_ASSERT_THROWS_CONTAINS(small_mesh.GetNode(1), "does not belong to processor");
                TS_ASSERT_THROWS_NOTHING(small_mesh.GetAnyNode(1));//It's a halo
                TS_ASSERT_DELTA(small_mesh.GetAnyNode(1)->rGetLocation()[0], 10.0, 1e-5);
                TS_ASSERT_DELTA(small_mesh.GetAnyNode(2)->rGetLocation()[0], 20.0, 1e-5);
            }
        }
        if (PetscTools::GetMyRank() < 3)
        {
            TS_ASSERT_EQUALS(small_mesh.GetNumLocalElements(), expected_elements);
        }
        else
        {
            //This process owns nothing
            TS_ASSERT_EQUALS(small_mesh.GetNumLocalNodes(), 0u);
            TS_ASSERT_EQUALS(owned, 0u);
            TS_ASSERT_EQUALS(small_mesh.GetNumLocalElements(), 0u);
            TS_ASSERT_EQUALS(halo_indices.size(), 0u);
        }    
    }

    void TestConstructLinearMeshSmall()
    {
        unsigned width=2;
        //Works well with exactly 3 processors
        if (PetscTools::GetNumProcs() != width + 1)
        {
            TS_TRACE("This test works with exactly 3 processes.");
            return;
        }
        TetrahedralMesh<1,1> base_mesh;
        base_mesh.ConstructLinearMesh(width);
        TrianglesMeshWriter<1,1> mesh_writer("", "linear");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<1,1> mesh_reader(output_dir+"linear");
        DistributedTetrahedralMesh<1,1> read_mesh(DistributedTetrahedralMesh<1,1>::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<1,1> constructed_mesh;
        constructed_mesh.ConstructLinearMesh(width);

        //Double check
        TS_ASSERT_EQUALS(constructed_mesh.GetNumBoundaryNodes(), read_mesh.GetNumBoundaryNodes());
        CompareParallelMeshOwnership(read_mesh, constructed_mesh);
    }
    void TestConstructRetangularMeshSmall()
    {
        unsigned width=1;
        unsigned height=2;
        //Works well with exactly 3 processors
        if (PetscTools::GetNumProcs() != height + 1)
        {
            TS_TRACE("This test works with exactly 3 processes.");
            return;
        }
        TetrahedralMesh<2,2> base_mesh;
        base_mesh.ConstructRectangularMesh(width, height);
        TrianglesMeshWriter<2,2> mesh_writer("", "rectangle");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader(output_dir+"rectangle");
        DistributedTetrahedralMesh<2,2> read_mesh(DistributedTetrahedralMesh<2,2>::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<2,2> constructed_mesh;
        constructed_mesh.ConstructRectangularMesh(width, height, false);

        TS_ASSERT_EQUALS(constructed_mesh.GetNumBoundaryNodes(), read_mesh.GetNumBoundaryNodes());
        CompareParallelMeshOwnership(read_mesh, constructed_mesh);
    }

    void TestConstructRetangularMesh()
    {
        unsigned width=5;
        unsigned height=4*PetscTools::GetNumProcs()-1; //4*NumProcs layers of nodes (ensure dumb partition works in slices)

        TetrahedralMesh<2,2> base_mesh;
        base_mesh.ConstructRectangularMesh(width, height);
        TrianglesMeshWriter<2,2> mesh_writer("", "rectangle");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader(output_dir+"rectangle");
        DistributedTetrahedralMesh<2,2> read_mesh(DistributedTetrahedralMesh<2,2>::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<2,2> constructed_mesh;
        //Coverage
        TS_ASSERT_THROWS_THIS(constructed_mesh.ConstructRectangularMesh(width, 0, false),
                            "There aren't enough nodes to make parallelisation worthwhile");

        //Real mesh construction
        constructed_mesh.ConstructRectangularMesh(width, height, false);

        CompareParallelMeshOwnership(read_mesh, constructed_mesh);

        if (PetscTools::AmTopMost())
        {
            //Verify some element indices -- top left diagonal goes NW-SE (normal)
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[0],            2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[1], (height-1)+2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[0],            1.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[1], (height-1)+1.0/3.0, 1e-5);
        }
        if (PetscTools::AmMaster())
        {
            //Verify some element indices -- bottom left diagonal goes NW-SE (normal)
            TS_ASSERT_DELTA(constructed_mesh.GetElement(0)->CalculateCentroid()[0], 2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(0)->CalculateCentroid()[1], 2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(1)->CalculateCentroid()[0], 1.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(1)->CalculateCentroid()[1], 1.0/3.0, 1e-5);
        }

    }


    void TestConstructRetangularMeshStagger()
    {
        unsigned width=4;
        unsigned height=4*PetscTools::GetNumProcs()-1; //4*NumProcs layers of nodes (ensure dumb partition works in slices)

        TetrahedralMesh<2,2> base_mesh;
        base_mesh.ConstructRectangularMesh(width, height, true);
        TrianglesMeshWriter<2,2> mesh_writer("", "rectangle");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<2,2> mesh_reader(output_dir+"rectangle");
        DistributedTetrahedralMesh<2,2> read_mesh(DistributedTetrahedralMesh<2,2>::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<2,2> constructed_mesh;
        constructed_mesh.ConstructRectangularMesh(width, height, true);

        CompareParallelMeshOwnership(read_mesh, constructed_mesh);

        if (PetscTools::AmTopMost())
        {
            //Verify some element indices -- top left diagonal goes NW-SE (normal)
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[0],              2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1))->CalculateCentroid()[1],   (height-1)+2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[0],            1.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2*(width)*(height-1)+1)->CalculateCentroid()[1], (height-1)+1.0/3.0, 1e-5);
        }
        if (PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(height%2, 1u);//If height is odd the bottom left is not staggered - next one is
            //Verify some element indices -- bottom left diagonal goes SW-NE (stagger)
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2)->CalculateCentroid()[0], 1 + 1.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(2)->CalculateCentroid()[1], 2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(3)->CalculateCentroid()[0], 1 + 2.0/3.0, 1e-5);
            TS_ASSERT_DELTA(constructed_mesh.GetElement(3)->CalculateCentroid()[1], 1.0/3.0, 1e-5);
        }
    }

    void TestConstructCuboidMesh()
    {
        unsigned width=2;
        unsigned height=3;
        unsigned depth=4*PetscTools::GetNumProcs()-1;

        TetrahedralMesh<3,3> base_mesh;
        base_mesh.ConstructCuboid(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer("", "cuboid");
        mesh_writer.WriteFilesUsingMesh(base_mesh);
        PetscTools::Barrier();
        std::string output_dir = mesh_writer.GetOutputDirectory();
        TrianglesMeshReader<3,3> mesh_reader(output_dir+"cuboid");
        DistributedTetrahedralMesh<3,3> read_mesh(DistributedTetrahedralMesh<3,3>::DUMB);
        read_mesh.ConstructFromMeshReader(mesh_reader);

        DistributedTetrahedralMesh<3,3> constructed_mesh;
        //Coverage
        TS_ASSERT_THROWS_THIS(constructed_mesh.ConstructCuboid(width, height, 0),
                            "There aren't enough nodes to make parallelisation worthwhile");
        constructed_mesh.ConstructCuboid(width, height, depth);

        CompareParallelMeshOwnership(read_mesh, constructed_mesh);
    }

    void TestParallelWriting1D()
    {
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");
        TetrahedralMesh<1,1> sequential_mesh;
        sequential_mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<1,1> mesh_writer1("TestDistributedMeshWriter", "seq_line_10_elements");
        mesh_writer1.WriteFilesUsingMesh(sequential_mesh);

        DistributedTetrahedralMesh<1,1> distributed_mesh(DistributedTetrahedralMesh<1,1>::DUMB); //Makes sure that there is no permutation
        AbstractTetrahedralMesh<1,1> *p_distributed_mesh = &distributed_mesh; //Hide the fact that it's distributed from the compiler
        distributed_mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<1,1> mesh_writer2("TestDistributedMeshWriter", "par_line_10_elements", false);
        mesh_writer2.WriteFilesUsingMesh(*p_distributed_mesh);

        std::string output_dir = mesh_writer1.GetOutputDirectory();
        /* Compare
         grep -ab "#" /tmp/$USER/testoutput/TestDistributedMeshWriter/seq_line_10_elements.node
         219:# Generated by Chaste mesh file writer
         258:# Created by Chaste version 1.1.7899 on Thu, 04 Feb 2010 12:55:33 +0000.  Chaste was built on Wed, 03 Feb 2010 18:14:17 +0000 by machine (uname) 'Linux chaste-bob 2.6.24-25-server #1 SMP Tue Oct 20 07:20:02 UTC 2009 x86_64' using settings: default, shared libraries.
         i.e. Bytes after 258 are provenance data and will change -- one second clock difference will be spotted
         grep -ab "#" /tmp/$USER/testoutput/TestDistributedMeshWriter/seq_line_10_elements.ele
         */
        TS_ASSERT_EQUALS(system(("cmp -n 258 " + output_dir + "/par_line_10_elements.node "+ output_dir + "/seq_line_10_elements.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 127 " + output_dir + "/par_line_10_elements.ele "+ output_dir + "/seq_line_10_elements.ele").c_str()), 0);
    }

    void TestParallelWriting3D()
    {
        TrianglesMeshReader<3,3> reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> sequential_mesh;
        sequential_mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<3,3> mesh_writer1("TestDistributedMeshWriter", "seq_cube_2mm_12_elements", false);
        mesh_writer1.WriteFilesUsingMesh(sequential_mesh);

        DistributedTetrahedralMesh<3,3> distributed_mesh(DistributedTetrahedralMesh<3,3>::DUMB); //Makes sure that there is no permutation
        AbstractTetrahedralMesh<3,3> *p_distributed_mesh = &distributed_mesh; //Hide the fact that it's distributed from the compiler

        distributed_mesh.ConstructFromMeshReader(reader);
        TrianglesMeshWriter<3,3> mesh_writer2("TestDistributedMeshWriter", "par_cube_2mm_12_elements", false);
        mesh_writer2.WriteFilesUsingMesh(*p_distributed_mesh);

        std::string output_dir = mesh_writer1.GetOutputDirectory();

        TS_ASSERT_EQUALS(system(("cmp -n 553 " + output_dir + "/par_cube_2mm_12_elements.node "+ output_dir + "/seq_cube_2mm_12_elements.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 206 " + output_dir + "/par_cube_2mm_12_elements.ele "+ output_dir + "/seq_cube_2mm_12_elements.ele").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp -n 226 " + output_dir + "/par_cube_2mm_12_elements.face "+ output_dir + "/seq_cube_2mm_12_elements.face").c_str()), 0);
    }


    void TestArchiveOfConstructedMesh() throw(Exception)
    {
        std::string archive_dir = "archive";
        std::string archive_file = "distributed_rectangle.arch";
        ArchiveLocationInfo::SetMeshFilename("distributed_rectangle");

        DistributedTetrahedralMesh<2,2>* p_mesh = new DistributedTetrahedralMesh<2,2>;
        //std::vector<unsigned> halo_node_indices;
        std::vector<Node<2>*> halo_nodes;
        unsigned num_nodes;
        unsigned local_num_nodes;
        unsigned num_elements;

        unsigned width=4;
        unsigned height=4*PetscTools::GetNumProcs()-1; //4*NumProcs layers of nodes (ensure dumb partition works in slices)
        // archive
        {
            p_mesh->ConstructRectangularMesh(width, height);
            num_nodes = p_mesh->GetNumNodes();
            local_num_nodes = p_mesh->GetNumLocalNodes();
            num_elements = p_mesh->GetNumElements();

            halo_nodes = p_mesh->mHaloNodes;

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractTetrahedralMesh<2,2>* const p_mesh_abstract = static_cast<AbstractTetrahedralMesh<2,2>* >(p_mesh);
            (*p_arch) << p_mesh_abstract;
        }

        // restore
        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<2,2>* p_mesh_abstract2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_mesh_abstract2;
            // Check we have the right number of nodes & elements
            DistributedTetrahedralMesh<2,2>* p_mesh2 = static_cast<DistributedTetrahedralMesh<2,2>*>(p_mesh_abstract2);

            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(), num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumLocalNodes(), local_num_nodes);
            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), num_elements);

            CompareParallelMeshOwnership(*p_mesh, *p_mesh2);
            // Check some node co-ordinates
            try
            {
                Node<2>* p_node1 =p_mesh->GetNode(0);
                Node<2>* p_node2 = p_mesh2->GetNode(0);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], p_node2->GetPoint()[1], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], 0.0, 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Node<2>* p_node1 = p_mesh->GetNode((width+1)*(height+1)-1);
                Node<2>* p_node2 = p_mesh2->GetNode((width+1)*(height+1)-1);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], p_node2->GetPoint()[0], 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[0], (double) width, 1e-6);
                TS_ASSERT_DELTA(p_node1->GetPoint()[1], (double) height, 1e-6);
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check first element has the right nodes
            try
            {
                Element<2,2>* p_element = p_mesh->GetElement(0);
                Element<2,2>* p_element2 = p_mesh2->GetElement(0);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            try
            {
                Element<2,2>* p_element = p_mesh->GetElement(width*height);
                Element<2,2>* p_element2 = p_mesh2->GetElement(500);
                TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), p_element2->GetNodeGlobalIndex(0));
            }
            catch(Exception& e)
            {
                TS_ASSERT_DIFFERS((int)e.GetShortMessage().find("does not belong to processor"),-1);
            }

            // Check the halo nodes are right
            std::vector<Node<2>*> halo_nodes2 = p_mesh2->mHaloNodes;
            TS_ASSERT_EQUALS(halo_nodes2.size(), halo_nodes.size());
            delete p_mesh2;
        }
        delete p_mesh;
    }
    
//    {
//
//        TetrahedralMesh<3,3> cuboid_mesh;
//        cuboid_mesh.ConstructCuboid(20, 20, 20);
//        TS_ASSERT_EQUALS(cuboid_mesh.GetNumNodes(), 9261u); // 21x21x21 nodes
//        TS_ASSERT_EQUALS(cuboid_mesh.GetNumElements(), 48000u);
//        TS_ASSERT_EQUALS(cuboid_mesh.GetNumBoundaryElements(), 4800u);
//        
//        cuboid_mesh.Translate(-10,-10,-10);
//        cuboid_mesh.Scale(0.25/10.0,0.25/10.0,0.25/10.0);
//        
//        TrianglesMeshWriter<3,3> mesh_writer("", "Cube21");
//        mesh_writer.WriteFilesUsingMesh(cuboid_mesh);
//    }


    void TestLoadBadFacesException() throw (Exception)
    {


        DistributedTetrahedralMesh<3,3> distributed_mesh_bad;
        TrianglesMeshReader<3,3> mesh_reader_bad("mesh/test/data/cube_21_nodes_side/Cube21_bad_faces"); // 5x5x5mm cube (internode distance = 0.25mm)
        if (PetscTools::IsSequential())
        {
            distributed_mesh_bad.ConstructFromMeshReader(mesh_reader_bad);
            TS_ASSERT_EQUALS(distributed_mesh_bad.GetNumNodes(), 9261u); // 21x21x21 nodes
            TS_ASSERT_EQUALS(distributed_mesh_bad.GetNumElements(), 48000u);
            TS_ASSERT_EQUALS(distributed_mesh_bad.GetNumBoundaryElements(), 4800u);
        }
        else
        {
            TS_ASSERT_THROWS_CONTAINS(distributed_mesh_bad.ConstructFromMeshReader(mesh_reader_bad), "Face does not appear in element file (Face ");
        }
        
        DistributedTetrahedralMesh<3,3> distributed_mesh_good;
        TrianglesMeshReader<3,3> mesh_reader_good("mesh/test/data/cube_21_nodes_side/Cube21"); // 5x5x5mm cube (internode distance = 0.25mm)
        distributed_mesh_good.ConstructFromMeshReader(mesh_reader_good);
        TS_ASSERT_EQUALS(distributed_mesh_good.GetNumNodes(), 9261u); // 21x21x21 nodes
        TS_ASSERT_EQUALS(distributed_mesh_good.GetNumElements(), 48000u);
        TS_ASSERT_EQUALS(distributed_mesh_good.GetNumBoundaryElements(), 4800u);
    }    

};
#endif /*TESTDISTRIBUTEDTETRAHEDRALMESH_HPP_*/
