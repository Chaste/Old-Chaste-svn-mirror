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
#ifndef _TESTQUADRATICMESH_HPP_
#define _TESTQUADRATICMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "QuadraticMesh.hpp"
#include "OutputFileHandler.hpp"
#include "ArchiveOpener.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestQuadraticMesh : public CxxTest::TestSuite
{
public:

    void TestQuadraticMesh1d() throw(Exception)
    {
        QuadraticMesh<1> mesh;
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);

        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11u);

        // Node 2 (ie middle) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 11u);
        TS_ASSERT_DELTA(mesh.GetNode(11)->rGetLocation()[0], 0.05, 1e-12);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            bool is_boundary_node = mesh.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node, ((x==0)||(x==1)));
        }
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            // Check internal nodes have corrent element associated with them
            std::set<unsigned> internal_node_elems;
            internal_node_elems.insert(mesh.GetElement(i)->GetIndex());
            TS_ASSERT_EQUALS(internal_node_elems,mesh.GetElement(i)->GetNode(2)->rGetContainingElementIndices());
        }
    }


    // identical to above, except mesh is generated not read
    void TestQuadraticMesh1dAutomaticallyGenerated() throw(Exception)
    {
        QuadraticMesh<1> mesh(0.1, 1.0);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);

        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11u);

        // Node 2 (ie middle) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 11u);
        TS_ASSERT_DELTA(mesh.GetNode(11)->rGetLocation()[0], 0.05, 1e-12);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            bool is_boundary_node = mesh.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node, ((x==0)||(x==1)));
        }
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            // Check internal nodes have corrent element associated with them
            std::set<unsigned> internal_node_elems;
            internal_node_elems.insert(mesh.GetElement(i)->GetIndex());
            TS_ASSERT_EQUALS(internal_node_elems,mesh.GetElement(i)->GetNode(2)->rGetContainingElementIndices());
        }
    }

    void TestQuadraticMesh2d() throw(Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1, false);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 289u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 128u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 81u);

        // Each element should have 6 nodes
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 6u);

            for (unsigned j=0; j<2; j++)
            {
                // Check internal nodes have corrent element associated with them
                TS_ASSERT(mesh.GetElement(i)->GetNode(j+3)->GetNumContainingElements() <= 2u);
                TS_ASSERT(mesh.GetElement(i)->GetNode(j+3)->GetNumContainingElements() > 0u);

                std::set<unsigned> current_node_indices = mesh.GetElement(i)->GetNode(j)->rGetContainingElementIndices();
                TS_ASSERT_EQUALS(current_node_indices.count(mesh.GetElement(i)->GetIndex()), 1u);

                current_node_indices = mesh.GetElement(i)->GetNode(j+3)->rGetContainingElementIndices();
                TS_ASSERT_EQUALS(current_node_indices.count(mesh.GetElement(i)->GetIndex()), 1u);
            }
        }

        // Node 3 (ie fourth) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(3), 82u);
        // Node 4 (ie fifth) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(4), 83u);
        // Node 5 (ie last) of element 0
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(5), 81u);

        // Each boundary element should have three nodes
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            TS_ASSERT_EQUALS((*iter)->GetNumNodes(), 3u);
        }

        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();

        // The first edge has nodes 53 and 0, according to the edge file...
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(0), 53u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(1), 0u);
        // ...the midnode has to be computed (found) by the QuadraticMesh class
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(2), 81u);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            bool is_boundary_node = mesh.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node, ((x==0)||(x==1)||(y==0)||(y==1)));
        }
    }

    void TestQuadraticMesh3d() throw(Exception)
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1, false);
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);

        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 4u);

        // Check getting global numbers of nodes 4-9 (in non-vertices)
        for (unsigned i=4; i<10; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), i);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(mesh.GetNode(i)->GetNumContainingElements(), 1u);
        }

        // Lots of internal and boundary nodes in this mesh..
        QuadraticMesh<3> mesh2;
        TrianglesMeshReader<3,3> mesh_reader2("mesh/test/data/cube_1626_elements_quadratic",2,1, false);
        mesh2.ConstructFromMeshReader(mesh_reader2);

        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), 2570u);
        TS_ASSERT_EQUALS(mesh2.GetNumElements(), 1626u);
        TS_ASSERT_EQUALS(mesh2.GetNumVertices(), 375u);

        // Each element should have 10 nodes
        for (unsigned i=0; i<mesh2.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh2.GetElement(i)->GetNumNodes(), 10u);

            for (unsigned j=3; j<9; j++)
            {
                // Check internal nodes have corrent element associated with them
                TS_ASSERT(mesh2.GetElement(i)->GetNode(j)->GetNumContainingElements() > 0u);
                std::set<unsigned> current_node_indices = mesh2.GetElement(i)->GetNode(j)->rGetContainingElementIndices();
                TS_ASSERT_EQUALS(current_node_indices.count(mesh2.GetElement(i)->GetIndex()),1u);
            }
        }

        // Each boundary element should have 6 nodes
        for (TetrahedralMesh<3,3>::BoundaryElementIterator iter= mesh2.GetBoundaryElementIteratorBegin();
             iter != mesh2.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_EQUALS((*iter)->GetNumNodes(), 6u);
        }

        TetrahedralMesh<3,3>::BoundaryElementIterator iter = mesh2.GetBoundaryElementIteratorBegin();

        // The first boundary elem has these nodes, according to the edge file..
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(0), 177u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(1), 43u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(2), 85u);
        // .. the internal nodes have to be computed (found) by the QuadraticMesh.
        // The nodes 177,43,85 are all in the third element in the ele file, and
        // they are nodes 1,3,2 respectively. Therefore, the internals are the local
        // nodes 9,5,8 respectively (look the the ordering picture), so..
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(3), 392u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(4), 388u);
        TS_ASSERT_EQUALS( (*iter)->GetNodeGlobalIndex(5), 391u);

        for (unsigned i=0; i<mesh2.GetNumNodes(); i++)
        {
            double x = mesh2.GetNode(i)->rGetLocation()[0];
            double y = mesh2.GetNode(i)->rGetLocation()[1];
            double z = mesh2.GetNode(i)->rGetLocation()[2];
            bool is_boundary_node = mesh2.GetNode(i)->IsBoundaryNode();
            TS_ASSERT_EQUALS(is_boundary_node,  ((x==0)||(x==1)||(y==0)||(y==1)||(z==0)||(z==1)));
        }
    }

    void TestAutomaticallyGenerated2dMesh1() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1.0);

        TS_ASSERT_THROWS_CONTAINS(QuadraticMesh<2> bad_mesh(0.645, 1.0, 1.0), "does not divide");

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 8u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 4u);

        // Each element should have 6 nodes and a valid Jacobian
        double det;
        c_matrix<double, 2, 2> jacob;
        c_matrix<double, 2, 2> inv;

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 6u);
            mesh.GetInverseJacobianForElement(i, jacob, det, inv);
            TS_ASSERT_EQUALS(det, 1.0);
        }

        //Test vertex containment
        TS_ASSERT_EQUALS(mesh.GetNode(0)->GetNumContainingElements(), 1u); //(0,0)
        TS_ASSERT_EQUALS(mesh.GetNode(1)->GetNumContainingElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->GetNumContainingElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNode(3)->GetNumContainingElements(), 1u); //(1,1)
        
       //Test internal node containment
        TS_ASSERT_EQUALS(mesh.GetNode(4)->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNode(5)->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNode(6)->GetNumContainingElements(), 2u); //(.5,.5)
        TS_ASSERT_EQUALS(mesh.GetNode(7)->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNode(8)->GetNumContainingElements(), 1u);
        
        TS_ASSERT_DELTA( mesh.GetNode(3)->rGetLocation()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA( mesh.GetNode(3)->rGetLocation()[1], 1.0, 1e-6);

        // Test boundary elements
        unsigned num_boundary_elements=0;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_EQUALS((*iter)->GetNumNodes(), 3u);

            bool all_x_zero =     (fabs((*iter)->GetNode(0)->rGetLocation()[0])<1e-6)
                               && (fabs((*iter)->GetNode(1)->rGetLocation()[0])<1e-6)
                               && (fabs((*iter)->GetNode(2)->rGetLocation()[0])<1e-6);

            bool all_x_one  =     (fabs((*iter)->GetNode(0)->rGetLocation()[0] - 1.0)<1e-6)
                               && (fabs((*iter)->GetNode(1)->rGetLocation()[0] - 1.0)<1e-6)
                               && (fabs((*iter)->GetNode(2)->rGetLocation()[0] - 1.0)<1e-6);

            bool all_y_zero =     (fabs((*iter)->GetNode(0)->rGetLocation()[1])<1e-6)
                               && (fabs((*iter)->GetNode(1)->rGetLocation()[1])<1e-6)
                               && (fabs((*iter)->GetNode(2)->rGetLocation()[1])<1e-6);

            bool all_y_one  =     (fabs((*iter)->GetNode(0)->rGetLocation()[1] - 1.0)<1e-6)
                               && (fabs((*iter)->GetNode(1)->rGetLocation()[1] - 1.0)<1e-6)
                               && (fabs((*iter)->GetNode(2)->rGetLocation()[1] - 1.0)<1e-6);

            TS_ASSERT_EQUALS(true, all_x_zero || all_x_one || all_y_zero || all_y_one);
            num_boundary_elements++;
        }
        TS_ASSERT_EQUALS(num_boundary_elements, 4u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4u);
        
    }

    void TestAutomaticallyGenerated2dMesh2() throw(Exception)
    {
        QuadraticMesh<2> mesh(3.14159/10,  3.14159, 3.14159/2);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21*11u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 60u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 100u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11*6u);

        // Each element should have 6 nodes
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 6u);
        }
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetBoundaryElement(i)->GetNumNodes(), 3u);
        }

        TS_ASSERT_DELTA( mesh.GetNode(65)->rGetLocation()[0], 3.14159, 1e-4);
        TS_ASSERT_DELTA( mesh.GetNode(65)->rGetLocation()[1], 3.14159/2, 1e-5);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 30u);
    }

    void TestAutomaticallyGenerated3dMeshSimple() throw(Exception)
    {
        double h = 3.14159;
        double width = h;

        QuadraticMesh<3> mesh(h, width, 2*width, 3*width);

        TS_ASSERT_THROWS_CONTAINS(QuadraticMesh<3> bad_mesh(0.645, 1.0, 1.0, 1.0), "does not divide");


        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3*5*7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6*1*2*3u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 2*3*4u); 
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 90u);

        for (unsigned i=1; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> x = mesh.GetNode(i)->rGetLocation();

            // Check the extra nodes aren't (0,0,0).
            // This fails with 32bit outdated binary.
            TS_ASSERT_LESS_THAN(1e-12, norm_2(x)); // assert x not equal to 0
        }
        
        TS_ASSERT_DELTA( mesh.GetNode(23)->rGetLocation()[0], width, 1e-8);
        TS_ASSERT_DELTA( mesh.GetNode(23)->rGetLocation()[1], 2*width, 1e-8);
        TS_ASSERT_DELTA( mesh.GetNode(23)->rGetLocation()[2], 3*width, 1e-8);
        
        
        // second 1 by 1 by 1 mesh
        QuadraticMesh<3> mesh2(h, width, width, width);

        for (unsigned i=1; i<mesh2.GetNumNodes(); i++)
        {
            //Check that all nodes have containg elements
            TS_ASSERT_LESS_THAN(0u, mesh2.GetNode(i)->GetNumContainingElements());
            //Mid-point of cube will have access to all 6 elements
            TS_ASSERT_LESS_THAN_EQUALS(mesh2.GetNode(i)->GetNumContainingElements(), 6u); 
        }

        TS_ASSERT_EQUALS(mesh2.CalculateMaximumContainingElementsPerProcess(), 6U); //The midpoint, as given above
        //There are  8 vertex nodes in the cube
        //There are 12 internal nodes on the cube edges
        //There are  6 internal nodes on the diagonals to the cube faces
        //There is   1 interal node on the separating diagonal
        // 8V + 19I = 27 nodes 
        TS_ASSERT_EQUALS(mesh2.CalculateMaximumNodeConnectivityPerProcess(), 27U); //The midpoint, as given above
    }

    void TestAutomaticallyGenerated3dMesh() throw(Exception)
    {
        QuadraticMesh<3> mesh(0.5,  2.5, 2.5, 2.5);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 1331u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 750u); // 5 cubes in each direction = 125 cubes => 125 x 6 tetrahedra per cube = 750
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 216u); // 6^3 = 216
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 602u);

        // Each element should have 10 nodes
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 10u);
        }

        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[0], 2.5, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[1], 2.5, 1e-5);
        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[2], 2.5, 1e-4);
        TS_ASSERT_EQUALS(mesh.CalculateMaximumContainingElementsPerProcess(), 24U); //Four surrounding cubes may have all 6 tetrahedra meeting at a node
        TS_ASSERT_EQUALS(mesh.CalculateMaximumNodeConnectivityPerProcess(), 65U);
    }



    void TestWritingReadingBoundaryElementsFile2d() throw(Exception)
    {
        // read in a quadratic mesh with linear boundary elements
        QuadraticMesh<2> mesh1;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh1.ConstructFromMeshReader(mesh_reader);

        // write the computed quadratic boundary elements
        mesh1.WriteBoundaryElementFile("TestQuadraticMesh","generated2d.face");

        OutputFileHandler handler("TestQuadraticMesh", false);
        std::string file = handler.GetOutputDirectoryFullPath() + "/generated2d.face";
        TS_ASSERT_EQUALS(system(("diff " + file + " mesh/test/data/square_128_elements_fully_quadratic.edge").c_str()), 0);

        // read in the same quadratic mesh with /quadratic/ boundary elements
        QuadraticMesh<2> mesh2;
        TrianglesMeshReader<2,2> mesh_reader2("mesh/test/data/square_128_elements_fully_quadratic",2,2,false);
        mesh2.ConstructFromMeshReader(mesh_reader2);

        // compare the boundary elements of both meshes, should be identical (as one was created from the other)
        QuadraticMesh<2>::BoundaryElementIterator iter1
               = mesh1.GetBoundaryElementIteratorBegin();

        for (QuadraticMesh<2>::BoundaryElementIterator iter2
               = mesh2.GetBoundaryElementIteratorBegin();
             iter2 != mesh2.GetBoundaryElementIteratorEnd();
             ++iter2)
        {
            TS_ASSERT_EQUALS( (*iter1)->GetNumNodes(), 3u );
            TS_ASSERT_EQUALS( (*iter2)->GetNumNodes(), 3u );

            for (unsigned i=0; i<3; i++)
            {
               TS_ASSERT_EQUALS( (*iter1)->GetNodeGlobalIndex(i), (*iter2)->GetNodeGlobalIndex(i));
            }
            iter1++;
        }
    }

    void TestWritingReadingBoundaryElementsFile3d() throw(Exception)
    {
        // read in a quadratic mesh with linear boundary elements
        QuadraticMesh<3> mesh1;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements_quadratic",2,1,false);
        mesh1.ConstructFromMeshReader(mesh_reader);

        // write the computed quadratic boundary elements
        mesh1.WriteBoundaryElementFile("TestQuadraticMesh","generated3d.face");

        OutputFileHandler handler("TestQuadraticMesh", false);
        std::string file = handler.GetOutputDirectoryFullPath() + "/generated3d.face";
        TS_ASSERT_EQUALS(system(("diff " + file + " mesh/test/data/cube_1626_elements_fully_quadratic.face").c_str()), 0);

        // read in the same quadratic mesh with /quadratic/ boundary elements
        QuadraticMesh<3> mesh2;
        TrianglesMeshReader<3,3> mesh_reader2("mesh/test/data/cube_1626_elements_fully_quadratic",2,2,false);
        mesh2.ConstructFromMeshReader(mesh_reader2);

        // compare the boundary elements of both meshes, should be identical (as one was created from the other)
        QuadraticMesh<3>::BoundaryElementIterator iter1
               = mesh1.GetBoundaryElementIteratorBegin();

        for (QuadraticMesh<3>::BoundaryElementIterator iter2
               = mesh2.GetBoundaryElementIteratorBegin();
             iter2 != mesh2.GetBoundaryElementIteratorEnd();
             ++iter2)
        {
            TS_ASSERT_EQUALS( (*iter1)->GetNumNodes(), 6u );
            TS_ASSERT_EQUALS( (*iter2)->GetNumNodes(), 6u );

            for(unsigned i=0; i<6; i++)
            {
               TS_ASSERT_EQUALS( (*iter1)->GetNodeGlobalIndex(i), (*iter2)->GetNodeGlobalIndex(i));
            }
            iter1++;
        }
    }

    void TestWritingReadingBoundaryElementsWithContainingElementInfo() throw(Exception)
    {
        // This mesh has quadratic node and ele files, a linear face file that has containing element
        // info
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_152_elements_v3",2,1,true);
        mesh.ConstructFromMeshReader(mesh_reader);

        for (QuadraticMesh<3>::BoundaryElementIterator iter
               = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            TS_ASSERT_EQUALS( (*iter)->GetNumNodes(), 6u );
        }
    }


    void TestExceptions() throw(Exception)
    {

        QuadraticMesh<1> mesh;

        //Bad data
        TrianglesMeshReader<1,1> mesh_reader1("mesh/test/data/baddata/bad_1D_0_to_1_10_elements_quadratic",2,1, false);
        TS_ASSERT_THROWS_THIS(mesh.ConstructFromMeshReader(mesh_reader1),
                "The quadratic mesh doesn\'t appear to have all vertices before the rest of the nodes");


        //Linear mesh
        TrianglesMeshReader<1,1> mesh_reader2("mesh/test/data/1D_0_to_1_10_elements");
        TS_ASSERT_THROWS_THIS(mesh.ConstructFromMeshReader(mesh_reader2),
                "Supplied mesh reader is reading a linear mesh into quadratic mesh");

    }

    void TestArchiving() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "quadratic_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("quadratic_mesh");

        AbstractTetrahedralMesh<3,3>* const p_mesh = new QuadraticMesh<3>;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements_fully_quadratic", 2, 2, false);
        static_cast<QuadraticMesh<3>*>(p_mesh)->ConstructFromMeshReader(mesh_reader);

        {
            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) << p_mesh;
        }

        {
            // Should archive the most abstract class you can to check boost knows what individual classes are.
            // (but here AbstractMesh doesn't have the methods below).
            AbstractTetrahedralMesh<3,3>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_mesh2;

            // compare the boundary elements of both meshes, should be identical (as one was created from the other)
            QuadraticMesh<3>::BoundaryElementIterator iter1
                   = p_mesh->GetBoundaryElementIteratorBegin();

            for (QuadraticMesh<3>::BoundaryElementIterator iter2
                   = p_mesh2->GetBoundaryElementIteratorBegin();
                 iter2 != p_mesh2->GetBoundaryElementIteratorEnd();
                 ++iter2)
            {
                TS_ASSERT_EQUALS( (*iter1)->GetNumNodes(), 6u );
                TS_ASSERT_EQUALS( (*iter2)->GetNumNodes(), 6u );

                for(unsigned i=0; i<6; i++)
                {
                   TS_ASSERT_EQUALS( (*iter1)->GetNodeGlobalIndex(i), (*iter2)->GetNodeGlobalIndex(i));
                }
                iter1++;
            }

            delete p_mesh2;
        }
        delete p_mesh;
    }

    void TestConstructRegularSlabMesh_Directly_1d() throw(Exception)
    {
        QuadraticMesh<1> mesh;
        TS_ASSERT_THROWS_THIS(mesh.ConstructRegularSlabMesh(0.75, 1.0), "Space step does not divide the size of the mesh");
        mesh.ConstructRegularSlabMesh(0.1, 1.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2u);
        
        for(unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], (i+0.0)/10, 1e-8);

            bool is_boundary = (i==0 || i+1==mesh.GetNumVertices());
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsBoundaryNode(), is_boundary);

            std::set<unsigned> containing_elems = mesh.GetNode(i)->rGetContainingElementIndices();
            TS_ASSERT_EQUALS(containing_elems.size(), (is_boundary ? 1u : 2u));
        }


        for(unsigned i=mesh.GetNumVertices(); i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], (i-11.0)/10 + 0.05, 1e-8);

            std::set<unsigned> containing_elems = mesh.GetNode(i)->rGetContainingElementIndices();
            TS_ASSERT_EQUALS(containing_elems.size(), 1u);
        }
    }

    
    void TestConstructRegularSlabMesh_Directly_2d() throw(Exception)
    {
        QuadraticMesh<2> mesh;
        mesh.ConstructRegularSlabMesh(0.1, 1.0, 2.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21*41u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 11*21u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2*10*20u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 120u);
    }

    void TestConstructRegularSlabMesh_Directly_3d() throw(Exception)
    {
        QuadraticMesh<3> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0, 2.0, 3.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3*5*7u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 2*3*4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6*1*2*3u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 90u);
    }
};

#endif // _TESTQUADRATICMESH_HPP_
