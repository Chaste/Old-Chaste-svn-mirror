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
#ifndef _TESTQUADRATICMESH_HPP_
#define _TESTQUADRATICMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "QuadraticMesh.hpp"

class TestQuadraticMesh : public CxxTest::TestSuite
{
public:

    void TestQuadraticMesh1d() throw(Exception)
    {
        QuadraticMesh<1> mesh("mesh/test/data/1D_0_to_1_10_elements_quadratic", false);
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
        QuadraticMesh<2> mesh("mesh/test/data/square_128_elements_quadratic", false);

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
        QuadraticMesh<3> mesh("mesh/test/data/3D_Single_tetrahedron_element_quadratic", false);
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
        }

        // Lots of internal and boundary nodes in this mesh..
        QuadraticMesh<3> mesh2("mesh/test/data/cube_1626_elements_quadratic", false);

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
        QuadraticMesh<2> mesh(1.0, 1.0, 1, 1);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);
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
        QuadraticMesh<2> mesh(3.14159, 2.71828183, 10, 10);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 21*21u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 200u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 121u);

        // Each element should have 6 nodes
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 6u);
        }
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetBoundaryElement(i)->GetNumNodes(), 3u);
        }

        TS_ASSERT_DELTA( mesh.GetNode(120)->rGetLocation()[0], 3.14159, 1e-4);
        TS_ASSERT_DELTA( mesh.GetNode(120)->rGetLocation()[1], 2.71828183, 1e-5);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 40u);
    }

    void TestAutomaticallyGenerated3dMeshSimple() throw(Exception)
    {
        QuadraticMesh<3> mesh(3.14159, 2.71828183, 2.99792 /* c! */, 1, 1, 1);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 27u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6u);
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 8u); // 6^3 = 216

        // Each element should have 10 nodes
        for (unsigned i=1; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> x = mesh.GetNode(i)->rGetLocation();

            // Check the extra nodes aren't (0,0,0).
            // This fails with 32bit outdated binary.
            TS_ASSERT_LESS_THAN(1e-12, norm_2(x)); // assert x not equal to 0
        }
    }

    void TestAutomaticallyGenerated3dMesh() throw(Exception)
    {
        QuadraticMesh<3> mesh(3.14159, 2.71828183, 2.99792 /* c! */, 5, 5, 5);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 1331u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 750u); // 5 cubes in each direction = 125 cubes => 125 x 6 tetrahedra per cube = 750
        TS_ASSERT_EQUALS(mesh.GetNumVertices(), 216u); // 6^3 = 216

        // Each element should have 10 nodes
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetNumNodes(), 10u);
        }

        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[0], 3.14159, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[1], 2.71828183, 1e-5);
        TS_ASSERT_DELTA(mesh.GetNode(215)->rGetLocation()[2], 2.99792, 1e-4);
    }
    

    
    void TestWritingReadingBoundaryElementsFile2d() throw(Exception)
    {
        // read in a quadratic mesh with linear boundary elements
        QuadraticMesh<2> mesh1("mesh/test/data/square_128_elements_quadratic", false);

        // write the computed quadratic boundary elements
        mesh1.WriteBoundaryElementFile("TestQuadraticMesh","generated2d.face");

        OutputFileHandler handler("TestQuadraticMesh", false);
        std::string file = handler.GetOutputDirectoryFullPath() + "/generated2d.face";
        TS_ASSERT_EQUALS(system(("diff " + file + " mesh/test/data/square_128_elements_fully_quadratic.edge").c_str()), 0);
 
        // read in the same quadratic mesh with /quadratic/ boundary elements 
        QuadraticMesh<2> mesh2("mesh/test/data/square_128_elements_fully_quadratic", true);
        
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
            
            for(unsigned i=0; i<3; i++)
            {
               TS_ASSERT_EQUALS( (*iter1)->GetNodeGlobalIndex(i), (*iter2)->GetNodeGlobalIndex(i));
            }
            iter1++;
        }
    }  

    void TestWritingReadingBoundaryElementsFile3d() throw(Exception)
    {
        // read in a quadratic mesh with linear boundary elements
        QuadraticMesh<3> mesh1("mesh/test/data/cube_1626_elements_quadratic", false);

        // write the computed quadratic boundary elements
        mesh1.WriteBoundaryElementFile("TestQuadraticMesh","generated3d.face");

        OutputFileHandler handler("TestQuadraticMesh", false);
        std::string file = handler.GetOutputDirectoryFullPath() + "/generated3d.face";
        TS_ASSERT_EQUALS(system(("diff " + file + " mesh/test/data/cube_1626_elements_fully_quadratic.face").c_str()), 0);
 
        // read in the same quadratic mesh with /quadratic/ boundary elements 
        QuadraticMesh<3> mesh2("mesh/test/data/cube_1626_elements_fully_quadratic", true);
        
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

    void TestExceptions() throw(Exception)
    {
        TS_ASSERT_THROWS_ANYTHING(QuadraticMesh<1> mesh("mesh/test/data/baddata/bad_1D_0_to_1_10_elements_quadratic", false));
    }
};

#endif // _TESTQUADRATICMESH_HPP_
