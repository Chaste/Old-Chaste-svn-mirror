/*

Copyright (C) University of Oxford, 2005-2011

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
#ifndef TESTNODESONLYMESH_HPP_
#define TESTNODESONLYMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <algorithm>

#include "UblasCustomFunctions.hpp"
#include "NodesOnlyMesh.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "ArchiveOpener.hpp"
#include "TrianglesMeshWriter.hpp"
#include "VtkMeshWriter.hpp"

class TestNodesOnlyMesh : public CxxTest::TestSuite
{
public:
    void TestConstructNodesWithoutMesh()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 8u);
        
        //When the mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
    void TestWriteNodesWithoutMesh()
    {
 #ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);
        
        //Note in my version of Paraview, you need data on points before you
        //can view with Glyphs.
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "just_nodes", false);



        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add boundary node "point" data
        std::vector<double> boundary;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            boundary.push_back(mesh.GetNode(i)->IsBoundaryNode());
        }
        writer.AddPointData("Boundary", boundary);
        
        // Add fibre type to "point" data
        std::vector< c_vector<double, 3> > location;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            location.push_back(mesh.GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        writer.WriteFilesUsingMesh(mesh);

        
        //When the mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
#endif //CHASTE_VTK
    }


};

#endif /*TESTNODESONLYMESH_HPP_*/
