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

#include "MixedTetrahedralMesh.hpp"

#include <sstream>
#include "Exception.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MixedTetrahedralMesh()
    : TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>()
{
    mpCoarseFineNodeMap = NULL;
    mpFineCoarseNodeMap = NULL;
    mpFineMesh = NULL;
    mAllocatedFineMeshMemory = false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~MixedTetrahedralMesh()
{
    delete mpCoarseFineNodeMap;
    delete mpFineCoarseNodeMap;
    if (mAllocatedFineMeshMemory)
    {
        delete mpFineMesh;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetFineMesh(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pFineMesh)
{
    if (mpFineMesh != NULL)
    {
        EXCEPTION("SetFineMesh can be called at most once.");
    }
    mpFineMesh=pFineMesh;

    assert(this->GetNumNodes() > 0);
    mpCoarseFineNodeMap = new NodeMap(this->GetNumNodes());
    mpFineCoarseNodeMap = new NodeMap(mpFineMesh->GetNumNodes());

    //Construct sorted lists of nodes from each mesh
    std::vector<Node<SPACE_DIM> * > coarse_nodes;
    for (unsigned i=0;i<this->GetNumNodes();i++)
    {
        coarse_nodes.push_back(this->GetNode(i));
    }
    std::sort(coarse_nodes.begin(), coarse_nodes.end(), CompareNodesLex());

    std::vector<Node<SPACE_DIM> * > fine_nodes;
    for (unsigned i=0;i<mpFineMesh->GetNumNodes();i++)
    {
        fine_nodes.push_back(mpFineMesh->GetNode(i));
    }
    std::sort(fine_nodes.begin(), fine_nodes.end(), CompareNodesLex());


    unsigned fine_mesh_index=0;
    //.. update node map
    for (unsigned coarse_mesh_index=0;coarse_mesh_index<this->GetNumNodes();coarse_mesh_index++)
    {
        while (!EqualNodes(fine_nodes[fine_mesh_index], coarse_nodes[coarse_mesh_index]) )
        {
            fine_mesh_index++;
            if (fine_mesh_index==fine_nodes.size())
            {
                std::stringstream exception_string;
                exception_string << "Coarse mesh node: "
                                 << coarse_nodes[coarse_mesh_index]->GetIndex()
                                 << " doesn't have a partner in the fine mesh.";
                EXCEPTION(exception_string.str());
            }
        }
        //Same node, set maps
        mpCoarseFineNodeMap->SetNewIndex(coarse_mesh_index, fine_mesh_index);
        mpFineCoarseNodeMap->SetNewIndex(fine_mesh_index, coarse_mesh_index);
    }

    //Calculate a map from fine nodes to coarse elements
    mFineNodeToCoarseElementMap.resize(mpFineMesh->GetNumNodes());
    for (unsigned i=0; i<mpFineMesh->GetNumNodes(); i++)
    {
        //Find a representative coarse element for this node and put it in the map
        unsigned coarse_element_index;
        try
        {
            coarse_element_index = this->GetContainingElementIndex(mpFineMesh->GetNode(i)->GetPoint());
        }
        catch (Exception &e)
        {
            // find nearest coarse element
            coarse_element_index = this->GetNearestElementIndex(mpFineMesh->GetNode(i)->GetPoint());
        }
        mFineNodeToCoarseElementMap[i]=this->GetElement(coarse_element_index);
    }

//    //Calculate the map from coarse elements to fine ones
//    mCoarseFineElementsMap.resize(this->GetNumElements());
//    typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator i_fine_element;
//    for (i_fine_element = mpFineMesh->GetElementIteratorBegin();
//         i_fine_element != mpFineMesh->GetElementIteratorEnd();
//         i_fine_element++)
//    {
//        bool coarse_elements_found=false;
//        for (unsigned node_local_index = 0; node_local_index <= ELEMENT_DIM; node_local_index++)
//        {
//            ChastePoint<SPACE_DIM> vertex = (*i_fine_element)->GetNode(node_local_index)->GetPoint();
//            try
//            {
//                unsigned coarse_element_index = this->GetContainingElementIndex(vertex, true);
//                coarse_elements_found=true;
//                mCoarseFineElementsMap[coarse_element_index].insert(*i_fine_element);
//            }
//            catch (Exception &e)
//            {
//                // vertex must coincide with coarse node
//            }
//        }
//        if (!coarse_elements_found)
//        {
//            // fine element must coincide with coarse element
//            ChastePoint<SPACE_DIM> centroid = ChastePoint<SPACE_DIM>((*i_fine_element)->CalculateCentroid());
//            try
//            {
//                unsigned coarse_element_index = this->GetContainingElementIndex(centroid, true);
//                mCoarseFineElementsMap[coarse_element_index].insert(*i_fine_element);
//            }
//            #define COVERAGE_IGNORE
//            catch (Exception &e)
//            {
//                //EXCEPTION("Fine mesh contains an element which does not overlap any coarse mesh element");
//            }
//            #undef COVERAGE_IGNORE
//        }
//
//    }
//    std::cout<<"Done\n";
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMeshes(double width,
                                                                              double height,
                                                                              unsigned numCoarseElemInEachDirection,
                                                                              unsigned numFineElemInEachDirection)
{
    assert(ELEMENT_DIM==2);
    assert(SPACE_DIM==2);
    assert(numCoarseElemInEachDirection <= numFineElemInEachDirection); // probably entered wrong way round if not

    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_fine_mesh = new TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>;
    p_fine_mesh->ConstructRectangularMesh(numFineElemInEachDirection, numFineElemInEachDirection, false);
    p_fine_mesh->Scale(width/numFineElemInEachDirection, height/numFineElemInEachDirection, 0.0);

    this->ConstructRectangularMesh(numCoarseElemInEachDirection, numCoarseElemInEachDirection, false);
    this->Scale(width/numCoarseElemInEachDirection, height/numCoarseElemInEachDirection, 0.0);

    SetFineMesh(p_fine_mesh);

    mAllocatedFineMeshMemory = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MixedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructCuboidMeshes(double width,
                                                                         double height,
                                                                         double depth,
                                                                         unsigned numCoarseElemInEachDirection,
                                                                         unsigned numFineElemInEachDirection)
{
    assert(ELEMENT_DIM==3);
    assert(SPACE_DIM==3);
    assert(numCoarseElemInEachDirection <= numFineElemInEachDirection); // probably entered wrong way round if not

    TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_fine_mesh = new TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>;

    p_fine_mesh->ConstructCuboid(numFineElemInEachDirection, numFineElemInEachDirection, numFineElemInEachDirection, false);
    p_fine_mesh->Scale(width/numFineElemInEachDirection, height/numFineElemInEachDirection, depth/numFineElemInEachDirection);

    this->ConstructCuboid(numCoarseElemInEachDirection, numCoarseElemInEachDirection, numCoarseElemInEachDirection, false);
    this->Scale(width/numCoarseElemInEachDirection, height/numCoarseElemInEachDirection, depth/numCoarseElemInEachDirection);

    SetFineMesh(p_fine_mesh);

    mAllocatedFineMeshMemory = true;
}



/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

//template class MixedTetrahedralMesh<1,1>;
template class MixedTetrahedralMesh<2,2>;
template class MixedTetrahedralMesh<3,3>;
