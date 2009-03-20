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

#include "QuadraticMesh.hpp"

template<unsigned DIM>
QuadraticMesh<DIM>::QuadraticMesh(const std::string& fileName)
{
    LoadFromFile(fileName);
}

template<unsigned DIM>
QuadraticMesh<DIM>::QuadraticMesh(double xEnd, double yEnd, unsigned numElemX, unsigned numElemY)
{
    assert(DIM==2);
    
    assert(xEnd>0);
    assert(yEnd>0);
    assert(numElemX>0);
    assert(numElemY>0);

    std::string tempfile_name_stem = "temp_quadmesh";

    ////////////////////////////////////////
    // write the node file (vertices only)
    ////////////////////////////////////////
    OutputFileHandler handler("");
    out_stream p_file = handler.OpenOutputFile(tempfile_name_stem+".node");
    
    *p_file << (numElemX+1)*(numElemY+1) << " 2 0 1\n";
    unsigned node_index = 0;
    for (unsigned j=0; j<=numElemY; j++)
    {
        for (unsigned i=0; i<=numElemX; i++)
        {
            double x = xEnd*i/numElemX;
            double y = yEnd*j/numElemY;
            
            bool on_boundary = ( (i==0) || (i==numElemX) || (j==0) || (j==numElemX) );
            *p_file << node_index++ << " " << x << " " << y << " " << (on_boundary?1:0) << "\n";
        }
    }
    p_file->close();
    
    ////////////////////////////////////////////////////////////
    // create the quadratic mesh files using triangle and load
    ////////////////////////////////////////////////////////////
     
    RunMesherAndReadMesh("triangle", handler.GetOutputDirectoryFullPath(), tempfile_name_stem);
}


template<unsigned DIM>
QuadraticMesh<DIM>::QuadraticMesh(double xEnd, double yEnd, double zEnd, 
                                  unsigned numElemX, unsigned numElemY, unsigned numElemZ)
{
    assert(DIM==3);
    
    assert(xEnd>0);
    assert(yEnd>0);
    assert(zEnd>0);
    assert(numElemX>0);
    assert(numElemY>0);
    assert(numElemZ>0);

    std::string tempfile_name_stem = "temp_quadmesh3d";

    ////////////////////////////////////////
    // write the node file (vertices only)
    ////////////////////////////////////////
    OutputFileHandler handler("");
    out_stream p_file = handler.OpenOutputFile(tempfile_name_stem+".node");
    
    *p_file << (numElemX+1)*(numElemY+1)*(numElemZ+1) << " 3 0 0\n";
    unsigned node_index = 0;
    for (unsigned k=0; k<=numElemZ; k++)
    {
        for (unsigned j=0; j<=numElemY; j++)
        {
            for (unsigned i=0; i<=numElemX; i++)
            {
                double x = xEnd*i/numElemX;
                double y = yEnd*j/numElemY;
                double z = zEnd*k/numElemZ; //Not yEnd!
                
                //bool on_boundary = ( (i==0) || (i==numElemX) || (j==0) || (j==numElemY) || (k==0) || (k==numElemZ) );
                *p_file << node_index++ << " " << x << " " << y << " " << z << "\n"; // << (on_boundary?1:0) << "\n";
            }
        }
    }
    p_file->close();
    
    ////////////////////////////////////////////////////////////
    // create the quadratic mesh files using triangle and load
    ////////////////////////////////////////////////////////////
    

    RunMesherAndReadMesh("tetgen", handler.GetOutputDirectoryFullPath(), tempfile_name_stem);
}



template<unsigned DIM>
void QuadraticMesh<DIM>::RunMesherAndReadMesh(std::string binary, 
                                              std::string outputDir, 
                                              std::string fileStem)
{
    // Q = quiet, e = make edge data, o2 = order of elements is 2, ie quadratics
    std::string args = "-Qeo2";
    
    // In 2D we need an edge file. In 3D we need a face file (which is written automatically in Tetgen)
    if (DIM == 3)
    {
        args = "-Qo2";   
    }
        
    std::string command =  binary + " " + args + " " + outputDir
                           + "/" + fileStem + ".node";
    
    if (DIM == 3)
    {
        // Tetgen's quiet mode isn't as quiet as Triangle's
        command += " > /dev/null";
    }
                 
    int return_value = system(command.c_str());
 
    if (return_value != 0)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("Remeshing (by calling " + binary + ") failed.  Do you have it in your path?\n"+
        "The quadratic mesh relies on functionality from triangle (http://www.cs.cmu.edu/~quake/triangle.html) and tetgen (http://tetgen.berlios.de/).");
        #undef COVERAGE_IGNORE 
    }
    
    // move the output files to the chaste directory
    command =   "mv " + outputDir + "/" 
              + fileStem + ".1.* .";
    
    // NOTE: we don't check whether the return value here is zero, because if CHASTE_TESTOUTPUT
    // is "." (ie if it hasn't been exported), then the mv will fail (source and destination files
    // are the same), but this isn't a problem.
    return_value = system(command.c_str());
    
    // load
    LoadFromFile( fileStem + ".1");
    
    // delete the temporary files
    command = "rm -f " + outputDir + "/" + fileStem + ".node";
    EXPECT0(system, command);
    EXPECT0(system, "rm -f " + fileStem + ".1.node");
    EXPECT0(system, "rm -f " + fileStem + ".1.ele");
    
    if (DIM==2) 
    {
        EXPECT0(system, "rm -f " + fileStem + ".1.edge");
    }
    if (DIM==3) 
    {
        EXPECT0(system, "rm -f " + fileStem + ".1.face");
    }
}


template<unsigned DIM>
void QuadraticMesh<DIM>::LoadFromFile(const std::string& fileName)
{
    TrianglesMeshReader<DIM,DIM> mesh_reader(fileName, 2); // 2=quadratic mesh

    ConstructFromMeshReader(mesh_reader);

    // set up the information on whether a node is an internal node or not (if not,
    // it'll be a vertex)
    mIsInternalNode.resize(this->GetNumNodes(), true);
    for (unsigned elem_index=0; elem_index<this->GetNumElements(); elem_index++)
    {
        for (unsigned i=0; i<DIM+1 /*num vertices*/; i++)
        {
            unsigned node_index = this->GetElement(elem_index)->GetNodeGlobalIndex(i);
            mIsInternalNode[ node_index ] = false;
        }
    }
    
    // count the number of vertices, and also check all vertices come before the 
    // rest of the nodes (as this is assumed in other parts of the code)
    mNumVertices = 0;
    bool vertices_mode = true;
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        if (mIsInternalNode[i]==false)
        {
            mNumVertices++;
        }
        if ((vertices_mode == false)  && (mIsInternalNode[i]==false ) )
        {
            EXCEPTION("The quadratic mesh doesn't appear to have all vertices before the rest of the nodes");
        }
        if ( (vertices_mode == true)  && (mIsInternalNode[i]==true) )
        {
            vertices_mode = false;
        }
    }
        
    
    mesh_reader.Reset();

    // add the extra nodes (1 extra node in 1D, 3 in 2D, 6 in 3D) to the element
    // data.
    for (unsigned i=0; i<this->GetNumElements(); i++)
    {
        std::vector<unsigned> nodes = mesh_reader.GetNextElementData().NodeIndices;
        assert(nodes.size()==(DIM+1)*(DIM+2)/2);
        for (unsigned j=DIM+1; j<(DIM+1)*(DIM+2)/2; j++)
        {
            this->GetElement(i)->AddNode( this->GetNode(nodes[j]) );
            this->GetNode(nodes[j])->AddElement(this->GetElement(i)->GetIndex());
        }
    }
    
    // Loop over all boundary elements, find the equivalent face from all
    // the elements, and add the extra nodes to the boundary element
    if (DIM>1)
    {
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
              = this->GetBoundaryElementIteratorBegin();
            iter != this->GetBoundaryElementIteratorEnd();
            ++iter)
        {
            // collect the nodes of this boundary element in a set        
            std::set<unsigned> boundary_element_node_indices;
            for (unsigned i=0; i<DIM; i++)
            {
                boundary_element_node_indices.insert( (*iter)->GetNodeGlobalIndex(i) );
            }
    
            bool found_this_boundary_element = false;

            // loop over elements
            for (unsigned i=0; i<this->GetNumElements(); i++)
            {
                Element<DIM,DIM>* p_element = this->GetElement(i);
                
                // for each element, loop over faces (the opposites to a node)
                for (unsigned face=0; face<DIM+1; face++)
                {
                    // collect the node indices for this face
                    std::set<unsigned> node_indices;
                    for (unsigned local_node_index=0; local_node_index<DIM+1; local_node_index++)
                    {  
                        if (local_node_index!=face)
                        {
                            node_indices.insert( p_element->GetNodeGlobalIndex(local_node_index) );
                        }
                    }
    
                    assert(node_indices.size()==DIM);

                    // see if this face matches the boundary element,
                    // and call AddExtraBoundaryNodes() if so
                    if (node_indices==boundary_element_node_indices)
                    {
                        AddExtraBoundaryNodes(*iter, p_element, face);
                        
                        found_this_boundary_element = true;
                        break;
                    }
                }
    
                if (found_this_boundary_element)
                {
                    break;
                }
            }
            
            if (!found_this_boundary_element)
            {
                #define COVERAGE_IGNORE
                EXCEPTION("Unable to find a face of an element which matches one of the boundary elements");
                #undef COVERAGE_IGNORE
            }
        }
    }
}


template<unsigned DIM>
void QuadraticMesh<DIM>::AddNodeToBoundaryElement(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                                  Element<DIM,DIM>* pElement,
                                                  unsigned internalNode)
{
    assert(DIM>1);
    assert(internalNode >= DIM+1);
    assert(internalNode < (DIM+1)*(DIM+2)/2);
    Node<DIM>* p_internal_node = pElement->GetNode(internalNode);

    // add node to the boundary node list   
    if (!p_internal_node->IsBoundaryNode())
    {
        p_internal_node->SetAsBoundaryNode();
        this->mBoundaryNodes.push_back(p_internal_node);
    }

    pBoundaryElement->AddNode( p_internal_node );        
}


template<unsigned DIM>
void QuadraticMesh<DIM>::AddExtraBoundaryNodes(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                               Element<DIM,DIM>* pElement,
                                               unsigned nodeIndexOppositeToFace)
{
    assert(DIM!=1);
    if (DIM==2)
    {
        assert(nodeIndexOppositeToFace<3);
        // the single internal node of the elements face will be numbered 'face+3'
        AddNodeToBoundaryElement(pBoundaryElement, pElement, nodeIndexOppositeToFace+3);
    }        
    else
    {
        assert(DIM==3);        

        unsigned b_elem_n0 = pBoundaryElement->GetNodeGlobalIndex(0);
        unsigned b_elem_n1 = pBoundaryElement->GetNodeGlobalIndex(1);

        unsigned offset;
        bool reverse;

        if (nodeIndexOppositeToFace==0)
        {
            // face opposite to node 0 = {1,2,3}, with corresponding internals {9,8,5}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 1, 2, 3, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 9, 8, 5, offset, reverse);
        }
        else if (nodeIndexOppositeToFace==1)
        {
            // face opposite to node 1 = {2,0,3}, with corresponding internals {7,9,6}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 2, 0, 3, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 7, 9, 6, offset, reverse);
        }
        else if (nodeIndexOppositeToFace==2)
        {
            // face opposite to node 2 = {0,1,3}, with corresponding internals {8,7,4}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 0, 1, 3, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 8, 7, 4, offset, reverse);
        }
        else
        {
            assert(nodeIndexOppositeToFace==3);
            // face opposite to node 3 = {0,1,2}, with corresponding internals {5,6,4}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 0, 1, 2, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 5, 6, 4, offset, reverse);
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// two unpleasant helper methods for AddExtraBoundaryNodes()
///////////////////////////////////////////////////////////////////////////////

/** This methods takes in the three vertices of a face which match the given boundary
 *  element, and figure out if the order of the nodes in the face is reversed in
 *  the boundary element (returned in the bool 'rReverse'). Also, the offset between
 *  the first node in the face (as given to this method) and the first node in
 *  the boundary element is computed (returned in the variable 'rOffset'). Offset
 *  should then be applied before reverse to match the face nodes to the boundary
 *  element nodes.
 */

#define COVERAGE_IGNORE /// \todo These helper methods aren't properly covered
template<unsigned DIM>
void QuadraticMesh<DIM>::HelperMethod1(unsigned boundaryElemNode0, unsigned boundaryElemNode1,
                                       Element<DIM,DIM>* pElement,
                                       unsigned node0, unsigned node1, unsigned node2,
                                       unsigned& rOffset,
                                       bool& rReverse)
{
    if (pElement->GetNodeGlobalIndex(node0)==boundaryElemNode0)
    {
        rOffset = 0;
        if (pElement->GetNodeGlobalIndex(node1)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else
        {
            assert(pElement->GetNodeGlobalIndex(node2)==boundaryElemNode1);
            rReverse = true;
        }
    }
    else if (pElement->GetNodeGlobalIndex(node1)==boundaryElemNode0)
    {
        rOffset = 1;
        if (pElement->GetNodeGlobalIndex(node2)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else 
        {
            assert(pElement->GetNodeGlobalIndex(node0)==boundaryElemNode1);
            rReverse = true;
        }
    }
    else
    {
        assert(pElement->GetNodeGlobalIndex(node2)==boundaryElemNode0);
        rOffset = 2;
        if (pElement->GetNodeGlobalIndex(node0)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else
        {
            assert(pElement->GetNodeGlobalIndex(node1)==boundaryElemNode1);
            rReverse = true;
        }
    }
}
#undef COVERAGE_IGNORE /// \todo These helper methods aren't properly covered



/**
 *  This method takes the three internal nodes for some face in some element,
 *  applies the given offset and reverse (see HelperMethod1) to them, to get 
 *  the ordered internal nodes which should given to the boundary element.
 *  It then calls AddNodeToBoundaryElement with each of the three internal nodes.
 */
#define COVERAGE_IGNORE /// \todo These helper methods aren't properly covered
template<unsigned DIM>
void QuadraticMesh<DIM>::HelperMethod2(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                       Element<DIM,DIM>* pElement,
                                       unsigned internalNode0, unsigned internalNode1, unsigned internalNode2,
                                       unsigned offset,
                                       bool reverse)
{
    if (offset==1)
    {
        unsigned temp = internalNode0;
        internalNode0 = internalNode1;
        internalNode1 = internalNode2;
        internalNode2 = temp;
    }
    else if (offset == 2)
    {
        unsigned temp = internalNode0;
        internalNode0 = internalNode2;
        internalNode2 = internalNode1;
        internalNode1 = temp;
    }
    
    if (reverse)
    {
        unsigned temp = internalNode1;
        internalNode1 = internalNode2;
        internalNode2 = temp;
    }
    
    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode0);
    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode1);
    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode2);
}
#undef COVERAGE_IGNORE /// \todo These helper methods aren't properly covered

template class QuadraticMesh<1>;
template class QuadraticMesh<2>;
template class QuadraticMesh<3>;
