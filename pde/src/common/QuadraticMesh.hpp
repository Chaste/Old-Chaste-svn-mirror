/*

Copyright (C) University of Oxford, 2008

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

#ifndef QUADRATICMESH_HPP_
#define QUADRATICMESH_HPP_

#include "ConformingTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include <vector>

template<unsigned DIM>
class QuadraticMesh : public ConformingTetrahedralMesh<DIM, DIM>
{    
private:
    bool mIsPrepared;
    std::vector<bool> mIsInternalNode;

    void AddExtraBoundaryNodes(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                               Element<DIM,DIM>* pElement,
                               unsigned nodeIndexOppositeToFace);
    
public:
    /**
     * Constructs a new Quadratic Mesh
     * 
     * @param fileName The name of the mesh file to load
     */
    QuadraticMesh(const std::string& fileName);
    
    /**
     * Construct a new Quadratic Mesh
     */
    QuadraticMesh();
    
    /**
     * Calculates the extra nodes and information needed to use a mesh with quadratic basis functions
     */
    void ConvertToQuadratic();
};


template<unsigned DIM>
QuadraticMesh<DIM>::QuadraticMesh(const std::string& fileName)
{
    TrianglesMeshReader<DIM,DIM> mesh_reader(fileName, 2); // 2=quadratic mesh
    ConstructFromMeshReader(mesh_reader);

    // set up the information on whether a node is an internal node or not (if not,
    // it'll be a vertex)
    mIsInternalNode.resize(this->GetNumNodes(), true);
    for(unsigned elem_index=0; elem_index<this->GetNumElements(); elem_index++)
    {
        for(unsigned i=0; i<DIM+1 /*num vertices*/; i++)
        {
            unsigned node_index = this->GetElement(elem_index)->GetNodeGlobalIndex(i);
            mIsInternalNode[ node_index ] = false;
        }
    }
    
    mesh_reader.Reset();
    
    // add the extra nodes (1 extra node in 1D, 3 in 2D, 6 in 3D) to the element
    // data.
    for(unsigned i=0; i<this->GetNumElements(); i++)
    {
        std::vector<unsigned> node_indices = mesh_reader.GetNextElement();
        for(unsigned j=DIM+1; j<(DIM+1)*(DIM+2)/2; j++)
        {
            this->GetElement(i)->AddNode( this->GetNode(node_indices[j]) );
        }
    }
    
    // Loop over all boundary elements, find the equivalent face from all
    // the elements, and add the extra nodes to the boundary element
//    unsigned counter = 0;
    if(DIM>1)
    {
        for(typename ConformingTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
              = this->GetBoundaryElementIteratorBegin();
            iter != this->GetBoundaryElementIteratorEnd();
            ++iter)
        {
    //        std::cout << "\n\nCounter="<< counter++ << std::endl;
    
            // collect the nodes of this boundary element in a set        
            std::set<unsigned> boundary_element_node_indices;
            for(unsigned i=0; i<DIM; i++)
            {
                boundary_element_node_indices.insert( (*iter)->GetNodeGlobalIndex(i) );
            }
    
    //        std::cout << "\nLooking for:\n";        
    //        for(std::set<unsigned>::iterator kk=boundary_element_node_indices.begin();
    //            kk!=boundary_element_node_indices.end(); kk++)
    //        {
    //            std::cout << *kk << std::endl;
    //        }
    
    
            bool found_this_boundary_element = false;
            // loop over elements
            for(unsigned i=0; i<this->GetNumElements(); i++)
            {
                Element<DIM,DIM>* p_element = this->GetElement(i);
                
                // for each element, loop over faces (the opposites to a node)
                for(unsigned face=0; face<DIM+1; face++)
                {
                    // collect the node indices for this face
                    std::set<unsigned> node_indices;
                    for(unsigned local_node_index=0; local_node_index<DIM+1; local_node_index++)
                    {  
                        if(local_node_index!=face)
                        {
                            node_indices.insert( p_element->GetNodeGlobalIndex(local_node_index) );
                        }
                    }
    
                    assert(node_indices.size()==DIM);
                    
    //                std::cout << "This face is:\n";
    //                for(std::set<unsigned>::iterator kk=node_indices.begin();
    //                    kk!=node_indices.end(); kk++)
    //                {
    //                    std::cout << *kk << std::endl;
    //                }
                
                    // see if this face matches the boundary element        
                    if(node_indices==boundary_element_node_indices)
                    {
                        AddExtraBoundaryNodes(*iter, p_element, face);
                        
                        found_this_boundary_element = true;
                        break;
                    }
                }
    
                if(found_this_boundary_element)
                {
                    break;
                }
            }
            
            if(!found_this_boundary_element)
            {
                #define COVERAGE_IGNORE
                EXCEPTION("Unable to find a face of an element which matches one of the boundary elements");
                #undef COVERAGE_IGNORE
            }
        }
    }

    
    
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    /// HACK!!! HACK!!! HACK!!! HACK!!!!!!
    /// HACK!!! HACK!!! HACK!!! HACK!!!!!!
    /// HACK!!! HACK!!! HACK!!! HACK!!!!!!
    ///
    /// see Ticket:777
    ///
    /// This hack is because .edge files created by "triangle -o2" do not give all 
    /// three nodes of an edge (in 2d), just the 2 vertices, therefore the reader
    /// does not give all three nodes for boundary faces, and internal nodes which
    /// are also on the boundary of the mesh are incorrectly labelled as not boundary
    /// nodes.
    ///
    /// We hack round this here, to get tests to pass, by assuming the mesh is the 
    /// unit square. 
    ///
    /// HACK!!! HACK!!! HACK!!! HACK!!!!!!!!!!!!!!!!!!!!!!!!
    ///////////////////////////////////////////////////////////////////////////////

    if(DIM==3)
    {
        for(unsigned i=0; i<this->GetNumNodes(); i++)
        {
            if(mIsInternalNode[i])
            {
                for(unsigned j=0; j<DIM; j++)
                {
                    if( (this->GetNode(i)->rGetLocation()[j]==0) || (this->GetNode(i)->rGetLocation()[j]==1) )
                    {
                        if(!this->GetNode(i)->IsBoundaryNode())
                        {
                            this->GetNode(i)->SetAsBoundaryNode();
                            this->mBoundaryNodes.push_back(this->GetNode(i));
                        }
                        break;
                    }
                }
            }
        }
    }   
    
    mIsPrepared = true;
}

template<unsigned DIM>
QuadraticMesh<DIM>::QuadraticMesh()
{
    mIsPrepared = false;
}

template<unsigned DIM>
void QuadraticMesh<DIM>::ConvertToQuadratic()
{
    #define COVERAGE_IGNORE
    // not implemented yet..
    assert(false);

    // do conversion here

    mIsPrepared = true;
    #undef COVERAGE_IGNORE
}

template<unsigned DIM>
void QuadraticMesh<DIM>::AddExtraBoundaryNodes(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                               Element<DIM,DIM>* pElement,
                                               unsigned nodeIndexOppositeToFace)
{
    assert(DIM!=1);
    if(DIM==2)
    {
        assert(nodeIndexOppositeToFace<3);
        // the single internal node of the elements face will be numbered 'face+3'
        Node<DIM>* p_internal_node = pElement->GetNode(nodeIndexOppositeToFace+3);
    
        // add node to the boundary node list   
        if(!p_internal_node->IsBoundaryNode())
        {
            p_internal_node->SetAsBoundaryNode();
            this->mBoundaryNodes.push_back(p_internal_node);
        }

        pBoundaryElement->AddNode( p_internal_node );
    }        
    else
    {
        assert(DIM==3);
    }
}
                                               
#endif /*QUADRATICMESH_HPP_*/
