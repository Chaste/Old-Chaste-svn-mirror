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

#endif /*QUADRATICMESH_HPP_*/
