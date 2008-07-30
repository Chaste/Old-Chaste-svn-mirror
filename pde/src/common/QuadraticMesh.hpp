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

// currently just 2d
class QuadraticMesh : public ConformingTetrahedralMesh<2,2>
{
private:
    bool mIsPrepared;
    std::vector<std::vector<unsigned> > mLnods;
    
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
    
    unsigned GetElementNode(unsigned elemIndex, unsigned nodeIndex)
    {
        assert(mIsPrepared);
        assert(elemIndex<this->GetNumElements());
        assert(nodeIndex>=3 && nodeIndex<6);
        return mLnods[elemIndex][(unsigned)(nodeIndex-3)];
    }
};



QuadraticMesh::QuadraticMesh(const std::string& fileName)
{
    TrianglesMeshReader<2,2> mesh_reader(fileName, 2);
    ConstructFromMeshReader(mesh_reader);
    
    mesh_reader.Reset();
    
    mLnods.resize(this->GetNumElements());
    for(unsigned i=0; i<mLnods.size(); i++)
    {
        std::vector<unsigned> node_indices = mesh_reader.GetNextElement();
        
        mLnods[i].push_back( node_indices[3] );
        mLnods[i].push_back( node_indices[4] );
        mLnods[i].push_back( node_indices[5] );
    }
    
    mIsPrepared = true;
}


QuadraticMesh::QuadraticMesh()
{
    mIsPrepared = false;
}

void QuadraticMesh::ConvertToQuadratic()
{
    #define COVERAGE_IGNORE
    // not implemented yet..
    assert(false);

    // do conversion here

    mIsPrepared = true;
    #undef COVERAGE_IGNORE

}

#endif /*QUADRATICMESH_HPP_*/
