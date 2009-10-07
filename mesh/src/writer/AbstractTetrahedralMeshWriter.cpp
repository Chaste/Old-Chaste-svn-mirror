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

#include "AbstractTetrahedralMeshWriter.hpp"
#include "AbstractTetrahedralMesh.hpp"

#include "ParallelTetrahedralMesh.hpp"
#include "Debug.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralMeshWriter(const std::string &rDirectory,
                   const std::string &rBaseName,
                   const bool clearOutputDir)
    : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir)
{
}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(
     const AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    //Have we got a parallel mesh?
    const ParallelTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh = dynamic_cast<const ParallelTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);
    if (p_mesh != NULL)
    {
        //It's a parallel mesh
        WriteFilesUsingParallelMesh(*p_mesh);
        return;
    }
    NodeMap node_map(rMesh.GetNumAllNodes());
    unsigned new_index = 0;
    for (unsigned i=0; i<(unsigned)rMesh.GetNumAllNodes(); i++)
    {
        Node<SPACE_DIM>* p_node = rMesh.GetNode(i);

        if (p_node->IsDeleted() == false)
        {
            std::vector<double> coords(SPACE_DIM);
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                coords[j] = p_node->GetPoint()[j];
            }
            this->SetNextNode(coords);
            node_map.SetNewIndex(i, new_index++);
        }
        else
        {
            node_map.SetDeleted(i);
        }
    }
    assert(new_index==(unsigned)rMesh.GetNumNodes());

    for (unsigned i=0; i<(unsigned)rMesh.GetNumAllElements(); i++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = rMesh.GetElement(i);

        if (p_element->IsDeleted() == false)
        {
            std::vector<unsigned> indices(p_element->GetNumNodes());

            for (unsigned j=0; j<indices.size(); j++)
            {
                unsigned old_index = p_element->GetNodeGlobalIndex(j);
                indices[j] = node_map.GetNewIndex(old_index);
            }
            this->SetNextElement(indices);
        }
    }

    for (unsigned i=0; i<(unsigned)rMesh.GetNumAllBoundaryElements(); i++)
    {
        BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = rMesh.GetBoundaryElement(i);
        if (p_boundary_element->IsDeleted() == false)
        {
            std::vector<unsigned> indices(ELEMENT_DIM);
            for (unsigned j=0; j<ELEMENT_DIM; j++)
            {
                unsigned old_index = p_boundary_element->GetNodeGlobalIndex(j);
                indices[j] = node_map.GetNewIndex(old_index);
            }
            this->SetNextBoundaryFace(indices);
        }
    }
    this->WriteFiles();
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingParallelMesh(
     const ParallelTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    //Concentrate node information to the master
    
    MPI_Status status;
    double raw_coords[SPACE_DIM];
    for (unsigned i=0; i<(unsigned)rMesh.GetNumNodes(); i++)
    {
        try 
        {
            Node<SPACE_DIM>* p_node = rMesh.GetNode(i);
            assert (p_node->IsDeleted() == false);
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                raw_coords[j] = p_node->GetPoint()[j];
            }
            MPI_Send(raw_coords, SPACE_DIM, MPI_DOUBLE, 0, i, PETSC_COMM_WORLD);
        }
        catch (Exception e)
        {
        }
        if (PetscTools::AmMaster())
        {
            MPI_Recv(raw_coords, SPACE_DIM, MPI_DOUBLE, MPI_ANY_SOURCE, i, PETSC_COMM_WORLD, &status);
            std::vector<double> coords(SPACE_DIM);
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                coords[j] = raw_coords[j];
            }
            this->SetNextNode(coords);
        }
        else
        {
            assert(this->GetNumNodes() == 0);
        }            
    }
    
    TRACE("Only node file can be written at this stage");
    //Master writes as usual
    if (PetscTools::AmMaster())
    {
        this->WriteFiles();
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractTetrahedralMeshWriter<1,1>;
template class AbstractTetrahedralMeshWriter<1,2>;
template class AbstractTetrahedralMeshWriter<1,3>;
template class AbstractTetrahedralMeshWriter<2,2>;
template class AbstractTetrahedralMeshWriter<2,3>;
template class AbstractTetrahedralMeshWriter<3,3>;
