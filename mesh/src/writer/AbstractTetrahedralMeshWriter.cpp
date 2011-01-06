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

// Disable PETSc logging of MPI calls (we don't use this, anyway) to fix
// "right-hand operand of comma has no effect" warnings when building with
// PETSc 2.2.1.
#define PETSC_HAVE_BROKEN_RECURSIVE_MACRO

#include "AbstractTetrahedralMeshWriter.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "Version.hpp"

#include <mpi.h> // For MPI_Send, MPI_Recv

/**
 * Convenience collection of iterators, primarily to get compilation to happen.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MeshWriterIterators
{
    /** Iterator over nodes */
    typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator* pNodeIter;
    /** Iterator over elements */
    typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator* pElemIter;
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralMeshWriter(const std::string &rDirectory,
                   const std::string &rBaseName,
                   const bool clearOutputDir)
    : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
      mpMesh(NULL),
      mpNodeMap(NULL),
      mNodesPerElement(ELEMENT_DIM+1),
      mpDistributedMesh(NULL),
      mpIters(new MeshWriterIterators<ELEMENT_DIM,SPACE_DIM>),
      mNodeCounterForParallelMesh(0),
      mElementCounterForParallelMesh(0),
      mFilesAreBinary(false)
{
    mpIters->pNodeIter = NULL;
    mpIters->pElemIter = NULL;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::~AbstractTetrahedralMeshWriter()
{
    if (mpIters->pNodeIter)
    {
        delete mpIters->pNodeIter;
    }
    if (mpIters->pElemIter)
    {
        delete mpIters->pElemIter;
    }

    delete mpIters;

    if (mpNodeMap)
    {
        delete mpNodeMap;
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    // if we are writing from a mesh..
    assert(PetscTools::AmMaster());
    if (mpMesh)
    {
        std::vector<double> coords(SPACE_DIM);
        double raw_coords[SPACE_DIM];

        //Iterate over the locally-owned nodes
        if ( (*(mpIters->pNodeIter)) != mpMesh->GetNodeIteratorEnd())
        {
            //Either this a sequential mesh (and we own it all)
            // or it's parallel (and the master owns the first chunk)
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                coords[j] = (*(mpIters->pNodeIter))->GetPoint()[j];
            }

            mNodeCounterForParallelMesh=(*(mpIters->pNodeIter))->GetIndex() + 1;//Ready for when we run out of local nodes

            ++(*(mpIters->pNodeIter));
            return coords;
        }

        //If we didn't return then the iterator has reached the end of the local nodes.
        // It must be a parallel mesh and we are expecting messages...

        assert( mpDistributedMesh != NULL );

        MPI_Status status;
        // do receive, convert to std::vector on master
        MPI_Recv(raw_coords, SPACE_DIM, MPI_DOUBLE, MPI_ANY_SOURCE, mNodeCounterForParallelMesh, PETSC_COMM_WORLD, &status);
        assert(status.MPI_ERROR == MPI_SUCCESS);
        for (unsigned j=0; j<coords.size(); j++)
        {
            coords[j] = raw_coords[j];
        }


        mNodeCounterForParallelMesh++;
        return coords;
    }
    else
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextNode();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElement()
{
    assert(PetscTools::AmMaster());
    // if we are writing from a mesh..
    if (mpMesh)
    {
        ElementData elem_data;
        elem_data.NodeIndices.resize(mNodesPerElement);

        if ( mpDistributedMesh == NULL ) // not using parallel mesh
        {
            // Use the iterator
            assert(this->mNumElements==mpMesh->GetNumElements());

            for (unsigned j=0; j<elem_data.NodeIndices.size(); j++)
            {
                unsigned old_index = (*(mpIters->pElemIter))->GetNodeGlobalIndex(j);
                elem_data.NodeIndices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
            }
            // Set attribute
            elem_data.AttributeValue = (*(mpIters->pElemIter))->GetRegion();

            ++(*(mpIters->pElemIter));

            return elem_data;
        }
        else //Parallel mesh
        {
            //Use the mElementCounterForParallelMesh variable to identify next element
            if ( mpDistributedMesh->CalculateDesignatedOwnershipOfElement( mElementCounterForParallelMesh ) == true )
            {
                //Master owns this element
                Element<ELEMENT_DIM, SPACE_DIM>* p_element = mpDistributedMesh->GetElement(mElementCounterForParallelMesh);
                assert(elem_data.NodeIndices.size() == ELEMENT_DIM+1);
                assert( ! p_element->IsDeleted() );
                //Master can use the local data to recover node indices & attribute
                for (unsigned j=0; j<ELEMENT_DIM+1; j++)
                {
                    elem_data.NodeIndices[j] = p_element->GetNodeGlobalIndex(j);
                }
                elem_data.AttributeValue = p_element->GetRegion();
            }
            else
            {
                //Master doesn't own this element.
                // +2 to allow for attribute value too
                unsigned raw_indices[ELEMENT_DIM+2];
                MPI_Status status;
                //Get it from elsewhere
                MPI_Recv(raw_indices, ELEMENT_DIM+2, MPI_UNSIGNED, MPI_ANY_SOURCE,
                         this->mNumNodes + mElementCounterForParallelMesh,
                         PETSC_COMM_WORLD, &status);
                // Convert to std::vector
                for (unsigned j=0; j< elem_data.NodeIndices.size(); j++)
                {
                    elem_data.NodeIndices[j] = raw_indices[j];
                }
                // Attribute value
                elem_data.AttributeValue = raw_indices[ELEMENT_DIM+1];
            }
            // increment element counter
            mElementCounterForParallelMesh++;

            return elem_data; // non-master processors will return garbage here - but they should never write to file
        }
    }
    else // not writing from a mesh
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextElement();
    }
}


///\todo #1322 Mesh should be const
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(
      AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
      bool keepOriginalElementIndexing)
{
    this->mpMeshReader = NULL;
    mpMesh = &rMesh;

    this->mNumNodes = mpMesh->GetNumNodes();
    this->mNumElements = mpMesh->GetNumElements();
    this->mNumBoundaryElements = mpMesh->GetNumBoundaryElements();

    typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
    mpIters->pNodeIter = new NodeIterType(mpMesh->GetNodeIteratorBegin());

    typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator ElemIterType;
    mpIters->pElemIter = new ElemIterType(mpMesh->GetElementIteratorBegin());

    //Use this processes first element to gauge the size of all the elements
    if ( (*(mpIters->pElemIter)) != mpMesh->GetElementIteratorEnd())
    {
        mNodesPerElement = (*(mpIters->pElemIter))->GetNumNodes();
    }
    //Connectivity file is written when we write to a binary file (only available for TrianglesMeshWriter) and if we are preserving the element order
    if (this->mFilesAreBinary && keepOriginalElementIndexing)
    {
        unsigned max_elements_all;
        if (PetscTools::IsSequential())
        {
            max_elements_all = mpMesh->CalculateMaximumContainingElementsPerProcess();
        }
        else
        {
            unsigned max_elements_per_process = mpMesh->CalculateMaximumContainingElementsPerProcess();
            MPI_Allreduce(&max_elements_per_process, &max_elements_all, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD);
        }


        for (unsigned writing_process=0; writing_process<PetscTools::GetNumProcs(); writing_process++)
        {
            if (PetscTools::GetMyRank() == writing_process)
            {
                std::string node_connect_list_file_name = this->mBaseName + ".ncl";
                out_stream p_ncl_file=out_stream(NULL);
                
                if (PetscTools::AmMaster())
                {
                    assert(writing_process==0);
                    //Open the file for the first time                    
                    p_ncl_file = this->mpOutputFileHandler->OpenOutputFile(node_connect_list_file_name);
            
                    // Write the ncl header            
                    *p_ncl_file << this->mNumNodes << "\t";
                    *p_ncl_file << max_elements_all << "\t";
                    *p_ncl_file << "\tBIN\n";
                }
                else
                {
                    //Append to the existing file
                    p_ncl_file = this->mpOutputFileHandler->OpenOutputFile(node_connect_list_file_name, std::ios::app);                    
                }
                
                // Write each node's data
                unsigned default_marker = UINT_MAX;
    
                for (NodeIterType iter = mpMesh->GetNodeIteratorBegin();
                     iter != mpMesh->GetNodeIteratorEnd();
                     ++iter)
                {
                    //Get the containing element indices from the node's set and sort them
                    std::set<unsigned>& r_elem_set = iter->rGetContainingElementIndices();
                    std::vector<unsigned> elem_vector(r_elem_set.begin(),r_elem_set.end()); 
                    std::sort(elem_vector.begin(), elem_vector.end());
                    //Pad the vector with unsigned markers
                    for (unsigned elem_index=elem_vector.size();  elem_index<max_elements_all; elem_index++)
                    {
                        elem_vector.push_back(default_marker);
                    }
                    assert (elem_vector.size() == max_elements_all);
                    //Write raw data out of std::vector into the file
                    p_ncl_file->write((char*)&elem_vector[0], elem_vector.size()*sizeof(unsigned));
                }
                
                if (PetscTools::AmTopMost())
                {
                    *p_ncl_file << "#\n# " + ChasteBuildInfo::GetProvenanceString();
                }
                
                p_ncl_file->close();
            }
            
            PetscTools::Barrier("AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh");
        }
    }    


    //Have we got a parallel mesh?
    ///\todo #1322 This should be const too
    mpDistributedMesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);

    if (mpDistributedMesh != NULL)
    {
        //It's a parallel mesh
        WriteFilesUsingParallelMesh(keepOriginalElementIndexing);
        return;
    }

    if (!PetscTools::AmMaster())
    {
        PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingMesh"); //Paired with Master process writing files
        return;
    }

    // Set up node map if we might have deleted nodes
    unsigned node_map_index = 0;
    if (mpMesh->IsMeshChanging())
    {
        mpNodeMap = new NodeMap(mpMesh->GetNumAllNodes());
        for (NodeIterType it = mpMesh->GetNodeIteratorBegin(); it != mpMesh->GetNodeIteratorEnd(); ++it)
        {
            mpNodeMap->SetNewIndex(it->GetIndex(), node_map_index++);
        }
    }

    // Cache all of the BoundaryElements
    for (unsigned i=0; i<(unsigned)mpMesh->GetNumAllBoundaryElements(); i++)
    {
        BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = mpMesh->GetBoundaryElement(i);
        if (p_boundary_element->IsDeleted() == false)
        {
            std::vector<unsigned> indices(p_boundary_element->GetNumNodes());
            for (unsigned j=0; j<p_boundary_element->GetNumNodes(); j++)
            {
                unsigned old_index = p_boundary_element->GetNodeGlobalIndex(j);
                indices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
            }
            this->SetNextBoundaryFace(indices);
        }
    }
    this->WriteFiles();
    PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingMesh"); //Paired with waiting Slave processes
    delete mpIters->pNodeIter;
    mpIters->pNodeIter=NULL;
    delete mpIters->pElemIter;
    mpIters->pElemIter=NULL;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingParallelMesh(bool keepOriginalElementIndexing)
{
    if (keepOriginalElementIndexing)
    {
        //Concentrate node information to the master
    
        MPI_Status status;
        double raw_coords[SPACE_DIM];
        //Concentrate and cache all of the BoundaryElements
        unsigned raw_face_indices[ELEMENT_DIM];//Assuming that we don't have parallel quadratic meshes
        for (unsigned index=0; index<(unsigned)mpDistributedMesh->GetNumBoundaryElements(); index++)
        {
            try
            {
                if ( mpDistributedMesh->CalculateDesignatedOwnershipOfBoundaryElement( index ) == true )
                {
                    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = mpDistributedMesh->GetBoundaryElement(index);
                    assert(p_boundary_element->IsDeleted() == false);
                    for (unsigned j=0; j<ELEMENT_DIM; j++)
                    {
                        raw_face_indices[j] = p_boundary_element->GetNodeGlobalIndex(j);
                    }
                    MPI_Send(raw_face_indices, ELEMENT_DIM, MPI_UNSIGNED, 0, index, PETSC_COMM_WORLD);
                }
            }
            catch (Exception e)
            {
            }
            if (PetscTools::AmMaster())
            {
                MPI_Recv(raw_face_indices, ELEMENT_DIM, MPI_UNSIGNED, MPI_ANY_SOURCE, index, PETSC_COMM_WORLD, &status);
                std::vector<unsigned> indices(ELEMENT_DIM);
                for (unsigned j=0; j<indices.size(); j++)
                {
                    indices[j] = raw_face_indices[j];
                }
                this->SetNextBoundaryFace(indices);
            }
        }
        PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingParallelMesh 1");
    
        //Master goes on to write as usual
        if (PetscTools::AmMaster())
        {
            this->WriteFiles();
        }
        else
        {
            //Slaves concentrate the Nodes and Elements
            typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
            for (NodeIterType it = mpMesh->GetNodeIteratorBegin(); it != mpMesh->GetNodeIteratorEnd(); ++it)
            {
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    raw_coords[j] = it->GetPoint()[j];
                }
                MPI_Send(raw_coords, SPACE_DIM, MPI_DOUBLE, 0, it->GetIndex(), PETSC_COMM_WORLD);//Nodes sent with positive tags
            }
    
            // +2 allows for attribute value
            unsigned raw_indices[ELEMENT_DIM+2]; //Assuming that we don't have parallel quadratic elements
            typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator ElementIterType;
            for (ElementIterType it = mpMesh->GetElementIteratorBegin(); it != mpMesh->GetElementIteratorEnd(); ++it)
            {
                unsigned index =it->GetIndex();
                if ( mpDistributedMesh->CalculateDesignatedOwnershipOfElement( index ) == true )
                {
                    for (unsigned j=0; j<ELEMENT_DIM+1; j++)
                    {
                        raw_indices[j] = it->GetNodeGlobalIndex(j);
                    }
                    // Attribute value
                    raw_indices[ELEMENT_DIM+1] = it->GetRegion();
                    MPI_Send(raw_indices, ELEMENT_DIM+2, MPI_UNSIGNED, 0,
                             this->mNumNodes + (it->GetIndex()), //Elements sent with tags offset
                             PETSC_COMM_WORLD);
                }
            }
        }
        PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingParallelMesh 2");
    }
    else
    {
        for (unsigned writing_process=0; writing_process<PetscTools::GetNumProcs(); writing_process++)
        {
            if(PetscTools::GetMyRank() == writing_process)
            {
                if (PetscTools::AmMaster())
                {
                    // Make sure headers are written first
                    assert(PetscTools::GetMyRank() == 0);                   
                    CreateFilesWithHeaders();
                }                

                AppendLocalDataToFiles();
                
                if (PetscTools::AmTopMost())
                {
                    // Make sure footers are written last
                    assert(PetscTools::GetMyRank() == PetscTools::GetNumProcs()-1);
                    WriteFilesFooter();
                }        
            }
            //Process i+1 waits for process i to close the file
            PetscTools::Barrier();
        }//Loop in writing_process
                
    }    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::CreateFilesWithHeaders()
{
    // If control reaches this line you haven't implemented the optimised 
    // parallel write for whichever visualiser you are writing for.
    NEVER_REACHED;
}
   
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::AppendLocalDataToFiles()
{
    // If control reaches this line you haven't implemented the optimised 
    // parallel write for whichever visualiser you are writing for.
    NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesFooter()
{
    // If control reaches this line you haven't implemented the optimised 
    // parallel write for whichever visualiser you are writing for.
    NEVER_REACHED;
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
