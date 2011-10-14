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
#include "MixedDimensionMesh.hpp"
#include "Version.hpp"
#include "Exception.hpp"

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
    /** Iterator over boundary elements */
    typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator* pBoundaryElemIter;
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
      mNodesPerBoundaryElement(ELEMENT_DIM),
      mpDistributedMesh(NULL),
      mpMixedMesh(NULL),
      mpIters(new MeshWriterIterators<ELEMENT_DIM,SPACE_DIM>),
      mNodeCounterForParallelMesh(0),
      mElementCounterForParallelMesh(0),
      mBoundaryElementCounterForParallelMesh(0),
      mCableElementCounterForParallelMesh(0),
      mFilesAreBinary(false)
{
    mpIters->pNodeIter = NULL;
    mpIters->pElemIter = NULL;
    mpIters->pBoundaryElemIter = NULL;
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
    if (mpIters->pBoundaryElemIter)
    {
        delete mpIters->pBoundaryElemIter;
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
            // Either this is a sequential mesh (and we own it all)
            // or it's parallel (and the master owns the first chunk)
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                coords[j] = (*(mpIters->pNodeIter))->GetPoint()[j];
            }

            mNodeCounterForParallelMesh=(*(mpIters->pNodeIter))->GetIndex() + 1;//Ready for when we run out of local nodes

            ++(*(mpIters->pNodeIter));
            return coords;
        }

        // If we didn't return then the iterator has reached the end of the local nodes.
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
                assert(elem_data.NodeIndices.size() == mNodesPerElement);
                assert( ! p_element->IsDeleted() );
                //Master can use the local data to recover node indices & attribute
                for (unsigned j=0; j<mNodesPerElement; j++)
                {
                    elem_data.NodeIndices[j] = p_element->GetNodeGlobalIndex(j);
                }
                elem_data.AttributeValue = p_element->GetRegion();
            }
            else
            {
                //Master doesn't own this element.
                // +1 to allow for attribute value too
                unsigned raw_indices[mNodesPerElement+1];
                MPI_Status status;
                //Get it from elsewhere
                MPI_Recv(raw_indices, mNodesPerElement+1, MPI_UNSIGNED, MPI_ANY_SOURCE,
                         this->mNumNodes + mElementCounterForParallelMesh,
                         PETSC_COMM_WORLD, &status);
                // Convert to std::vector
                for (unsigned j=0; j< elem_data.NodeIndices.size(); j++)
                {
                    elem_data.NodeIndices[j] = raw_indices[j];
                }

                // Attribute value
                elem_data.AttributeValue = raw_indices[mNodesPerElement];
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextBoundaryElement()
{
    assert(PetscTools::AmMaster());
    // if we are writing from a mesh..
    if (mpMesh)
    {
        ElementData boundary_elem_data;
        boundary_elem_data.NodeIndices.resize(mNodesPerBoundaryElement);

        if ( mpDistributedMesh == NULL ) // not using parallel mesh
        {
            // Use the iterator
            assert(this->mNumBoundaryElements==mpMesh->GetNumBoundaryElements());

            for (unsigned j=0; j<boundary_elem_data.NodeIndices.size(); j++)
            {
                unsigned old_index = (*(*(mpIters->pBoundaryElemIter)))->GetNodeGlobalIndex(j);
                boundary_elem_data.NodeIndices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
            }

            ++(*(mpIters->pBoundaryElemIter));

            return boundary_elem_data;
        }
        else //Parallel mesh
        {
            //Use the mElementCounterForParallelMesh variable to identify next element
            if ( mpDistributedMesh->CalculateDesignatedOwnershipOfBoundaryElement( mBoundaryElementCounterForParallelMesh ) == true )
            {
                //Master owns this boundary element
                BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* p_boundary_element = mpDistributedMesh->GetBoundaryElement(mBoundaryElementCounterForParallelMesh);
                assert(boundary_elem_data.NodeIndices.size() == ELEMENT_DIM);
                assert( ! p_boundary_element->IsDeleted() );
                //Master can use the local data to recover node indices & attribute
                for (unsigned j=0; j<ELEMENT_DIM; j++)
                {
                    boundary_elem_data.NodeIndices[j] = p_boundary_element->GetNodeGlobalIndex(j);
                }
            }
            else
            {
                //Master doesn't own this boundary element.
                unsigned raw_indices[ELEMENT_DIM];
                MPI_Status status;
                //Get it from elsewhere
                MPI_Recv(raw_indices, ELEMENT_DIM, MPI_UNSIGNED, MPI_ANY_SOURCE,
                         this->mNumNodes + this->mNumElements + mBoundaryElementCounterForParallelMesh,
                         PETSC_COMM_WORLD, &status);
                // Convert to std::vector
                for (unsigned j=0; j< boundary_elem_data.NodeIndices.size(); j++)
                {
                    boundary_elem_data.NodeIndices[j] = raw_indices[j];
                }
            }
            // increment element counter
            mBoundaryElementCounterForParallelMesh++;

            return boundary_elem_data; // non-master processors will return garbage here - but they should never write to file
        }
    }
    else // not writing from a mesh
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextBoundaryElement();
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextCableElement()
{
    assert(PetscTools::AmMaster());

    // if we are writing from a mesh..
    if (mpMesh)
    {
        // Need to be using a MixedDimensionMesh or there will be no cable data
        assert(mpMixedMesh);

        ElementData elem_data;
        elem_data.NodeIndices.resize(2);

        //Use the mCableElementCounterForParallelMesh variable to identify next element
        if ( mpMixedMesh->CalculateDesignatedOwnershipOfCableElement( mCableElementCounterForParallelMesh ) == true )
        {
            //Master owns this element
            Element<1, SPACE_DIM>* p_element = mpMixedMesh->GetCableElement(mCableElementCounterForParallelMesh);
            assert( ! p_element->IsDeleted() );
            //Master can use the local data to recover node indices & attribute
            for (unsigned j=0; j<2; j++)
            {
                elem_data.NodeIndices[j] = p_element->GetNodeGlobalIndex(j);
            }
            elem_data.AttributeValue = p_element->GetRegion();
        }
        else
        {
            //Master doesn't own this element.
            // size is 3 to allow for 2 indices and attribute value too
            unsigned raw_indices[3];
            MPI_Status status;
            //Get it from elsewhere
            MPI_Recv(raw_indices, 3, MPI_UNSIGNED, MPI_ANY_SOURCE,
                     this->mNumNodes + this->mNumElements + this->mNumBoundaryElements + mCableElementCounterForParallelMesh,
                     PETSC_COMM_WORLD, &status);

            // Convert to ElementData (2 nodes plus an attribute value)
            for (unsigned j=0; j< 2; j++)
            {
                elem_data.NodeIndices[j] = raw_indices[j];
            }
            // Attribute value
            elem_data.AttributeValue = raw_indices[2];
        }
        // increment element counter
        mCableElementCounterForParallelMesh++;

        return elem_data; // non-master processors will return garbage here - but they should never write to file
    }
    else // not writing from a mesh
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextCableElement();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteNclFile(
        AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        bool invertMeshPermutation)
{
    unsigned max_elements_all;
    if (PetscTools::IsSequential())
    {
        max_elements_all = rMesh.CalculateMaximumContainingElementsPerProcess();
    }
    else
    {
        unsigned max_elements_per_process = rMesh.CalculateMaximumContainingElementsPerProcess();
        MPI_Allreduce(&max_elements_per_process, &max_elements_all, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD);
    }

    std::string node_connect_list_file_name = this->mBaseName + ".ncl";
    if (invertMeshPermutation && !rMesh.rGetNodePermutation().empty())
    {
        node_connect_list_file_name += "-tmp";
    }

    PetscTools::BeginRoundRobin();
    {
        out_stream p_ncl_file=out_stream(NULL);

        if (PetscTools::AmMaster())
        {
            //Open the file for the first time
            p_ncl_file = this->mpOutputFileHandler->OpenOutputFile(node_connect_list_file_name);

            // Write the ncl header
            *p_ncl_file << rMesh.GetNumNodes() << "\t";
            *p_ncl_file << max_elements_all << "\t";
            *p_ncl_file << "\tBIN\n";
        }
        else
        {
            // Append to the existing file
            p_ncl_file = this->mpOutputFileHandler->OpenOutputFile(node_connect_list_file_name, std::ios::app);
        }

        // Write each node's data
        unsigned default_marker = UINT_MAX;

        typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
        for (NodeIterType iter = rMesh.GetNodeIteratorBegin();
             iter != rMesh.GetNodeIteratorEnd();
             ++iter)
        {
            // Get the containing element indices from the node's set and sort them
            std::set<unsigned>& r_elem_set = iter->rGetContainingElementIndices();
            std::vector<unsigned> elem_vector(r_elem_set.begin(), r_elem_set.end());
            std::sort(elem_vector.begin(), elem_vector.end());
            // Pad the vector with unsigned markers
            for (unsigned elem_index=elem_vector.size();  elem_index<max_elements_all; elem_index++)
            {
                elem_vector.push_back(default_marker);
            }
            assert (elem_vector.size() == max_elements_all);
            // Write raw data out of std::vector into the file
            p_ncl_file->write((char*)&elem_vector[0], elem_vector.size()*sizeof(unsigned));
        }

        if (PetscTools::AmTopMost())
        {
            *p_ncl_file << "#\n# " + ChasteBuildInfo::GetProvenanceString();
        }

        p_ncl_file->close();
    }
    PetscTools::EndRoundRobin();

    if (invertMeshPermutation && PetscTools::AmMaster() && !rMesh.rGetNodePermutation().empty())
    {
        // Open files
        const std::string real_node_connect_list_file_name = this->mBaseName + ".ncl";
        out_stream p_ncl_file = this->mpOutputFileHandler->OpenOutputFile(real_node_connect_list_file_name);
        FileFinder temp_ncl_path = this->mpOutputFileHandler->FindFile(node_connect_list_file_name);
        std::ifstream temp_ncl_file(temp_ncl_path.GetAbsolutePath().c_str());
        // Copy the header
        std::string header_line;
        getline(temp_ncl_file, header_line);
        (*p_ncl_file) << header_line << "\n";
        const std::streampos data_start = temp_ncl_file.tellg();
        const std::streamoff item_width = max_elements_all * sizeof(unsigned);
        // Copy the binary data, permuted
        std::vector<unsigned> elem_vector(max_elements_all);
        for (unsigned node_index=0; node_index<rMesh.GetNumAllNodes(); node_index++)
        {
            unsigned permuted_index = rMesh.rGetNodePermutation()[node_index];
            temp_ncl_file.seekg(data_start + item_width * permuted_index, std::ios_base::beg);
            temp_ncl_file.read((char*)&elem_vector[0], max_elements_all*sizeof(unsigned));
            p_ncl_file->write((char*)&elem_vector[0], max_elements_all*sizeof(unsigned));
        }
        // Footer
        *p_ncl_file << "#\n# " + ChasteBuildInfo::GetProvenanceString();
        p_ncl_file->close();
        // Remove temp file
        remove(temp_ncl_path.GetAbsolutePath().c_str());
    }
    PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteNclFile");
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
    this->mNumCableElements = mpMesh->GetNumCableElements();

    typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
    mpIters->pNodeIter = new NodeIterType(mpMesh->GetNodeIteratorBegin());

    typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator ElemIterType;
    mpIters->pElemIter = new ElemIterType(mpMesh->GetElementIteratorBegin());

    typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator BoundaryElemIterType;
    mpIters->pBoundaryElemIter = new BoundaryElemIterType(mpMesh->GetBoundaryElementIteratorBegin());

    //Use this process's first element to gauge the size of all the elements
    if ( (*(mpIters->pElemIter)) != mpMesh->GetElementIteratorEnd())
    {
        mNodesPerElement = (*(mpIters->pElemIter))->GetNumNodes();
    }

    //Use this process's first boundary element to gauge the size of all the boundary elements
    if ( (*(mpIters->pBoundaryElemIter)) != mpMesh->GetBoundaryElementIteratorEnd())
    {
        mNodesPerBoundaryElement = (*(*(mpIters->pBoundaryElemIter)))->GetNumNodes();
    }

    //Connectivity file is written when we write to a binary file (only available for TrianglesMeshWriter) and if we are preserving the element order
    if (this->mFilesAreBinary && keepOriginalElementIndexing)
    {
        WriteNclFile(*mpMesh);
    }

    // Have we got a parallel mesh?
    ///\todo #1322 This should be const too
    mpDistributedMesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);

    // Have we got a MixedDimensionMesh?
    ///\todo #1322,  This should be const too
    mpMixedMesh = dynamic_cast<MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* >(this->mpMesh);

    if (mpDistributedMesh != NULL)
    {
        // It's a parallel mesh
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

    this->WriteFiles();
    PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingMesh"); // Paired with waiting Slave processes
    delete mpIters->pNodeIter;
    mpIters->pNodeIter = NULL;
    delete mpIters->pElemIter;
    mpIters->pElemIter = NULL;
    delete mpIters->pBoundaryElemIter;
    mpIters->pBoundaryElemIter = NULL;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReaderAndMesh(
        AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
        AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    WriteNclFile(rMesh, true);
    WriteFilesUsingMeshReader(rMeshReader);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingParallelMesh(bool keepOriginalElementIndexing)
{
    if (keepOriginalElementIndexing)
    {
        // Master goes on to write as usual
        if (PetscTools::AmMaster())
        {
            this->WriteFiles();
        }
        else
        {
            double raw_coords[SPACE_DIM];
            // Slaves concentrate the Nodes
            typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
            for (NodeIterType it = mpMesh->GetNodeIteratorBegin(); it != mpMesh->GetNodeIteratorEnd(); ++it)
            {
                for (unsigned j=0; j<SPACE_DIM; j++)
                {
                    raw_coords[j] = it->GetPoint()[j];
                }
                MPI_Send(raw_coords, SPACE_DIM, MPI_DOUBLE, 0, it->GetIndex(), PETSC_COMM_WORLD);//Nodes sent with positive tags
            }

            // Slaves concentrate the Elements for which they are owners
            // +1 allows for attribute value
            unsigned raw_indices[mNodesPerElement+1]; // Assuming that we don't have parallel quadratic elements
            typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator ElementIterType;
            for (ElementIterType it = mpMesh->GetElementIteratorBegin(); it != mpMesh->GetElementIteratorEnd(); ++it)
            {
                unsigned index =it->GetIndex();
                if ( mpDistributedMesh->CalculateDesignatedOwnershipOfElement( index ) == true )
                {
                    for (unsigned j=0; j<mNodesPerElement; j++)
                    {
                        raw_indices[j] = it->GetNodeGlobalIndex(j);
                    }
                    // Attribute value
                    raw_indices[mNodesPerElement] = it->GetRegion();

                    MPI_Send(raw_indices, mNodesPerElement+1, MPI_UNSIGNED, 0,
                             this->mNumNodes + index, //Elements sent with tags offset
                             PETSC_COMM_WORLD);
                }
            }

            // Slaves concentrate the Faces for which they are owners (not in 1-d)
            if (ELEMENT_DIM != 1)
            {
                unsigned raw_face_indices[ELEMENT_DIM]; // Assuming that we don't have parallel quadratic meshes
                typedef typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator BoundaryElementIterType;
                for (BoundaryElementIterType it = mpMesh->GetBoundaryElementIteratorBegin(); it != mpMesh->GetBoundaryElementIteratorEnd(); ++it)
                {
                    unsigned index =(*it)->GetIndex();
                    if ( mpDistributedMesh->CalculateDesignatedOwnershipOfBoundaryElement( index ) == true )
                    {
                        for (unsigned j=0; j<ELEMENT_DIM; j++)
                        {
                            raw_face_indices[j] = (*it)->GetNodeGlobalIndex(j);
                        }
                        MPI_Send(raw_face_indices, ELEMENT_DIM, MPI_UNSIGNED, 0,
                                 this->mNumNodes + this->mNumElements + index, //Faces sent with tags offset even more
                                 PETSC_COMM_WORLD);
                    }
                }
            }

            // Slaves concentrate the cable elements for which they are owners
            if (mpMixedMesh)
            {
                unsigned raw_cable_element_indices[3];
                typedef typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::CableElementIterator CableElementIterType;
                for (CableElementIterType it = mpMixedMesh->GetCableElementIteratorBegin(); it != mpMixedMesh->GetCableElementIteratorEnd(); ++it)
                {
                    unsigned index =(*it)->GetIndex();
                    if ( mpMixedMesh->CalculateDesignatedOwnershipOfCableElement( index ) == true )
                    {
                        for (unsigned j=0; j<2; j++)
                        {
                            raw_cable_element_indices[j] = (*it)->GetNodeGlobalIndex(j);
                        }
                        raw_cable_element_indices[2] = (*it)->GetRegion();
                        MPI_Send(raw_cable_element_indices, 3, MPI_UNSIGNED, 0,
                                 this->mNumNodes + this->mNumElements + this->mNumBoundaryElements + index, //Cable elements sent with tags offset even more
                                 PETSC_COMM_WORLD);
                    }
                }
            }
        }
        PetscTools::Barrier("AbstractTetrahedralMeshWriter::WriteFilesUsingParallelMesh");
    }
    else
    {
        PetscTools::BeginRoundRobin();

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

        PetscTools::EndRoundRobin();
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
