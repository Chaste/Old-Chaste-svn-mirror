/*

Copyright (C) University of Oxford, 2005-2010

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

#include "DistributedTetrahedralMesh.hpp"

#include <cassert>
#include <sstream>
#include <string>

#include "Exception.hpp"
#include "Element.hpp"
#include "BoundaryElement.hpp"

#include "PetscTools.hpp"
#include "DistributedVectorFactory.hpp"
#include "OutputFileHandler.hpp"

/////////////////////////////////////////////////////////////////////////////////////
//   IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DistributedTetrahedralMesh(PartitionType metisPartitioning)
    : mMetisPartitioning((SPACE_DIM!=1)?metisPartitioning:DUMB)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~DistributedTetrahedralMesh()
{
    for (unsigned i=0; i<this->mHaloNodes.size(); i++)
    {
        delete this->mHaloNodes[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ComputeMeshPartitioning(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
    std::set<unsigned>& rNodesOwned,
    std::set<unsigned>& rHaloNodesOwned,
    std::set<unsigned>& rElementsOwned,
    std::vector<unsigned>& rProcessorsOffset)
{
    ///\todo: add a timing event for the partitioning

    if (mMetisPartitioning==METIS_BINARY && !PetscTools::IsSequential())
    {
        MetisBinaryNodePartitioning(rMeshReader, rNodesOwned, rProcessorsOffset);
    }
    else if (mMetisPartitioning==METIS_LIBRARY && !PetscTools::IsSequential())
    {
        MetisLibraryNodePartitioning(rMeshReader, rNodesOwned, rProcessorsOffset);
    }
    else
    {
        DumbNodePartitioning(rMeshReader, rNodesOwned);
    }

    for (unsigned element_number = 0; element_number < mTotalNumElements; element_number++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();

        bool element_owned = false;
        std::set<unsigned> temp_halo_nodes;

        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            if (rNodesOwned.find(element_data.NodeIndices[i]) != rNodesOwned.end())
            {
                element_owned = true;
                rElementsOwned.insert(element_number);
            }
            else
            {
                temp_halo_nodes.insert(element_data.NodeIndices[i]);
            }
        }

        if (element_owned)
        {
            rHaloNodesOwned.insert(temp_halo_nodes.begin(), temp_halo_nodes.end());
        }
    }

    rMeshReader.Reset();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    this->mMeshFileBaseName = rMeshReader.GetMeshFileBaseName();
    mTotalNumElements = rMeshReader.GetNumElements();
    mTotalNumNodes = rMeshReader.GetNumNodes();
    mTotalNumBoundaryElements = rMeshReader.GetNumFaces();

    std::set<unsigned> nodes_owned;
    std::set<unsigned> halo_nodes_owned;
    std::set<unsigned> elements_owned;
    std::vector<unsigned> proc_offsets;//(PetscTools::GetNumProcs());

    ComputeMeshPartitioning(rMeshReader, nodes_owned, halo_nodes_owned, elements_owned, proc_offsets);

    // Reserve memory
    this->mElements.reserve(elements_owned.size());
    this->mNodes.reserve(nodes_owned.size());
    
    if ( rMeshReader.IsFileFormatBinary() )
    {
        std::vector<double> coords;
        // Binary : load only the nodes which are needed
        for (std::set<unsigned>::const_iterator it=nodes_owned.begin(); it!=nodes_owned.end(); it++)
        { 
            //Loop over wholey-owned nodes
            unsigned global_node_index=*it;
            coords = rMeshReader.GetNode(global_node_index);
            RegisterNode(global_node_index);
            this->mNodes.push_back(new Node<SPACE_DIM>(global_node_index, coords, false));
        }
        for (std::set<unsigned>::const_iterator it=halo_nodes_owned.begin(); it!=halo_nodes_owned.end(); it++)
        { 
            //Loop over halo-owned nodes
            unsigned global_node_index=*it;
            coords = rMeshReader.GetNode(global_node_index);
            RegisterHaloNode(global_node_index);
            mHaloNodes.push_back(new Node<SPACE_DIM>(global_node_index, coords, false));
        }
    }
    else
    {
        // Ascii : Sequentially load all the nodes and store those owned (or halo-owned) by the process
        for (unsigned node_index=0; node_index < mTotalNumNodes; node_index++)
        {
            std::vector<double> coords;
            /// \todo: assert the node is not considered both owned and halo-owned. Remove continue statement few lines below then.
            coords = rMeshReader.GetNextNode();
    
            // The node is owned by the processor
            if (nodes_owned.find(node_index) != nodes_owned.end())
            {
                RegisterNode(node_index);
                this->mNodes.push_back(new Node<SPACE_DIM>(node_index, coords, false));
            }
    
            // The node is a halo node in this processor
            if (halo_nodes_owned.find(node_index) != halo_nodes_owned.end())
            {
                RegisterHaloNode(node_index);
                mHaloNodes.push_back(new Node<SPACE_DIM>(node_index, coords, false));
            }
        }
    }
    
    if ( rMeshReader.IsFileFormatBinary() )
    {
        // Binary format, we loop only over the elements we have been assigned
        for (std::set<unsigned>::const_iterator elem_it=elements_owned.begin(); elem_it!=elements_owned.end(); elem_it++)
        {
            unsigned global_element_index=*elem_it;
            ElementData element_data;   
            element_data = rMeshReader.GetElementData(global_element_index);
            
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                //because we have populated mNodes and mHaloNodes above, we can now use this method, which should never throw
                nodes.push_back(this->GetAnyNode(element_data.NodeIndices[j]));
            }

            RegisterElement(global_element_index);
            Element<ELEMENT_DIM,SPACE_DIM>* p_element = new Element<ELEMENT_DIM,SPACE_DIM>(global_element_index, nodes);
            this->mElements.push_back(p_element);
            
            if (rMeshReader.GetNumElementAttributes() > 0)
            {
                assert(rMeshReader.GetNumElementAttributes() == 1);
                unsigned attribute_value = element_data.AttributeValue;
                p_element->SetRegion(attribute_value);
            }
        }
    }
    else
    {
        // Load the elements owned by the processor
        for (unsigned element_index=0; element_index < mTotalNumElements; element_index++)
        {
            ElementData element_data;
    
            element_data = rMeshReader.GetNextElementData();
            
            // The element is owned by the processor
            if (elements_owned.find(element_index) != elements_owned.end())
            {
                std::vector<Node<SPACE_DIM>*> nodes;
                for (unsigned j=0; j<ELEMENT_DIM+1; j++)
                {
                    //because we have populated mNodes and mHaloNodes above, we can now use this method, which should never throw
                    nodes.push_back(this->GetAnyNode(element_data.NodeIndices[j]));
                }
    
                RegisterElement(element_index);
    
                Element<ELEMENT_DIM,SPACE_DIM>* p_element = new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes);
                this->mElements.push_back(p_element);
    
                if (rMeshReader.GetNumElementAttributes() > 0)
                {
                    assert(rMeshReader.GetNumElementAttributes() == 1);
                    unsigned attribute_value = element_data.AttributeValue;
                    p_element->SetRegion(attribute_value);
                }
            }
        }
    }

    // Boundary nodes and elements
    try
    {
        for (unsigned face_index=0; face_index<mTotalNumBoundaryElements; face_index++)
        {
            ElementData face_data = rMeshReader.GetNextFaceData();
            std::vector<unsigned> node_indices = face_data.NodeIndices;
    
            bool own = false;
    
            for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
            {
                // if I own this node
                if (mNodesMapping.find(node_indices[node_index]) != mNodesMapping.end())
                {
                    own = true;
                    break;
                }
            }
    
            if (!own)
            {
                continue;
            }
    
            // Determine if this is a boundary face
            std::set<unsigned> containing_element_indices; // Elements that contain this face
            std::vector<Node<SPACE_DIM>*> nodes;
    
            for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
            {
                //because we have populated mNodes and mHaloNodes above, we can now use this method, 
                //which SHOULD never throw (but it does).
                try
                {
                    nodes.push_back(this->GetAnyNode(node_indices[node_index]));
                }
                catch (Exception &e)
                {
                    std::stringstream message;
                    message << "Face does not appear in element file (Face " << face_index << " in "<<this->mMeshFileBaseName<< ")";
                    EXCEPTION(message.str().c_str());
                }
            }
    
            // This is a boundary face
            // Ensure all its nodes are marked as boundary nodes
            for (unsigned j=0; j<nodes.size(); j++)
            {
                if (!nodes[j]->IsBoundaryNode())
                {
                    nodes[j]->SetAsBoundaryNode();
                    this->mBoundaryNodes.push_back(nodes[j]);
                }
                // Register the index that this bounday element will have with the node
                nodes[j]->AddBoundaryElement(face_index);
            }
    
            RegisterBoundaryElement(face_index);
            BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* p_boundary_element = new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(face_index, nodes);
            this->mBoundaryElements.push_back(p_boundary_element);
    
            if (rMeshReader.GetNumFaceAttributes() > 0)
            {
                assert(rMeshReader.GetNumFaceAttributes() == 1);
                unsigned attribute_value = face_data.AttributeValue;
                p_boundary_element->SetRegion(attribute_value);
            }
        }
        }
    catch (Exception &e)
    {
        PetscTools::ReplicateException(true); //Bad face exception
        throw e;
    }
    PetscTools::ReplicateException(false);

    if (mMetisPartitioning != DUMB && !PetscTools::IsSequential())
    {
        assert(this->mNodesPermutation.size() != 0);
        // We reorder so that each process owns a contiguous set of the indices and we can then build a distributed vector factory. 
        ReorderNodes();

        unsigned num_owned;
        unsigned rank = PetscTools::GetMyRank();
        if ( !PetscTools::AmTopMost() )
        {
            num_owned =  proc_offsets[rank+1]-proc_offsets[rank];
        }
        else
        {
            num_owned = mTotalNumNodes - proc_offsets[rank];
        }

        assert(!this->mpDistributedVectorFactory);
        this->mpDistributedVectorFactory = new DistributedVectorFactory(this->GetNumNodes(), num_owned);
    }
    else
    {
        // Dumb or sequential partition
        assert(this->mpDistributedVectorFactory);
    }
    rMeshReader.Reset();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalNodes() const
{
    return this->mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalElements() const
{
    return this->mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalBoundaryElements() const
{
    return this->mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mTotalNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes() const
{
    return mTotalNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mTotalNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PartitionType DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetPartitionType() const
{
    return mMetisPartitioning;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return mTotalNumBoundaryElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetHaloNodeIndices(std::vector<unsigned>& rHaloIndices) const
{
    //Make sure the output vector is empty
    rHaloIndices.clear();
    for (unsigned i=0; i<mHaloNodes.size(); i++)
    {
        rHaloIndices.push_back(mHaloNodes[i]->GetIndex());
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    // This method exists just to keep compatibility with TetrahedralMesh.
    // In parallel, you only create the data structures for the elements you have been assigned.
    // Therefore, all the local elements are owned by the processor (obviously...)
    // We don't care about "hi" and "lo"
    assert(hi>=lo);
    for (unsigned element_index=0; element_index<this->mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element=this->mElements[element_index];
        p_element->SetOwnership(true);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateDesignatedOwnershipOfElement( unsigned elementIndex ) const
{
    try
    {
        unsigned tie_break_index = this->GetElement(elementIndex)->GetNodeGlobalIndex(0); // throws an exception if we don't own the element
        SolveNodeMapping(tie_break_index);      // throws an exception if we don't own node 0
        return true;
    }
    catch(Exception e)      // either we don't own the element or we don't own node 0 of a shared element
    {
        return false;
    }
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateDesignatedOwnershipOfBoundaryElement( unsigned faceIndex ) const
{
    try
    {
        unsigned tie_break_index = this->GetBoundaryElement(faceIndex)->GetNodeGlobalIndex(0); // throws an exception if we don't own the element
        SolveNodeMapping(tie_break_index);      // throws an exception if we don't own node 0
        return true;
    }
    catch(Exception e)      // either we don't own the element or we don't own node 0 of a shared element
    {
        return false;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterNode(unsigned index)
{
    mNodesMapping[index] = this->mNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterHaloNode(unsigned index)
{
    mHaloNodesMapping[index] = mHaloNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterElement(unsigned index)
{
    mElementsMapping[index] = this->mElements.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterBoundaryElement(unsigned index)
{
    mBoundaryElementsMapping[index] = this->mBoundaryElements.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    std::map<unsigned, unsigned>::const_iterator node_position = mNodesMapping.find(index);

    if (node_position == mNodesMapping.end())
    {
        std::stringstream message;
        message << "Requested node " << index << " does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());
    }
    return node_position->second;
}

//template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveHaloNodeMapping(unsigned index)
//{
//    assert(mHaloNodesMapping.find(index) != mHaloNodesMapping.end());
//    return mHaloNodesMapping[index];
//}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    std::map<unsigned, unsigned>::const_iterator element_position = mElementsMapping.find(index);

    if (element_position == mElementsMapping.end())
    {
        std::stringstream message;
        message << "Requested element " << index << " does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());
    }

    return element_position->second;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    std::map<unsigned, unsigned>::const_iterator boundary_element_position = mBoundaryElementsMapping.find(index);

    if (boundary_element_position == mBoundaryElementsMapping.end())
    {
        std::stringstream message;
        message << "Requested boundary element " << index << " does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());
    }

    return boundary_element_position->second;
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM> * DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetAnyNode(unsigned index) const
{
    std::map<unsigned, unsigned>::const_iterator node_position;
    //First search the halo
    if ((node_position=mHaloNodesMapping.find(index)) != mHaloNodesMapping.end())
    {
        return mHaloNodes[node_position->second];
    }
    //First search the owned node
    if ((node_position=mNodesMapping.find(index)) != mNodesMapping.end())
    {
        //Found an owned node
        return this->mNodes[node_position->second];
    }
    //Not here
    std::stringstream message;
    message << "Requested node/halo " << index << " does not belong to processor " << PetscTools::GetMyRank();
    EXCEPTION(message.str().c_str());
}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DumbNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                                                           std::set<unsigned>& rNodesOwned)
{
    assert(!this->mpDistributedVectorFactory);
    this->mpDistributedVectorFactory = new DistributedVectorFactory(mTotalNumNodes);
    for (unsigned node_index = this->mpDistributedVectorFactory->GetLow();
         node_index < this->mpDistributedVectorFactory->GetHigh();
         node_index++)
    {
         rNodesOwned.insert(node_index);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MetisBinaryNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                                                                  std::set<unsigned>& rNodesOwned,
                                                                                  std::vector<unsigned>& rProcessorsOffset)
{
    assert(!PetscTools::IsSequential());
    #define COVERAGE_IGNORE
    assert( ELEMENT_DIM==2 || ELEMENT_DIM==3 ); // Metis works with triangles and tetras
    #undef COVERAGE_IGNORE

    unsigned num_procs = PetscTools::GetNumProcs();

    // Open a file for the elements
    OutputFileHandler handler("");

    // Filenames
    std::string basename = "metis.mesh";
    std::stringstream output_file;
    output_file << basename << ".npart." << num_procs;
    std::string nodes_per_proc_file = basename + ".nodesperproc";

    // Only the master process should do IO and call METIS
    std::string full_path = handler.GetOutputDirectoryFullPath();

    if (PetscTools::AmMaster())
    {
        /*
         *  Create input file for METIS
         */
        out_stream metis_file=handler.OpenOutputFile(basename);

        // File header
        (*metis_file) << this->GetNumElements() << "\t";
        if (ELEMENT_DIM==2)
        {
            (*metis_file) << 1 << "\n"; //1 is Metis speak for triangles
        }
        else
        {
            (*metis_file) << 2 << "\n"; //2 is Metis speak for tetrahedra
        }

        // Graph representation of the mesh
        for (unsigned element_number = 0; element_number < mTotalNumElements; element_number++)
        {
            ElementData element_data = rMeshReader.GetNextElementData();

            for (unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                    (*metis_file) << element_data.NodeIndices[i] + 1 << "\t";
            }
            (*metis_file) << "\n";
        }
        metis_file->close();

        rMeshReader.Reset();

        /*
         * Call METIS binary to perform the partitioning.
         * It will output a file called metis.mesh.npart.numProcs
         */
        std::stringstream permute_command;
        permute_command <<  "./bin/partdmesh "
                        <<  full_path
                        <<  basename << " "
                        <<  num_procs
                        <<  " > /dev/null";

        // METIS doesn't return 0 after a successful execution
        IGNORE_RET(system, permute_command.str());
    }

    /*
     * Wait for the permutation to be available
     */
    PetscTools::Barrier("DistributedTetrahedralMesh::MetisBinaryNodePartitioning");

    /*
     *  Read partition file and compute local node ownership and processors offset
     */
    std::ifstream partition_stream;
    full_path +=  output_file.str();

    partition_stream.open(full_path.c_str());
    assert(partition_stream.is_open());

    assert(rProcessorsOffset.size() == 0); // Making sure the vector is empty. After calling resize() only newly created memory will be initialised to 0.
    rProcessorsOffset.resize(PetscTools::GetNumProcs(), 0);

    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        unsigned part_read;

        partition_stream >> part_read;

        // METIS output says I own this node
        if (part_read == PetscTools::GetMyRank())
        {
            rNodesOwned.insert(node_index);
        }

        // Offset is defined as the first node owned by a processor. We compute it incrementally.
        // i.e. if node_index belongs to proc 3 (of 6) we have to shift the processors 4, 5, and 6
        // offset a position
        for (unsigned proc=part_read+1; proc<PetscTools::GetNumProcs(); proc++)
        {
            rProcessorsOffset[proc]++;
        }
    }

    /*
     *  Once we know the offsets we can compute the permutation vector
     */
    // Move stream pointer to the beginning of the file
    partition_stream.seekg (0, std::ios::beg);

    std::vector<unsigned> local_index(PetscTools::GetNumProcs(), 0);
    this->mNodesPermutation.resize(this->GetNumNodes());

    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        unsigned part_read;

        partition_stream >> part_read;

        this->mNodesPermutation[node_index] = rProcessorsOffset[part_read] + local_index[part_read];

        local_index[part_read]++;
    }

    partition_stream.close();

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MetisLibraryNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                                                                  std::set<unsigned>& rNodesOwned,
                                                                                  std::vector<unsigned>& rProcessorsOffset)
{
    assert(!PetscTools::IsSequential());

    assert(ELEMENT_DIM==2 || ELEMENT_DIM==3); // Metis works with triangles and tetras

    int nn = rMeshReader.GetNumNodes();
    idxtype* npart = new idxtype[nn];
    assert(npart != NULL);

    //Only the master process will access the element data and perform the partitioning
    if (PetscTools::AmMaster())
    {
        int ne = rMeshReader.GetNumElements();
        idxtype* elmnts = new idxtype[ne * (ELEMENT_DIM+1)];
        assert(elmnts != NULL);

        unsigned counter=0;
        for (unsigned element_number = 0; element_number < mTotalNumElements; element_number++)
        {
            ElementData element_data = rMeshReader.GetNextElementData();

            for (unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                elmnts[counter++] = element_data.NodeIndices[i];
            }
        }
        rMeshReader.Reset();

        int etype;

        switch (ELEMENT_DIM)
        {
            case 2:
                etype = 1; //1 is Metis speak for triangles
                break;
            case 3:
                etype = 2; //2 is Metis speak for tetrahedra
                break;
            default:
                NEVER_REACHED;
        }

        int numflag = 0; //0 means C-style numbering is assumed
        int nparts = PetscTools::GetNumProcs();
        int edgecut;
        idxtype* epart = new idxtype[ne];
        assert(epart != NULL);

        METIS_PartMeshNodal(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);//, wgetflag, vwgt);
        delete[] elmnts;
        delete[] epart;
    }

    //Here's the new bottle-neck: share all the node ownership data
    //idxtype is normally int (see metis-4.0/Lib/struct.h 17-22)
    assert(sizeof(idxtype) == sizeof(int));
    MPI_Bcast(npart /*data*/, nn /*size*/, MPI_INT, 0 /*From Master*/, PETSC_COMM_WORLD);

    assert(rProcessorsOffset.size() == 0); // Making sure the vector is empty. After calling resize() only newly created memory will be initialised to 0.
    rProcessorsOffset.resize(PetscTools::GetNumProcs(), 0);

    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        unsigned part_read = npart[node_index];

        // METIS output says I own this node
        if (part_read == PetscTools::GetMyRank())
        {
            rNodesOwned.insert(node_index);
        }

        // Offset is defined as the first node owned by a processor. We compute it incrementally.
        // i.e. if node_index belongs to proc 3 (of 6) we have to shift the processors 4, 5, and 6
        // offset a position.
        for (unsigned proc=part_read+1; proc<PetscTools::GetNumProcs(); proc++)
        {
            rProcessorsOffset[proc]++;
        }
    }

    /*
     *  Once we know the offsets we can compute the permutation vector
     */
    std::vector<unsigned> local_index(PetscTools::GetNumProcs(), 0);

    this->mNodesPermutation.resize(this->GetNumNodes());

    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        unsigned part_read = npart[node_index];

        this->mNodesPermutation[node_index] = rProcessorsOffset[part_read] + local_index[part_read];

        local_index[part_read]++;
    }

    delete[] npart;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReorderNodes()
{
    assert(!PetscTools::IsSequential());

    // Need to rebuild global-local maps
    mNodesMapping.clear();
    mHaloNodesMapping.clear();

    // Update indices
    for (unsigned index=0; index<this->mNodes.size(); index++)
    {
        unsigned old_index = this->mNodes[index]->GetIndex();
        unsigned new_index = this->mNodesPermutation[old_index];

        this->mNodes[index]->SetIndex(new_index);
        mNodesMapping[new_index] = index;
    }

    for (unsigned index=0; index<mHaloNodes.size(); index++)
    {
        unsigned old_index = mHaloNodes[index]->GetIndex();
        unsigned new_index = this->mNodesPermutation[old_index];

        mHaloNodes[index]->SetIndex(new_index);
        mHaloNodesMapping[new_index] = index;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructLinearMesh(unsigned width)
{
    assert(ELEMENT_DIM == 1);

     //Check that there are enough nodes to make the parallelisation worthwhile
    if (width<2 || width+1 < PetscTools::GetNumProcs())
    {
        EXCEPTION("There aren't enough nodes to make parallelisation worthwhile");
    }
    //Use dumb partition so that archiving doesn't permute anything
    mMetisPartitioning=DUMB;
    mTotalNumNodes=width+1;
    mTotalNumBoundaryElements=2u;
    mTotalNumElements=width;

    //Use DistributedVectorFactory to make a dumb partition of the nodes
    assert(!this->mpDistributedVectorFactory);
    this->mpDistributedVectorFactory = new DistributedVectorFactory(mTotalNumNodes);

    unsigned lo_node=this->mpDistributedVectorFactory->GetLow();
    unsigned hi_node=this->mpDistributedVectorFactory->GetHigh();

    if (!PetscTools::AmMaster())
    {
        //Allow for a halo node
        lo_node--;
    }
    if (!PetscTools::AmTopMost())
    {
        //Allow for a halo node
        hi_node++;
    }
    Node<SPACE_DIM>* p_old_node=NULL;
    for (unsigned node_index=lo_node; node_index<hi_node; node_index++)
    {
        // create node or halo-node
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, node_index==0 || node_index==width, node_index);
        if (node_index<this->mpDistributedVectorFactory->GetLow() ||
            node_index==this->mpDistributedVectorFactory->GetHigh() )
        {
            //Beyond left or right it's a halo node
            RegisterHaloNode(node_index);
            mHaloNodes.push_back(p_node);
        }
        else
        {
            RegisterNode(node_index);
            this->mNodes.push_back(p_node); // create node

            //A boundary face has to be wholely owned by the process
            //Since, when ELEMENT_DIM>1 we have *at least* boundary node as a non-halo
            if (node_index==0) // create left boundary node and boundary element
            {
                this->mBoundaryNodes.push_back(p_node);
                RegisterBoundaryElement(0);
                this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(0, p_node) );
            }
            if (node_index==width) // create right boundary node and boundary element
            {
                this->mBoundaryNodes.push_back(p_node);
                RegisterBoundaryElement(1);
                this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(1, p_node) );
            }
            }
        if (node_index>lo_node) // create element
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            nodes.push_back(p_old_node);
            nodes.push_back(p_node);
            RegisterElement(node_index-1);
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(node_index-1, nodes) );
        }
        //Keep track of the node which we've just created
        p_old_node=p_node;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMesh(unsigned width, unsigned height, bool stagger)
{
    assert(SPACE_DIM == 2);
    assert(ELEMENT_DIM == 2);
    //Check that there are enough nodes to make the parallelisation worthwhile
    if (height<2 || height+1 < PetscTools::GetNumProcs())
    {
        EXCEPTION("There aren't enough nodes to make parallelisation worthwhile");
    }
    //Use dumb partition so that archiving doesn't permute anything
    mMetisPartitioning=DUMB;

    mTotalNumNodes=(width+1)*(height+1);
    mTotalNumBoundaryElements=(width+height)*2;
    mTotalNumElements=width*height*2;

    //Use DistributedVectorFactory to make a dumb partition of space
    DistributedVectorFactory y_partition(height+1);
    unsigned lo_y = y_partition.GetLow();
    unsigned hi_y = y_partition.GetHigh();
    //Dumb partition of nodes has to be such that each process gets complete slices
    assert(!this->mpDistributedVectorFactory);
    this->mpDistributedVectorFactory = new DistributedVectorFactory(mTotalNumNodes, (width+1)*y_partition.GetLocalOwnership());

    if (!PetscTools::AmMaster())
    {
        //Allow for a halo node
        lo_y--;
    }
    if (!PetscTools::AmTopMost())
    {
        //Allow for a halo node
        hi_y++;
    }

    //Construct the nodes
    for (unsigned j=lo_y; j<hi_y; j++)
    {
        for (unsigned i=0; i<width+1; i++)
        {
            bool is_boundary=false;
            if (i==0 || j==0 || i==width || j==height)
            {
                is_boundary=true;
            }
            unsigned global_node_index=((width+1)*(j) + i); //Verified from sequential
            Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(global_node_index, is_boundary, i, j);
            if (j<y_partition.GetLow() || j==y_partition.GetHigh() )
            {
                //Beyond left or right it's a halo node
                RegisterHaloNode(global_node_index);
                mHaloNodes.push_back(p_node);
            }
            else
            {
                RegisterNode(global_node_index);
                this->mNodes.push_back(p_node);
            }
            if (is_boundary)
            {
                this->mBoundaryNodes.push_back(p_node);
            }
        }
    }

    //Construct the boundary elements
    unsigned belem_index;
    //Top
    if (PetscTools::AmTopMost())
    {
       for (unsigned i=0; i<width; i++)
       {
            std::vector<Node<SPACE_DIM>*> nodes;
            nodes.push_back(GetAnyNode( height*(width+1)+i ));
            nodes.push_back(GetAnyNode( height*(width+1)+i+1 ));
            belem_index=i;
            RegisterBoundaryElement(belem_index);
            this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index,nodes));
        }
    }

    //Right
    for (unsigned j=lo_y+1; j<hi_y; j++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(GetAnyNode( (width+1)*(j+1)-1 ));
        nodes.push_back(GetAnyNode( (width+1)*j-1 ));
        belem_index=width+j-1;
        RegisterBoundaryElement(belem_index);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index,nodes));
    }

    //Bottom
    if (PetscTools::AmMaster())
    {
        for (unsigned i=0; i<width; i++)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            nodes.push_back(GetAnyNode( i+1 ));
            nodes.push_back(GetAnyNode( i ));
            belem_index=width+height+i;
            RegisterBoundaryElement(belem_index);
            this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index,nodes));
        }
    }

    //Left
    for (unsigned j=lo_y; j<hi_y-1; j++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(GetAnyNode( (width+1)*(j+1) ));
        nodes.push_back(GetAnyNode( (width+1)*(j) ));
        belem_index=2*width+height+j;
        RegisterBoundaryElement(belem_index);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index,nodes));
    }


    //Construct the elements
    unsigned elem_index;
    for (unsigned j=lo_y; j<hi_y-1; j++)
    {
        for (unsigned i=0; i<width; i++)
        {
            unsigned parity=(i+(height-j))%2;//Note that parity is measured from the top-left (not bottom left) for historical reasons
            unsigned nw=(j+1)*(width+1)+i; //ne=nw+1
            unsigned sw=(j)*(width+1)+i;   //se=sw+1
            std::vector<Node<SPACE_DIM>*> upper_nodes;
            upper_nodes.push_back(GetAnyNode( nw ));
            upper_nodes.push_back(GetAnyNode( nw+1 ));
            if (stagger==false  || parity == 1)
            {
                upper_nodes.push_back(GetAnyNode( sw+1 ));
            }
            else
            {
                upper_nodes.push_back(GetAnyNode( sw ));
            }
            elem_index=2*(j*width+i);
            RegisterElement(elem_index);
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index,upper_nodes));
            std::vector<Node<SPACE_DIM>*> lower_nodes;
            lower_nodes.push_back(GetAnyNode( sw+1 ));
            lower_nodes.push_back(GetAnyNode( sw ));
            if (stagger==false  ||parity == 1)
            {
                lower_nodes.push_back(GetAnyNode( nw ));
            }
            else
            {
                lower_nodes.push_back(GetAnyNode( nw+1 ));
            }
            elem_index++;
            RegisterElement(elem_index);
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index,lower_nodes));
        }
    }

}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructCuboid(unsigned width,
        unsigned height,
        unsigned depth)
{
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 3);
    //Check that there are enough nodes to make the parallelisation worthwhile
    if (depth<2 || depth+1 < PetscTools::GetNumProcs())
    {
        EXCEPTION("There aren't enough nodes to make parallelisation worthwhile");
    }

    //Use dumb partition so that archiving doesn't permute anything
    mMetisPartitioning=DUMB;

    mTotalNumNodes=(width+1)*(height+1)*(depth+1);
    mTotalNumBoundaryElements=((width*height)+(width*depth)+(height*depth))*4;//*2 for top-bottom, *2 for tessellating each unit square
    mTotalNumElements=width*height*depth*6;

    //Use DistributedVectorFactory to make a dumb partition of space
    DistributedVectorFactory z_partition(depth+1);
    unsigned lo_z = z_partition.GetLow();
    unsigned hi_z = z_partition.GetHigh();
    //Dumb partition of nodes has to be such that each process gets complete slices
    assert(!this->mpDistributedVectorFactory);
    this->mpDistributedVectorFactory = new DistributedVectorFactory(mTotalNumNodes, (width+1)*(height+1)*z_partition.GetLocalOwnership());
    if (!PetscTools::AmMaster())
    {
        //Allow for a halo node
        lo_z--;
    }
    if (!PetscTools::AmTopMost())
    {
        //Allow for a halo node
        hi_z++;
    }

    //Construct the nodes
    unsigned global_node_index;
    for (unsigned k=lo_z; k<hi_z; k++)
    {
        for (unsigned j=0; j<height+1; j++)
        {
            for (unsigned i=0; i<width+1; i++)
            {
                bool is_boundary = false;
                if (i==0 || j==0 || k==0 || i==width || j==height || k==depth)
                {
                    is_boundary = true;
                }
                global_node_index = (k*(height+1)+j)*(width+1)+i;

                Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(global_node_index, is_boundary, i, j, k);

                if (k<z_partition.GetLow() || k==z_partition.GetHigh() )
                {
                    //Beyond left or right it's a halo node
                    RegisterHaloNode(global_node_index);
                    mHaloNodes.push_back(p_node);
                }
                else
                {
                    RegisterNode(global_node_index);
                    this->mNodes.push_back(p_node);
                }

                if (is_boundary)
                {
                    this->mBoundaryNodes.push_back(p_node);
                }
            }
        }
    }

    // Construct the elements

    unsigned element_nodes[6][4] = {{0, 1, 5, 7}, {0, 1, 3, 7},
                                        {0, 2, 3, 7}, {0, 2, 6, 7},
                                        {0, 4, 6, 7}, {0, 4, 5, 7}};
    std::vector<Node<SPACE_DIM>*> tetrahedra_nodes;

    for (unsigned k=lo_z; k<hi_z-1; k++)
    {
        unsigned belem_index = 0;
        if (k != 0)
        {
            // height*width squares on upper face, k layers of 2*height+2*width square aroun
            belem_index =   2*(height*width+k*2*(height+width));
        }

        for (unsigned j=0; j<height; j++)
        {
            for (unsigned i=0; i<width; i++)
            {
                // Compute the nodes' index
                unsigned global_node_indices[8];
                unsigned local_node_index = 0;

                for (unsigned z = 0; z < 2; z++)
                {
                    for (unsigned y = 0; y < 2; y++)
                    {
                        for (unsigned x = 0; x < 2; x++)
                        {
                            global_node_indices[local_node_index] = i+x+(width+1)*(j+y+(height+1)*(k+z));

                            local_node_index++;
                        }
                    }
                }

                for (unsigned m = 0; m < 6; m++)
                {
                    // Tetrahedra #m

                    tetrahedra_nodes.clear();

                    for (unsigned n = 0; n < 4; n++)
                    {
                        tetrahedra_nodes.push_back(GetAnyNode( global_node_indices[element_nodes[m][n]] ));
                    }
                    unsigned elem_index = 6 * ((k*height+j)*width+i)+m;
                    RegisterElement(elem_index);
                    this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index, tetrahedra_nodes));
                }

                //Are we at a boundary?
                std::vector<Node<SPACE_DIM>*> triangle_nodes;

                if (i == 0) //low face at x==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[0] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[2] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[6] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[0] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[6] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[4] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (i == width-1) //high face at x=width
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[1] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[5] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[7] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[1] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[7] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[3] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == 0) //low face at y==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[0] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[5] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[1] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[0] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[4] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[5] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == height-1) //high face at y=height
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[2] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[3] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[7] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[2] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[7] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[6] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == 0) //low face at z==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[0] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[3] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[2] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[0] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[1] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[3] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == depth-1) //high face at z=depth
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[4] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[7] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[5] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[4] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[6] ));
                    triangle_nodes.push_back(GetAnyNode( global_node_indices[7] ));
                    RegisterBoundaryElement(belem_index);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
            }//i
        }//j
    }//k
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Scale(const double xFactor, const double yFactor, const double zFactor)
{
    //Base class scale (scales node positions)
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Scale(xFactor, yFactor, zFactor);
    //Scales halos
    for (unsigned i=0; i<mHaloNodes.size(); i++)
    {
        c_vector<double, SPACE_DIM>& r_location = mHaloNodes[i]->rGetModifiableLocation();
        if (SPACE_DIM>=3)
        {
            r_location[2] *= zFactor;
        }
        if (SPACE_DIM>=2)
        {
            r_location[1] *= yFactor;
        }
        r_location[0] *= xFactor;
    }

}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class DistributedTetrahedralMesh<1,1>;
template class DistributedTetrahedralMesh<1,2>;
template class DistributedTetrahedralMesh<1,3>;
template class DistributedTetrahedralMesh<2,2>;
template class DistributedTetrahedralMesh<2,3>;
template class DistributedTetrahedralMesh<3,3>;
