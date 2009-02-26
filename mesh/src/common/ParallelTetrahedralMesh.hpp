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

#ifndef PARALLELTETRAHEDRALMESH_HPP_
#define PARALLELTETRAHEDRALMESH_HPP_

#include "AbstractMesh.hpp"
#include "AbstractMeshReader.hpp"
#include "Element.hpp"
#include "BoundaryElement.hpp"
#include "Node.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"

/*
 *  The following definition fixes an odd incompatibility of METIS 4.0 and Chaste. Since 
 * the library was compiled with a plain-C compiler, it fails to link using a C++ compiler.
 * Note that METIS 4.0 fails to compile with g++ or icpc, so a C compiler should be used.
 *  
 * Somebody had this problem before: http://www-users.cs.umn.edu/~karypis/.discus/messages/15/113.html?1119486445 
 *
 * Note that it is necessary to define the function header before the #include statement. 
*/ 
extern "C" {
extern void METIS_PartMeshNodal(int*, int*, int*, int*, int*, int*, int*, int*, int*);
};
#include "metis.h"

 
/**
 * \todo explicit instantiation
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ParallelTetrahedralMesh : public AbstractMesh< ELEMENT_DIM, SPACE_DIM>
{
    
    friend class TestParallelTetrahedralMesh;

public:

    typedef enum
    {
        DUMB=0,
        METIS_BINARY,
        METIS_LIBRARY
    } PartitionType;            
    
private: 

    unsigned mTotalNumElements;
    unsigned mTotalNumBoundaryElements; 
    unsigned mTotalNumNodes;
    
    std::vector<Node<SPACE_DIM>* > mHaloNodes;
    
    std::map<unsigned, unsigned> mNodesMapping;
    std::map<unsigned, unsigned> mHaloNodesMapping;    
    std::map<unsigned, unsigned> mElementsMapping;
    std::map<unsigned, unsigned> mBoundaryElementsMapping;        
    
    PartitionType mMetisPartitioning;
        
public:        

    ParallelTetrahedralMesh(PartitionType metisPartitioning=METIS_LIBRARY);

    virtual ~ParallelTetrahedralMesh();

    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                 bool cullInternalFaces=false);

    unsigned GetNumLocalNodes() const;
    unsigned GetNumLocalElements() const;

    unsigned GetNumNodes() const;
    unsigned GetNumElements() const;    
    unsigned GetNumBoundaryElements() const;    
    
    void SetElementOwnerships(unsigned lo, unsigned hi);

private:

    void RegisterNode(unsigned index);
    void RegisterHaloNode(unsigned index);
    void RegisterElement(unsigned index);
    void RegisterBoundaryElement(unsigned index);

    unsigned SolveNodeMapping(unsigned index) const;
    unsigned SolveHaloNodeMapping(unsigned index);
    unsigned SolveElementMapping(unsigned index) const;            
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    void ComputeMeshPartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                               std::set<unsigned>& rNodesOwned,
                               std::set<unsigned>& rHaloNodesOwned,
                               std::set<unsigned>& rElementsOwned,
                               std::vector<unsigned>& rProcessorsOffset,
                               std::vector<unsigned>& rNodePermutation) const;    
    
    void DumbNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                              std::set<unsigned>& rNodesOwned) const;

    void MetisBinaryNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                                     std::set<unsigned>& rNodesOwned, std::vector<unsigned>& rProcessorsOffset,
                                     std::vector<unsigned>& rNodePermutation) const;

    void MetisLibraryNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                                     std::set<unsigned>& rNodesOwned, std::vector<unsigned>& rProcessorsOffset,
                                     std::vector<unsigned>& rNodePermutation) const;
                                     
    void ReorderNodes(std::vector<unsigned>& rNodePermutation);    
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ParallelTetrahedralMesh(PartitionType metisPartitioning)
    : mMetisPartitioning(metisPartitioning)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~ParallelTetrahedralMesh()
{
    for (unsigned i=0; i<this->mHaloNodes.size(); i++)
    {
        delete this->mHaloNodes[i];
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ComputeMeshPartitioning(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    std::set<unsigned>& rNodesOwned,
    std::set<unsigned>& rHaloNodesOwned,
    std::set<unsigned>& rElementsOwned,
    std::vector<unsigned>& rProcessorsOffset,
    std::vector<unsigned>& rNodePermutation) const
{
    ///\todo: add a timing event for the partitioning
    
    if (mMetisPartitioning==METIS_BINARY && PetscTools::NumProcs() > 1)
    {        
        MetisBinaryNodePartitioning(rMeshReader, rNodesOwned, rProcessorsOffset, rNodePermutation);                 
    }
    else if (mMetisPartitioning==METIS_LIBRARY && PetscTools::NumProcs() > 1)
    {
        MetisLibraryNodePartitioning(rMeshReader, rNodesOwned, rProcessorsOffset, rNodePermutation);
    }
    else    
    {
        DumbNodePartitioning(rMeshReader, rNodesOwned);
    }
        
    for(unsigned element_number = 0; element_number < mTotalNumElements; element_number++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();

        bool element_owned = false;
        std::set<unsigned> temp_halo_nodes;
        
        for(unsigned i=0; i<ELEMENT_DIM+1; i++)
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
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    bool cullInternalFaces)
{
    
    if(ELEMENT_DIM==1)
    {
        cullInternalFaces = true;
    }
    
    // ticket #922: Face culling is not supported in parallel yet
    assert(!cullInternalFaces);    
    
    mTotalNumElements = rMeshReader.GetNumElements();
    mTotalNumNodes = rMeshReader.GetNumNodes();
    mTotalNumBoundaryElements = rMeshReader.GetNumFaces();   
        
    std::set<unsigned> nodes_owned;
    std::set<unsigned> halo_nodes_owned;
    std::set<unsigned> elements_owned;
    std::vector<unsigned> proc_offsets;//(PetscTools::NumProcs());
    std::vector<unsigned> node_permutation;//(this->GetNumNodes());    
    
    ComputeMeshPartitioning(rMeshReader, nodes_owned, halo_nodes_owned, elements_owned, proc_offsets, node_permutation);
    
    // Reserve memory
    this->mElements.reserve(elements_owned.size());
    this->mNodes.reserve(nodes_owned.size());
        
    // Load the nodes owned by the processor
    std::vector<double> coords;
    for (unsigned node_index=0; node_index < mTotalNumNodes; node_index++)
    {
        /// \todo: assert the node is not considered both owned and halo-owned. Remove continue statement few lines below then.
        coords = rMeshReader.GetNextNode();
        
        // The node is owned by the processor
        if (nodes_owned.find(node_index) != nodes_owned.end())
        {
            RegisterNode(node_index);
            this->mNodes.push_back(new Node<SPACE_DIM>(node_index, coords, false));
            continue;
        }

        // The node is a halo node in this processor
        if (halo_nodes_owned.find(node_index) != halo_nodes_owned.end())
        {
            RegisterHaloNode(node_index);
            mHaloNodes.push_back(new Node<SPACE_DIM>(node_index, coords, false));
        }

    }

    // Load the elements owned by the processor
    for (unsigned element_index=0; element_index < mTotalNumElements; element_index++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();

        // The element is owned by the processor
        if (elements_owned.find(element_index) != elements_owned.end())
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            unsigned node_local_index;
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                if (nodes_owned.find(element_data.NodeIndices[j]) != nodes_owned.end())
                {
                    node_local_index = SolveNodeMapping(element_data.NodeIndices[j]);
                    nodes.push_back(this->mNodes[node_local_index]);
                }
                else
                {
                    node_local_index = SolveHaloNodeMapping(element_data.NodeIndices[j]);
                    nodes.push_back(this->mHaloNodes[node_local_index]);
                }                    
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
    
    // Boundary nodes and elements    
    unsigned actual_face_index = 0; 
    for (unsigned face_index=0; face_index<mTotalNumBoundaryElements; face_index++)
    {
        ElementData face_data = rMeshReader.GetNextFaceData();
        std::vector<unsigned> node_indices = face_data.NodeIndices;

        bool own = false;

        for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
        {
            // if I own this node
            if(mNodesMapping.find(node_indices[node_index]) != mNodesMapping.end())
            {
                own = true;
            }
        }
        
        if (!own )
        {
            // ticket #922: If we are culling internal faces we need to check if this is an external one before incrementing the index. 
            //              This turned to be tricky in parallel... since you can't check it in faces you don't own.  
            assert(!cullInternalFaces);
            actual_face_index++; 
            continue;
        }
            
        
        bool is_boundary_face = true;

        // Determine if this is a boundary face
        std::set<unsigned> containing_element_indices; // Elements that contain this face
        std::vector<Node<SPACE_DIM>*> nodes;

        for (unsigned node_index=0; node_index<node_indices.size(); node_index++)
        {
            // if I own this node
            if(mNodesMapping.find(node_indices[node_index]) != mNodesMapping.end())
            {
                // Add Node pointer to list for creating an element
                unsigned node_local_index = SolveNodeMapping(node_indices[node_index]);
                nodes.push_back(this->mNodes[node_local_index]);
            }

            // if I halo-own this node
            if(mHaloNodesMapping.find(node_indices[node_index]) != mHaloNodesMapping.end())
            {
                // Add Node pointer to list for creating an element
                unsigned node_local_index = SolveHaloNodeMapping(node_indices[node_index]);                
                nodes.push_back(this->mHaloNodes[node_local_index]); 
            }
            
            if(cullInternalFaces)
            {
                // Work out what elements contain this face, by taking the intersection
                // of the sets of elements containing each node in the face.
                if (node_index == 0)
                {
                    containing_element_indices = nodes[node_index]->rGetContainingElementIndices();
                }
                else
                {
                    std::set<unsigned> temp;
                    std::set_intersection(nodes[node_index]->rGetContainingElementIndices().begin(),
                                          nodes[node_index]->rGetContainingElementIndices().end(),
                                          containing_element_indices.begin(), containing_element_indices.end(),
                                          std::inserter(temp, temp.begin()));
                    containing_element_indices = temp;
                }
            }
        }
       

        if(cullInternalFaces)
        {
            // only if not 1D as this assertion does not apply to quadratic 1D meshes
            if(ELEMENT_DIM!=1)
            {
                //If the following assertion is thrown, it means that the .edge/.face file does not
                //match the .ele file -- they were generated at separate times.  Simply remove the internal
                //edges/faces by hand.
                assert(containing_element_indices.size() != 0);
            }

            // if num_containing_elements is greater than 1, it is not an boundary face
            if(containing_element_indices.size() > 1)
            {
                is_boundary_face = false;
            }
            
            // in 1D QUADRATICS, all nodes are faces, so internal nodes which don't have any
            // containing elements must also be unmarked as a boundary face
            if( (ELEMENT_DIM==1) && (containing_element_indices.size()==0))
            {
                is_boundary_face = false;
            }
        }

        if (is_boundary_face)
        {
            // This is a boundary face
            // Ensure all its nodes are marked as boundary nodes
            for (unsigned j=0; j<nodes.size(); j++)
            {
                if (!nodes[j]->IsBoundaryNode())
                {
                    nodes[j]->SetAsBoundaryNode();
                    this->mBoundaryNodes.push_back(nodes[j]);
                }
                //Register the index that this bounday element will have
                //with the node
                nodes[j]->AddBoundaryElement(actual_face_index);
            }

            RegisterBoundaryElement(actual_face_index);
            BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* p_boundary_element = new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(actual_face_index, nodes);
            this->mBoundaryElements.push_back(p_boundary_element);
            
            if (rMeshReader.GetNumFaceAttributes() > 0)
            {
                assert(rMeshReader.GetNumFaceAttributes() == 1);
                unsigned attribute_value = face_data.AttributeValue;
                p_boundary_element->SetRegion(attribute_value);
            }         
            actual_face_index++;
        }
    }
   
    if (mMetisPartitioning != DUMB && PetscTools::NumProcs()>1)
    {
        ReorderNodes(node_permutation);

        // Compute nodes per processor based on offset information
        for (unsigned num_proc=0; num_proc<PetscTools::NumProcs()-1; num_proc++)
        {
            this->mNodesPerProcessor.push_back( proc_offsets[num_proc+1]-proc_offsets[num_proc] );
        }
        
        // Entry for the last processor
        this->mNodesPerProcessor.push_back( mTotalNumNodes - proc_offsets[PetscTools::NumProcs()-1] );
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalNodes() const
{
    return this->mNodes.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalElements() const
{
    return this->mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mTotalNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mTotalNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return mTotalNumBoundaryElements;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    // all the local elements are owned by the processor (obviously...)    
    assert(hi>=lo);
    for (unsigned element_index=0; element_index<this->mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element=this->mElements[element_index];
        p_element->SetOwnership(true);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterNode(unsigned index)
{
    mNodesMapping[index] = this->mNodes.size();
}    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterHaloNode(unsigned index)
{
    mHaloNodesMapping[index] = mHaloNodes.size();    
}    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterElement(unsigned index)
{
    mElementsMapping[index] = this->mElements.size();    
}    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RegisterBoundaryElement(unsigned index)
{
    mBoundaryElementsMapping[index] = this->mBoundaryElements.size();    
}    


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    std::map<unsigned, unsigned>::const_iterator node_position = mNodesMapping.find(index); 
        
    if(node_position == mNodesMapping.end())
    {
        std::stringstream message;
        message << "Requested node " << index << " does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());        
    }
    return node_position->second;    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveHaloNodeMapping(unsigned index)
{
    assert(mHaloNodesMapping.find(index) != mHaloNodesMapping.end());
    return mHaloNodesMapping[index];    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    std::map<unsigned, unsigned>::const_iterator element_position = mElementsMapping.find(index); 
    
    if(element_position == mElementsMapping.end())
    {
        std::stringstream message;
        message << "Requested element " << index << " does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());        
    }

    return element_position->second;    
}            

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    std::map<unsigned, unsigned>::const_iterator boundary_element_position = mBoundaryElementsMapping.find(index); 
    
    if(boundary_element_position == mBoundaryElementsMapping.end())
    {
        std::stringstream message;
        message << "Requested boundary element " << index << " does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());        
    }

    return boundary_element_position->second;    
}            

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::DumbNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                                                                           std::set<unsigned>& rNodesOwned) const
{
    DistributedVector::SetProblemSize(mTotalNumNodes);
    for(DistributedVector::Iterator node_number = DistributedVector::Begin(); node_number != DistributedVector::End(); ++node_number)
    {
         rNodesOwned.insert(node_number.Global);
    }
    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MetisBinaryNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                                                                                  std::set<unsigned>& rNodesOwned, 
                                                                                  std::vector<unsigned>& rProcessorsOffset,
                                                                                  std::vector<unsigned>& rNodePermutation) const
{
    assert(PetscTools::NumProcs() > 1);
    
    assert( ELEMENT_DIM==2 || ELEMENT_DIM==3 ); // Metis works with triangles and tetras

    unsigned num_procs = PetscTools::NumProcs();
    
    // Open a file for the elements
    OutputFileHandler handler("");

    // Filenames
    std::string basename = "metis.mesh";
    std::stringstream output_file;
    output_file << basename << ".npart." << num_procs;
    std::string nodes_per_proc_file = basename + ".nodesperproc";

    // Only the master process should do IO and call METIS
    if (handler.IsMaster())
    {
        try
        {
        
            /*
             *  Create input file for METIS
             */
            out_stream metis_file=handler.OpenOutputFile(basename);
    
            // File header
            (*metis_file)<<this->GetNumElements()<<"\t";
            if (ELEMENT_DIM==2)
            {
                (*metis_file)<<1<<"\n"; //1 is Metis speak for triangles
            }
            else
            {
                (*metis_file)<<2<<"\n"; //2 is Metis speak for tetrahedra
            }
    
    
            // Graph representation of the mesh
            for(unsigned element_number = 0; element_number < mTotalNumElements; element_number++)
            {
                ElementData element_data = rMeshReader.GetNextElementData();
        
                for(unsigned i=0; i<ELEMENT_DIM+1; i++)
                {
                        (*metis_file)<<element_data.NodeIndices[i] + 1<<"\t";
                }
                (*metis_file)<<"\n";
            }
            metis_file->close();
        
            rMeshReader.Reset();
    
            /*
             *  Call METIS binary to perform the partitioning.
             *  It will output a file called metis.mesh.npart.numProcs
             */
            std::stringstream permute_command;
            permute_command <<  "./bin/partdmesh "
                            <<  handler.GetOutputDirectoryFullPath("")
                            <<  basename << " "
                            <<  num_procs
                            <<  " > /dev/null";
    
            // METIS doesn't return 0 after a successful execution
            IGNORE_RET(system, permute_command.str());
        }
        catch (Exception &e)
        {
            PetscTools::ReplicateException(true);
            throw e;
        }
    }
    //Other processes need to see the exception happening
    PetscTools::ReplicateException(false);
     
     /*
     * Wait for the permutation to be available
     */
    PetscTools::Barrier();

    /*
     *  Read partition file and compute local node ownership and processors offset
     */
    std::ifstream partition_stream;
    std::string full_path = handler.GetOutputDirectoryFullPath("")
                            + output_file.str();

    partition_stream.open(full_path.c_str());
    assert(partition_stream.is_open());

    assert(rProcessorsOffset.size() == 0); // Making sure the vector is empty. After calling resize() only newly created memory will be initialised to 0.
    rProcessorsOffset.resize(PetscTools::NumProcs(), 0);

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
        for (unsigned proc=part_read+1; proc<PetscTools::NumProcs(); proc++)
        {
            rProcessorsOffset[proc]++;
        }
    }

    /*
     *  Once we know the offsets we can compute the permutation vector
     */
    // Move stream pointer to the beginning of the file
    partition_stream.seekg (0, std::ios::beg);
    
    std::vector<unsigned> local_index(this->GetNumNodes(), 0);    
    
    rNodePermutation.resize(this->GetNumNodes());
    
    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {   
        unsigned part_read;

        partition_stream >> part_read;

        rNodePermutation[node_index] = rProcessorsOffset[part_read] + local_index[part_read];
        
        local_index[part_read]++;
    }

    partition_stream.close();

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MetisLibraryNodePartitioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                                                                                  std::set<unsigned>& rNodesOwned, 
                                                                                  std::vector<unsigned>& rProcessorsOffset,
                                                                                  std::vector<unsigned>& rNodePermutation) const
{
    assert(PetscTools::NumProcs() > 1);
    
    assert( ELEMENT_DIM==2 || ELEMENT_DIM==3 ); // Metis works with triangles and tetras


    int ne = rMeshReader.GetNumElements(); 
    int nn = rMeshReader.GetNumNodes(); 
    idxtype* elmnts = new idxtype[ne * (ELEMENT_DIM+1)];
    assert(elmnts != NULL);   

    unsigned counter=0;    
    for(unsigned element_number = 0; element_number < mTotalNumElements; element_number++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();

        for(unsigned i=0; i<ELEMENT_DIM+1; i++)
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
    int nparts = PetscTools::NumProcs();
    int edgecut;
    idxtype* epart = new idxtype[ne];
    assert(epart != NULL); 
    idxtype* npart = new idxtype[nn];
    assert(epart != NULL);

    METIS_PartMeshNodal(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);//, wgetflag, vwgt);    

    assert(rProcessorsOffset.size() == 0); // Making sure the vector is empty. After calling resize() only newly created memory will be initialised to 0.
    rProcessorsOffset.resize(PetscTools::NumProcs(), 0);

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
        for (unsigned proc=part_read+1; proc<PetscTools::NumProcs(); proc++)
        {
            rProcessorsOffset[proc]++;
        }
    }

    /*
     *  Once we know the offsets we can compute the permutation vector
     */    
    std::vector<unsigned> local_index(this->GetNumNodes(), 0);    
    
    rNodePermutation.resize(this->GetNumNodes());
    
    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {   
        unsigned part_read = npart[node_index];

        rNodePermutation[node_index] = rProcessorsOffset[part_read] + local_index[part_read];
        
        local_index[part_read]++;
    }
    
    delete[] elmnts;
    delete[] epart;
    delete[] npart;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReorderNodes(std::vector<unsigned>& rNodePermutation)
{
    assert(PetscTools::NumProcs() > 1);
    
    // Need to rebuild global-local maps
    mNodesMapping.clear();
    mHaloNodesMapping.clear();
    
    //Update indices
    for (unsigned index=0; index<this->mNodes.size(); index++)
    {
        unsigned old_index = this->mNodes[index]->GetIndex();
        unsigned new_index = rNodePermutation[old_index];
        
        this->mNodes[index]->SetIndex(new_index);
        mNodesMapping[new_index] = index;
    }

    for (unsigned index=0; index<mHaloNodes.size(); index++)
    {
        unsigned old_index = mHaloNodes[index]->GetIndex();
        unsigned new_index = rNodePermutation[old_index];
        
        mHaloNodes[index]->SetIndex(new_index);
        mHaloNodesMapping[new_index] = index;
    }
}


#endif /*PARALLELTETRAHEDRALMESH_HPP_*/
