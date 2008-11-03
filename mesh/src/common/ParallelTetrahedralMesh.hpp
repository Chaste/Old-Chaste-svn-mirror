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

#ifndef PARALLELTETRAHEDRALMESH_HPP_
#define PARALLELTETRAHEDRALMESH_HPP_

#include "TetrahedralMesh.hpp"
#include "Element.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ParallelTetrahedralMesh : public TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
    
private: 

    unsigned mTotalNumElements; 
    unsigned mTotalNumNodes;
    
    std::map< unsigned, Element<ELEMENT_DIM, SPACE_DIM>* > mElements;
    //std::vector< Node<SPACE_DIM>* > mNodes; //defined in base class
    std::map< unsigned, BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* > mBoundaryElements;

    std::map< unsigned,Node<SPACE_DIM> *> mBoundaryNodes;
    

public:

//    ParallelTetrahedralMesh();

    virtual ~ParallelTetrahedralMesh();

    void ComputeMeshPartioning(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
                               std::vector<unsigned> &rElementsPartition);    

    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                 bool cullInternalFaces=false);

    // There are no local nodes yet
    //unsigned GetNumLocalNodes();

    unsigned GetNumLocalElements();

    unsigned GetNumNodes();

    unsigned GetNumElements();

    
};

//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ParallelTetrahedralMesh()
//{
//}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~ParallelTetrahedralMesh()
{
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ComputeMeshPartioning(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    std::vector<unsigned> &rElementsPartition)
{
    ///\todo: add a timing event for the partitioning
    
    // Not calling METIS yet, dumb partitioning (for the moment)
    unsigned num_procs = PetscTools::NumProcs();
    
    unsigned offset = 0;
    for (unsigned proc_index=0; proc_index<num_procs; proc_index++)
    {
        unsigned elements_assigned = mTotalNumElements/num_procs;
        if (proc_index < mTotalNumElements%num_procs)
        {
            elements_assigned++;
        }
                
        for(unsigned element_index=0; element_index<elements_assigned; element_index++)
        {
            rElementsPartition[offset++] = proc_index;            
        }        
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> &rMeshReader,
    bool cullInternalFaces)
{
    mTotalNumElements = rMeshReader.GetNumElements();
    mTotalNumNodes = rMeshReader.GetNumNodes();

    unsigned my_rank = PetscTools::GetMyRank();    
    
    std::vector<unsigned> elements_partition(mTotalNumElements);
    
    ComputeMeshPartioning(rMeshReader, elements_partition);

    // Load all the nodes (for the moment)
    std::vector<double> coords;
    for (unsigned node_index=0; node_index < mTotalNumNodes; node_index++)
    {
        coords = rMeshReader.GetNextNode();
        this->mNodes.push_back(new Node<SPACE_DIM>(node_index, coords, false));
    }

    // Load the elements owned by the processor
    std::vector<unsigned> node_indices;
    for (unsigned element_index=0; element_index < mTotalNumElements; element_index++)
    {
        node_indices = rMeshReader.GetNextElement();

        if (elements_partition[element_index] == my_rank)
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                nodes.push_back(this->mNodes[node_indices[j]]);
            }
    
            mElements[element_index] = new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes);                        
        }
    }
    
}


//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalNodes()
//{
//    return mNodes.size();
//}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumLocalElements()
{
    return mElements.size();
    //return mElementsPerProc[PetscTools::GetMyRank()];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mTotalNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mTotalNumElements;
}



#endif /*PARALLELTETRAHEDRALMESH_HPP_*/
