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
#ifndef MESHBASEDTISSUE_HPP_
#define MESHBASEDTISSUE_HPP_

#include "AbstractTissue.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "VoronoiTessellation.hpp"
#include "Exception.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>


/**
 * A facade class encapsulating a mesh-based 'tissue'
 * 
 * Contains a group of cells and maintains the associations between cells and
 * nodes in the mesh.
 * 
 */
template<unsigned DIM>
class MeshBasedTissue : public AbstractTissue<DIM>
{
protected:
    
    ConformingTetrahedralMesh<DIM, DIM>& mrMesh;
    
    VoronoiTessellation<DIM>* mpVoronoiTessellation;
    
    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this tissue has been de-serialized.
     */
    bool mDeleteMesh;
        
    /**
     * Special springs that we want to keep track of for some reason.
     * Currently used to track cells in the process of dividing
     * (which are represented as two cells joined by a shorter spring).
     */
    std::set<std::set<TissueCell*> > mMarkedSprings;
    
    /** Whether to print out cell area and perimeter info */
    bool mWriteVoronoiData;
    
    /** Whether to follow only the logged cell if writing voronoi data */
    bool mFollowLoggedCell;
    
    /** Whether to print out tissue areas */
    bool mWriteTissueAreas;
    
    /** Results file for elements */
    out_stream mpElementFile;
    
    /** Results file for Voronoi data */
    out_stream mpVoronoiFile;
    
    /** Results file for tissue area data */
    out_stream mpTissueAreasFile;
        
    /** Helper method used by the spring marking routines */
    std::set<TissueCell*> CreateCellPair(TissueCell&, TissueCell&);
        
    friend class boost::serialization::access;
    /**
     * Serialize the facade.
     * 
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     * 
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTissue<DIM> >(*this);
  
        // The Voronoi stuff can't be archived yet
        //archive & mpVoronoiTessellation
        delete mpVoronoiTessellation;
        mpVoronoiTessellation = NULL;
        
        archive & mMarkedSprings;
        archive & mWriteVoronoiData;
        archive & mFollowLoggedCell;
        archive & mWriteTissueAreas;
        
        // In its present form, a call to MeshBasedTissue::Validate() here
        // would result in a seg fault in the situation where we are actually loading 
        // a MeshBasedTissueWithGhostNodes. Commenting out this call breaks no tests.
           
        // Validate();
    }

public:

    /** Hack until meshes are fully archived using boost::serialization */
    static std::string meshPathname;

    /**
     * Create a new tissue facade from a mesh and collection of cells.
     * 
     * There must be precisely 1 cell for each node of the mesh.
     * 
     * @param rMesh a conforming tetrahedral mesh.
     * @param cells TissueCells corresponding to the nodes of the mesh.
     * @param deleteMesh set to true if you want the tissue to free the mesh memory on destruction
     */
    MeshBasedTissue(ConformingTetrahedralMesh<DIM, DIM>&, 
                    const std::vector<TissueCell>&,
                    bool deleteMesh=false, 
                    bool validate=true);

    /**
     * Constructor for use by the de-serializer.
     * 
     * @param rMesh a conforming tetrahedral mesh.
     */
    MeshBasedTissue(ConformingTetrahedralMesh<DIM, DIM>&);

    ~MeshBasedTissue();

    ConformingTetrahedralMesh<DIM, DIM>& rGetMesh();

    const ConformingTetrahedralMesh<DIM, DIM>& rGetMesh() const;

    bool GetWriteVoronoiData();

    bool GetWriteTissueAreas();

    void SetWriteVoronoiData(bool writeVoronoiData, bool followLoggedCell);

    void SetWriteTissueAreas(bool writeTissueAreas);

    /** 
     * Remove all cells labelled as dead. 
     * 
     * Note that this now calls 
     * ConformingTetrahedralMesh::DeleteNodePriorToReMesh() 
     * and therefore a ReMesh(map) must be called before
     * any element information is used.
     * 
     * Note also that after calling this method the tissue will be in an inconsistent state until a
     * ReMesh is performed!  So don't try iterating over cells or anything like that.
     * \todo weaken the data invariant in this class so it doesn't require an exact correspondance
     *  between nodes and cells.
     * 
     *  @return number of cells removed
     */
    unsigned RemoveDeadCells();

    void CreateOutputFiles(const std::string &rDirectory, 
                           bool rCleanOutputDirectory,
                           bool outputCellMutationStates,
                           bool outputCellTypes,
                           bool outputCellVariables,
                           bool outputCellCyclePhases,
                           bool outputCellAncestors);

    void CloseOutputFiles(bool outputCellMutationStates,
                          bool outputCellTypes,
                          bool outputCellVariables,
                          bool outputCellCyclePhases,
                          bool outputCellAncestors);

    /**
     * Move a cell to a new location.
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    void MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation);

    /**
     * Add a new cell to the tissue.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell*  AddCell(TissueCell cell, c_vector<double,DIM> newLocation);

    virtual void ReMesh();

    Node<DIM>* GetNode(unsigned index);

    unsigned GetNumNodes();

    /** 
     * Sets the Ancestor index of all the cells at the bottom in order,
     * can be used to trace clonal populations.
     */   
    void SetBottomCellAncestors();

    /**
     * Check consistency of our internal data structures. Each node must
     * have a cell associated with it.
     */
    virtual void Validate();

    void WriteResultsToFiles(bool outputCellMutationStates, 
                             bool outputCellTypes, 
                             bool outputCellVariables,
                             bool outputCellCyclePhases,
                             bool outputCellAncestors);

    void WriteVoronoiResultsToFile();

    void WriteTissueAreaResultsToFile();

    /** Get a reference to a Voronoi Tessellation of the mesh */                         
    void CreateVoronoiTessellation();

    VoronoiTessellation<DIM>& rGetVoronoiTessellation();
    
    /**
     * Update mIsGhostNode if required by a remesh.
     */ 
    virtual void UpdateGhostNodesAfterReMesh(NodeMap& rMap);

    /**
     * Iterator over edges in the mesh, which correspond to springs between cells.
     * 
     * This class takes care of the logic to make sure that you consider each edge exactly once.
     */
    class SpringIterator
    {
    public:
    
        /**
         * Get a pointer to the node in the mesh at end A of the spring.
         */
        Node<DIM>* GetNodeA();
        
        /**
         * Get a pointer to the node in the mesh at end B of the spring.
         */
        Node<DIM>* GetNodeB();
        
        /**
         * Get a *reference* to the cell at end A of the spring.
         */
        TissueCell& rGetCellA();
        
        /**
         * Get a *reference* to the cell at end B of the spring.
         */
        TissueCell& rGetCellB();
        
        bool operator!=(const SpringIterator& other);
        
        /**
         * Prefix increment operator.
         */
        SpringIterator& operator++();
        
        /**
         * Constructor for a new iterator.
         */
        SpringIterator(MeshBasedTissue& rTissue, typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edgeIter);
        
    private:
    
        /** Keep track of what edges have been visited */
        std::set<std::set<unsigned> > mSpringsVisited;
    
        MeshBasedTissue& mrTissue;
        
        typename ConformingTetrahedralMesh<DIM, DIM>::EdgeIterator mEdgeIter;
    };

    /**
     * @return iterator pointing to the first spring in the tissue
     */
    SpringIterator SpringsBegin();
    
    /**
     * @return iterator pointing to one past the last spring in the tissue
     */
    SpringIterator SpringsEnd();
    
    // For debugging
    void CheckTissueCellPointers();
    
    /**
     * Test whether the spring between 2 cells is marked.
     */
    bool IsMarkedSpring(TissueCell&, TissueCell&);
    
    /**
     * Mark the spring between the given cells.
     */
    void MarkSpring(TissueCell&, TissueCell&);
    
    /**
     * Stop marking the spring between the given cells.
     */
    void UnmarkSpring(TissueCell&, TissueCell&);

};

template<unsigned DIM>
MeshBasedTissue<DIM>::MeshBasedTissue(ConformingTetrahedralMesh<DIM, DIM>& rMesh,
                  const std::vector<TissueCell>& rCells,
                  bool deleteMesh,
                  bool validate)
             : AbstractTissue<DIM>(rCells),
               mrMesh(rMesh),
               mpVoronoiTessellation(NULL),
               mDeleteMesh(deleteMesh),
               mWriteVoronoiData(false),
               mFollowLoggedCell(false),
               mWriteTissueAreas(false)
{
    // This must always be true
    assert( this->mCells.size() <= mrMesh.GetNumNodes() );

    this->mTissueContainsMesh = true;
    
    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
MeshBasedTissue<DIM>::MeshBasedTissue(ConformingTetrahedralMesh<DIM, DIM>& rMesh)
             : mrMesh(rMesh)
{
    this->mTissueContainsMesh = true;
    mpVoronoiTessellation = NULL;
    mDeleteMesh = true;
}

template<unsigned DIM>
MeshBasedTissue<DIM>::~MeshBasedTissue()
{
    delete mpVoronoiTessellation;
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::Validate()
{
    std::vector<bool> validated_node = std::vector<bool>(this->GetNumNodes(), false);
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        validated_node[node_index] = true;
    }    
    
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str()); 
        }
    }
}

template<unsigned DIM>
ConformingTetrahedralMesh<DIM, DIM>& MeshBasedTissue<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
const ConformingTetrahedralMesh<DIM, DIM>& MeshBasedTissue<DIM>::rGetMesh() const
{
    return mrMesh;
}

template<unsigned DIM>
unsigned MeshBasedTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if (it->IsDead())
        {
            // Check if this cell is in a marked spring
            std::vector<const std::set<TissueCell*>*> pairs_to_remove; // Pairs that must be purged
            for (std::set<std::set<TissueCell*> >::iterator it1 = mMarkedSprings.begin();
                 it1 != mMarkedSprings.end();
                 ++it1)
            {
                const std::set<TissueCell*>& r_pair = *it1;
                for (std::set<TissueCell*>::iterator it2 = r_pair.begin();
                     it2 != r_pair.end();
                     ++it2)
                {
                    TissueCell* p_cell = *it2;
                    if (p_cell == &(*it))
                    {
                        // Remember to purge this spring
                        pairs_to_remove.push_back(&r_pair);
                        break;
                    }
                }
            }
            // Purge any marked springs that contained this cell
            for (std::vector<const std::set<TissueCell*>* >::iterator pair_it = pairs_to_remove.begin();
                 pair_it != pairs_to_remove.end();
                 ++pair_it)
            {
                mMarkedSprings.erase(**pair_it);
            }
            
            // Remove the node from the mesh
            num_removed++;
            mrMesh.DeleteNodePriorToReMesh(it->GetNodeIndex());
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    mrMesh.SetNode(index, rNewLocation, false);
}

template<unsigned DIM>  
TissueCell* MeshBasedTissue<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    Node<DIM>* p_new_node = new Node<DIM>(mrMesh.GetNumNodes(), newLocation, false);   // never on boundary
              
    unsigned new_node_index = mrMesh.AddNode(p_new_node);

    newCell.SetNodeIndex(new_node_index);
    this->mCells.push_back(newCell);
    
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mNodeCellMap[new_node_index] = p_created_cell;
    
    return p_created_cell;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::ReMesh()
{    
    NodeMap map(mrMesh.GetNumAllNodes());
    if (DIM==2)
    {
        mrMesh.ReMeshWithTriangleLibrary(map);
    }
    else
    {
        mrMesh.ReMesh(map);
    }

    if (!map.IsIdentityMap())
    {
        UpdateGhostNodesAfterReMesh(map);
        
        // Fix up the mappings between cells and nodes
        this->mNodeCellMap.clear();
        for (std::list<TissueCell>::iterator it = this->mCells.begin();
             it != this->mCells.end();
             ++it)
        {
            unsigned old_node_index = it->GetNodeIndex();
            
            // This shouldn't ever happen, as the cell vector only contains living cells
            assert(!map.IsDeleted(old_node_index));
            unsigned new_node_index = map.GetNewIndex(old_node_index);
            it->SetNodeIndex(new_node_index);
            this->mNodeCellMap[new_node_index] = &(*it);
        }
    }
    
    // Purge any marked springs that are no longer springs
    std::vector<const std::set<TissueCell*>*> springs_to_remove;
    for (std::set<std::set<TissueCell*> >::iterator spring_it = mMarkedSprings.begin();
         spring_it != mMarkedSprings.end();
         ++spring_it)
    {
        const std::set<TissueCell*>& r_pair = *spring_it;
        assert(r_pair.size() == 2);
        TissueCell* p_cell_1 = *(r_pair.begin());
        TissueCell* p_cell_2 = *(++r_pair.begin());
        Node<DIM>* p_node_1 = this->GetNodeCorrespondingToCell(*p_cell_1);
        Node<DIM>* p_node_2 = this->GetNodeCorrespondingToCell(*p_cell_2);
        
        bool joined = false;
        
        // For each element containing node1, if it also contains node2 then the cells are joined
        std::set<unsigned> node2_elements = p_node_2->rGetContainingElementIndices();
        for (typename Node<DIM>::ContainingElementIterator elt_it = p_node_1->ContainingElementsBegin();
             elt_it != p_node_1->ContainingElementsEnd();
             ++elt_it)
        {
            unsigned elt_index = *elt_it;
            if (node2_elements.find(elt_index) != node2_elements.end())
            {
                joined = true;
                break;
            }
        }
        
        // If no longer joined, remove this spring from the set
        if (!joined)
        {
            springs_to_remove.push_back(&r_pair);
        }
    }
    
    // Remove any springs necessary
    for (std::vector<const std::set<TissueCell*>* >::iterator spring_it = springs_to_remove.begin();
         spring_it != springs_to_remove.end();
         ++spring_it)
    {
        mMarkedSprings.erase(**spring_it);
    }
    
    Validate();
}

template<unsigned DIM>
Node<DIM>* MeshBasedTissue<DIM>::GetNode(unsigned index)
{
    return rGetMesh().GetNode(index);
}

template<unsigned DIM>
unsigned MeshBasedTissue<DIM>::GetNumNodes()
{
    return rGetMesh().GetNumAllNodes();
}

template<unsigned DIM> 
void MeshBasedTissue<DIM>::SetBottomCellAncestors()
{
    unsigned index = 0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        if (cell_iter.rGetLocation()[1] < 0.5)
        {
            cell_iter->SetAncestor(index++);
        }
    }
}

template<unsigned DIM> 
void MeshBasedTissue<DIM>::SetWriteVoronoiData(bool writeVoronoiData, bool followLoggedCell)
{
    assert(DIM == 2);
    mWriteVoronoiData = writeVoronoiData;
    mFollowLoggedCell = followLoggedCell;
}

template<unsigned DIM> 
void MeshBasedTissue<DIM>::SetWriteTissueAreas(bool writeTissueAreas)
{
    assert(DIM == 2);
    mWriteTissueAreas = writeTissueAreas;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{    
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void MeshBasedTissue<DIM>::CreateOutputFiles(const std::string &rDirectory, 
                                             bool rCleanOutputDirectory, 
                                             bool outputCellMutationStates,
                                             bool outputCellTypes,
                                             bool outputCellVariables,
                                             bool outputCellCyclePhases,
                                             bool outputCellAncestors)
{
    AbstractTissue<DIM>::CreateOutputFiles(rDirectory, 
                                           rCleanOutputDirectory, 
                                           outputCellMutationStates,
                                           outputCellTypes,
                                           outputCellVariables,
                                           outputCellCyclePhases,
                                           outputCellAncestors);
                                           
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpElementFile = output_file_handler.OpenOutputFile("results.vizelements");
    
    if (mWriteVoronoiData)
    {
        mpVoronoiFile = output_file_handler.OpenOutputFile("results.vizvoronoi");
    }
    if (mWriteTissueAreas)
    {
        mpTissueAreasFile = output_file_handler.OpenOutputFile("tissueareas.dat");
    }
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::CloseOutputFiles(bool outputCellMutationStates,
                                            bool outputCellTypes,
                                            bool outputCellVariables,
                                            bool outputCellCyclePhases,
                                            bool outputCellAncestors)
{
    AbstractTissue<DIM>::CloseOutputFiles(outputCellMutationStates,
                                          outputCellTypes,
                                          outputCellVariables,
                                          outputCellCyclePhases,
                                          outputCellAncestors);
    mpElementFile->close();
    
    if (mWriteVoronoiData)
    {
        mpVoronoiFile->close();
    }
    if (mWriteTissueAreas)
    {
        mpTissueAreasFile->close();
    }
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::GetWriteVoronoiData()
{
    return mWriteVoronoiData;        
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::GetWriteTissueAreas()
{
    return mWriteTissueAreas;        
}

template<unsigned DIM>  
void MeshBasedTissue<DIM>::WriteResultsToFiles(bool outputCellMutationStates, 
                                               bool outputCellTypes, 
                                               bool outputCellVariables,
                                               bool outputCellCyclePhases,
                                               bool outputCellAncestors)
{
    AbstractTissue<DIM>::WriteResultsToFiles(outputCellMutationStates, 
                                             outputCellTypes, 
                                             outputCellVariables,
                                             outputCellCyclePhases,
                                             outputCellAncestors);
    
    // Write element data to file
    
    *mpElementFile <<  SimulationTime::Instance()->GetDimensionalisedTime() << "\t";
    
    for (unsigned elem_index=0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
    {
        if (!mrMesh.GetElement(elem_index)->IsDeleted())
        {
            for(unsigned i=0; i<DIM+1; i++)
            {
                *mpElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i)<< " ";
            }
        }
    }
    
    *mpElementFile << "\n";    
    
    if (mpVoronoiTessellation!=NULL)
    {
        // Write Voronoi data to file if required    
        if (mWriteVoronoiData)
        {
            WriteVoronoiResultsToFile();
        }
        
        // Write tissue area data to file if required
        if (mWriteTissueAreas)
        {
            WriteTissueAreaResultsToFile();
        }
    }
}

template<unsigned DIM>  
void MeshBasedTissue<DIM>::WriteVoronoiResultsToFile()
{
    // Write time to file
    *mpVoronoiFile << SimulationTime::Instance()->GetDimensionalisedTime() << " ";
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        if ((!mFollowLoggedCell) || ((mFollowLoggedCell) && (cell_iter->IsLogged())))
        {
            unsigned node_index = cell_iter.GetNode()->GetIndex();
            double x = cell_iter.rGetLocation()[0];
            double y = cell_iter.rGetLocation()[1];
        
            double cell_area = rGetVoronoiTessellation().GetFaceArea(node_index);
            double cell_perimeter = rGetVoronoiTessellation().GetFacePerimeter(node_index);
        
            *mpVoronoiFile << node_index << " " << x << " " << y << " " << cell_area << " " << cell_perimeter << " ";
            
            if (mFollowLoggedCell)
            {
                break;
            }
        }
    }
    *mpVoronoiFile << "\n";
}

template<unsigned DIM>  
void MeshBasedTissue<DIM>::WriteTissueAreaResultsToFile()
{
    // Write time to file
    *mpTissueAreasFile << SimulationTime::Instance()->GetDimensionalisedTime() << " ";
        
    // Don't use the Voronoi tessellation to calculate the total area
    // because it gives huge areas for boundary cells
    double total_area = rGetMesh().CalculateVolume();    
        
    double necrotic_area = 0.0;
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Only bother calculating the cell area if it is necrotic
        if (cell_iter->GetCellType() == NECROTIC)
        {
            unsigned node_index = cell_iter.GetNode()->GetIndex();                
            double cell_area = rGetVoronoiTessellation().GetFace(node_index)->GetArea();
            necrotic_area += cell_area;
        }
    }       
    
    *mpTissueAreasFile << total_area << " " << necrotic_area << "\n";
}

//////////////////////////////////////////////////////////////////////////////
//                          Spring iterator class                           // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
Node<DIM>* MeshBasedTissue<DIM>::SpringIterator::GetNodeA()
{
    return mEdgeIter.GetNodeA();
}

template<unsigned DIM>
Node<DIM>* MeshBasedTissue<DIM>::SpringIterator::GetNodeB()
{
    return mEdgeIter.GetNodeB();
}

template<unsigned DIM>
TissueCell& MeshBasedTissue<DIM>::SpringIterator::rGetCellA()
{
    assert((*this) != mrTissue.SpringsEnd());
    return mrTissue.rGetCellAtNodeIndex(mEdgeIter.GetNodeA()->GetIndex());
}

template<unsigned DIM>
TissueCell& MeshBasedTissue<DIM>::SpringIterator::rGetCellB()
{
    assert((*this) != mrTissue.SpringsEnd());
    return mrTissue.rGetCellAtNodeIndex(mEdgeIter.GetNodeB()->GetIndex());
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::SpringIterator::operator!=(const MeshBasedTissue<DIM>::SpringIterator& other)
{
    return (mEdgeIter != other.mEdgeIter);
}

template<unsigned DIM>
typename MeshBasedTissue<DIM>::SpringIterator& MeshBasedTissue<DIM>::SpringIterator::operator++()
{
    bool edge_is_ghost = false;
    
    do
    {
        ++mEdgeIter;
        if (*this != mrTissue.SpringsEnd())
        {
            bool a_is_ghost = mrTissue.IsGhostNode(mEdgeIter.GetNodeA()->GetIndex());
            bool b_is_ghost = mrTissue.IsGhostNode(mEdgeIter.GetNodeB()->GetIndex());

            edge_is_ghost = (a_is_ghost || b_is_ghost);
        }
    }
    while (*this!=mrTissue.SpringsEnd() && edge_is_ghost); 

    return (*this);
}

template<unsigned DIM>
MeshBasedTissue<DIM>::SpringIterator::SpringIterator(MeshBasedTissue& rTissue,
                                           typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edgeIter)
    : mrTissue(rTissue),
      mEdgeIter(edgeIter)
{
    if (mEdgeIter!=mrTissue.mrMesh.EdgesEnd())
    {
        bool a_is_ghost = mrTissue.IsGhostNode(mEdgeIter.GetNodeA()->GetIndex());
        bool b_is_ghost = mrTissue.IsGhostNode(mEdgeIter.GetNodeB()->GetIndex());

        if (a_is_ghost || b_is_ghost)
        {
            ++(*this);
        }
    }
}

template<unsigned DIM>
typename MeshBasedTissue<DIM>::SpringIterator MeshBasedTissue<DIM>::SpringsBegin()
{
    return SpringIterator(*this, mrMesh.EdgesBegin());
}

template<unsigned DIM>
typename MeshBasedTissue<DIM>::SpringIterator MeshBasedTissue<DIM>::SpringsEnd()
{
    return SpringIterator(*this, mrMesh.EdgesEnd());
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::CreateVoronoiTessellation()
{
    delete mpVoronoiTessellation;
    mpVoronoiTessellation = new VoronoiTessellation<DIM>(mrMesh);
}

template<unsigned DIM>
VoronoiTessellation<DIM>& MeshBasedTissue<DIM>::rGetVoronoiTessellation()
{
    assert(mpVoronoiTessellation!=NULL);
    return *mpVoronoiTessellation;
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::CheckTissueCellPointers()
{
    bool res = true;
    for (std::list<TissueCell>::iterator it=this->mCells.begin();
         it!=this->mCells.end();
         ++it)
    {
        TissueCell* p_cell=&(*it);
        assert(p_cell);
        AbstractCellCycleModel *p_model = p_cell->GetCellCycleModel();
        assert(p_model);
        
        // Check cell exists in tissue
        unsigned node_index = p_cell->GetNodeIndex();
        std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
        TissueCell& r_cell = this->rGetCellAtNodeIndex(node_index);
#define COVERAGE_IGNORE //Debugging code.  Shouldn't fail under normal conditions 
        if (&r_cell != p_cell)
        {
            std::cout << "  Mismatch with tissue" << std::endl << std::flush;
            res = false;
        }
        
        // Check model links back to cell
        if (p_model->GetCell() != p_cell)
        {
            std::cout << "  Mismatch with cycle model" << std::endl << std::flush;
            res = false;
        }
    }
    assert(res);
#undef COVERAGE_IGNORE 
    
    res = true;
    for (std::set<std::set<TissueCell*> >::iterator it1 = mMarkedSprings.begin();
         it1 != mMarkedSprings.end();
         ++it1)
    {
        const std::set<TissueCell*>& r_pair = *it1;
        assert(r_pair.size() == 2);
        for (std::set<TissueCell*>::iterator it2 = r_pair.begin();
             it2 != r_pair.end();
             ++it2)
        {
            TissueCell* p_cell = *it2;
            assert(p_cell);
            AbstractCellCycleModel *p_model = p_cell->GetCellCycleModel();
            assert(p_model);
            unsigned node_index = p_cell->GetNodeIndex();
            std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
            
#define COVERAGE_IGNORE //Debugging code.  Shouldn't fail under normal conditions 
            // Check cell is alive
            if (p_cell->IsDead())
            {
                std::cout << "  Cell is dead" << std::endl << std::flush;
                res = false;
            }
            
            // Check cell exists in tissue
            TissueCell& r_cell = this->rGetCellAtNodeIndex(node_index);
            if (&r_cell != p_cell)
            {
                std::cout << "  Mismatch with tissue" << std::endl << std::flush;
                res = false;
            }
            
            // Check model links back to cell
            if (p_model->GetCell() != p_cell)
            {
                std::cout << "  Mismatch with cycle model" << std::endl << std::flush;
                res = false;
            }
        }
#undef COVERAGE_IGNORE
    }
    assert(res);
}

template<unsigned DIM>
std::set<TissueCell*> MeshBasedTissue<DIM>::CreateCellPair(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair;
    cell_pair.insert(&rCell1);
    cell_pair.insert(&rCell2);
    return cell_pair;
}

template<unsigned DIM>
bool MeshBasedTissue<DIM>::IsMarkedSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    return mMarkedSprings.find(cell_pair) != mMarkedSprings.end();
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::MarkSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    mMarkedSprings.insert(cell_pair);
}

template<unsigned DIM>
void MeshBasedTissue<DIM>::UnmarkSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    mMarkedSprings.erase(cell_pair);
}


template<unsigned DIM>
std::string MeshBasedTissue<DIM>::meshPathname = "";

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedTissue)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Tissue facade.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MeshBasedTissue<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const ConformingTetrahedralMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MeshBasedTissue<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    assert(MeshBasedTissue<DIM>::meshPathname.length() > 0);
    ConformingTetrahedralMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;
    
    // Re-initialise the mesh
    p_mesh->Clear();
    TrianglesMeshReader<DIM,DIM> mesh_reader(MeshBasedTissue<DIM>::meshPathname);
    p_mesh->ConstructFromMeshReader(mesh_reader);
    
    // Needed for cylindrical meshes at present; should be safe in any case.
    NodeMap map(p_mesh->GetNumNodes());
    if (DIM==2u)
    {
        p_mesh->ReMeshWithTriangleLibrary(map);
    }
    else
    {
        p_mesh->ReMesh(map);
    }
    
    // Invoke inplace constructor to initialize instance
    ::new(t)MeshBasedTissue<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*MESHBASEDTISSUE_HPP_*/
