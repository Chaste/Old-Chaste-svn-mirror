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
#ifndef TISSUECELL_HPP_
#define TISSUECELL_HPP_

#include <boost/serialization/access.hpp>

#include "Element.hpp"
#include "CellTypes.hpp"
#include "CellMutationStates.hpp"
#include "AbstractCellCycleModel.hpp"
#include "SimulationTime.hpp"
#include "CancerParameters.hpp"

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).

class TissueCell
{
private:
    /// Caches the result of ReadyToDivide() so Divide() can look at it
    bool mCanDivide;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        // These first 4 are also dealt with by {load,save}_construct_data
        archive & mCanDivide;
        archive & mCellType;
        archive & mMutationState;
        archive & mpCellCycleModel;
        archive & mNodeIndex;
        archive & mUndergoingApoptosis;
        archive & mDeathTime;
        archive & mIsDead;
        archive & mIsLogged;
        archive & mAncestor;
    }
    
protected:
    // NB - if you add any member variables make sure CommonCopy includes them.
    CellType mCellType;
    CellMutationState mMutationState;
    AbstractCellCycleModel *mpCellCycleModel;
    unsigned mNodeIndex;
    bool mUndergoingApoptosis;
    /// When the cell will/did die.
    double mDeathTime;
    bool mIsDead;
    /// Whether the cell is being tracked specially.
    bool mIsLogged;
    /** An index which is inherited by all children of this cell */
    unsigned mAncestor;

    /**
     * Contains code common to both the copy constructor and operator=.
     */
    void CommonCopy(const TissueCell &other_cell);
    
public:

    /**
     * Create a new tissue cell.
     * @param cellType  the type of cell this is
     * @param mutationState the mutation state of the cell
     * @param generation  its generation
     * @param pCellCycleModel  the cell cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     * @param archiving  whether this constructor is being called by the archiver - do things slightly differently!
     */       
    TissueCell(CellType cellType,
               CellMutationState mutationState,
               AbstractCellCycleModel *pCellCycleModel,
               bool archiving = false);
               
    /**
     * Destructor, which frees the memory allocated for our cell cycle model.
     */
    ~TissueCell();
    
    TissueCell(const TissueCell &other_cell);
    
    /**
     * Copy all the attributes of one cell to another
     * (used for periodic boundaries - does not copy node or position information)
     *
     * \todo Is this the way we should be doing it?  Seems a bit counter-intuitive if the comment above is correct!
     *
     * \todo Also, since cell cycle models don't have an operator=, only copies data members of AbstractCellCycleModel when the model is copied.
     */
    TissueCell& operator=(const TissueCell &other_cell);
    
    void SetBirthTime(double birthTime);
    
    /**
     * Change the cell cycle model used.  This takes effect immediately.
     */
    void SetCellCycleModel(AbstractCellCycleModel *pCellCycleModel);
    AbstractCellCycleModel *GetCellCycleModel() const;
    
    /**
     * Calls Initialise on the cell cycle model associated with this cell.
     */
    void InitialiseCellCycleModel();
    
    // This method doesn't appear to actually be defined anywhere...
    // void UpdateCellCycleModel();
    
    /**
     * Set the node at which this cell is positioned.
     * 
     * @param index Index of the node in the mesh
     */
    void SetNodeIndex(unsigned index);
    unsigned GetNodeIndex() const;
    
    double GetAge() const;
    double GetBirthTime() const;
    unsigned GetGeneration() const;
    
    CellType GetCellType() const;
    void SetCellType(CellType cellType);
    CellMutationState GetMutationState() const;
    void SetMutationState(CellMutationState mutationState);
    
    /**
     * Determine if this cell is ready to divide at the current simulation time.
     * MUST be called before calling Divide().
     */
    bool ReadyToDivide();
    
    /**
     * Divide this cell to produce a daughter cell.
     * ReadyToDivide MUST have been called at the current time, and returned true.
     *
     * @return the new daughter cell
     */
    TissueCell Divide();
    
    void StartApoptosis();
    void Kill();
    bool HasApoptosisBegun() const;
    double TimeUntilDeath() const;
    bool IsDead() const;
    
    void SetLogged();
    bool IsLogged();
    
    void SetAncestor(unsigned ancestorIndex);
    unsigned GetAncestor() const;
};


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Meineke cell.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const TissueCell * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const CellType cell_type = t->GetCellType();
    const CellMutationState mutation_state = t->GetMutationState();
    const AbstractCellCycleModel * const p_cell_cycle_model = t->GetCellCycleModel();
    ar << cell_type;
    ar << mutation_state;
    ar << p_cell_cycle_model;
}

/**
 * De-serialize constructor parameters and initialise Meineke cell.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, TissueCell * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    CellType cell_type;
    CellMutationState mutation_state;
    AbstractCellCycleModel *p_cell_cycle_model;
    ar >> cell_type;
    ar >> mutation_state;
    ar >> p_cell_cycle_model;
    bool archiving = true;
    // Invoke inplace constructor to initialize instance
    ::new(t)TissueCell(cell_type, mutation_state,p_cell_cycle_model,archiving);
}
}
} // namespace ...

#endif /*TISSUECELL_HPP_*/
