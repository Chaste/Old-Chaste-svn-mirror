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
#ifndef TISSUECELL_HPP_
#define TISSUECELL_HPP_

#include <boost/serialization/access.hpp>

#include "Element.hpp"
#include "CellTypes.hpp"
#include "CellMutationStates.hpp"
#include "AbstractCellCycleModel.hpp"
#include "SimulationTime.hpp"

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).


/**
 * Tissue cell is the basic container for all the biological information about a cell.
 * It contains the cell cycle model and all other biological properties such as mutation
 * state, cell type, whether it is undergoing apoptosis or not.
 *
 * This class should not store any spatial information - TissueCells are linked to space by the Tissue classes.
 */
class TissueCell
{
private:

    /** Caches the result of ReadyToDivide() so Divide() can look at it. */
    bool mCanDivide;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        // These first four are also dealt with by {load,save}_construct_data
        archive & mCanDivide;
        archive & mCellType;
        archive & mMutationState;
        archive & mpCellCycleModel;
        archive & mUndergoingApoptosis;
        archive & mDeathTime;
        archive & mStartOfApoptosisTime;
        archive & mIsDead;
        archive & mIsLogged;
        archive & mAncestor;
        archive & mCellId;
        archive & mMaxCellId;
    }

protected:

    // NB - if you add any member variables, make sure CommonCopy includes them.

    /** The cell type - defined in CellTypes.hpp */
    CellType mCellType;

    /** The cell's mutation state - defined in CellMutationStates.hpp */
    CellMutationState mMutationState;

    /** The cell's cell-cycle model */
    AbstractCellCycleModel *mpCellCycleModel;

    /** An index which is inherited by all children of this cell */
    unsigned mAncestor;

    /** An identifier which is unique to this cell */
    unsigned mCellId;

    /** maximum cell identifier */
    static unsigned mMaxCellId;

    /** When the cell will/did die */
    double mDeathTime;

    /** When the cell was commanded to start apoptosis */
    double mStartOfApoptosisTime;

    /** Whether the cell is currently in apoptosis - don't divide */
    bool mUndergoingApoptosis;

    /**
     * Whether the cell is dead or not (they exist in the Tissue until they are
     * removed by AbstractTissue::RemoveDeadCells()
     */
    bool mIsDead;

    /** Whether the cell is being tracked specially. */
    bool mIsLogged;

    /**
     * Contains code common to both the copy constructor and operator=.
     *
     * @param rOtherCell  An existing TissueCell
     */
    void CommonCopy(const TissueCell& rOtherCell);

public:

    /**
     * Create a new tissue cell.
     * @param cellType  the type of cell this is
     * @param mutationState the mutation state of the cell
     * @param pCellCycleModel  the cell cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     * @param archiving  whether this constructor is being called by the archiver - do things slightly differently!
     */
    TissueCell(CellType cellType,
               CellMutationState mutationState,
               AbstractCellCycleModel* pCellCycleModel,
               bool archiving = false);

    /**
     * Destructor, which frees the memory allocated for our cell cycle model.
     */
    ~TissueCell();

    /**
     * Create a new tissue cell that is a copy of an existing cell
     * @param rOtherCell  An existing TissueCell
     */
    TissueCell(const TissueCell& rOtherCell);

    /**
     * Copy all the attributes of one cell to another.
     *
     * @param rOtherCell  An existing TissueCell
     */
    TissueCell& operator=(const TissueCell& rOtherCell);

    /**
     * Set the birth time of the cell - can be negative so that your cells have an age when a simulation begins
     *
     * @param birthTime  The time the cell was born (in hours)
     */
    void SetBirthTime(double birthTime);

    /**
     * Change the cell cycle model used. This takes effect immediately.
     *
     * @param pCellCycleModel pointer to the cell cycle model to use
     */
    void SetCellCycleModel(AbstractCellCycleModel* pCellCycleModel);

    /**
     * Returns a pointer to the TissueCell's cell cycle model.
     */
    AbstractCellCycleModel* GetCellCycleModel() const;

    /**
     * Calls Initialise on the cell cycle model associated with this cell.
     */
    void InitialiseCellCycleModel();

    /**
     * Get the cell's age from its cell cycle model.
     */
    double GetAge() const;

    /**
     * Get the cell's birth time from its cell cycle model.
     */
    double GetBirthTime() const;

    /**
     * Get the time at which apoptosis was commanded to start.
     */
    double GetStartOfApoptosisTime() const;

    /**
     * Get method for mCellType.
     */
    CellType GetCellType() const;

    /**
     * Set method for mCellType.
     *
     * @param cellType the cell's type
     */
    void SetCellType(CellType cellType);

    /**
     * Get method for mMutationState.
     */
    CellMutationState GetMutationState() const;

    /**
     * Set method for mMutationState.
     *
     * @param mutationState the cell's mutation state
     */
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

    /**
     * Makes the cell enter apoptosis and sets #mDeathTime using the apoptosis
     * time from the cancer parameters.
     *
     * @param setDeathTime whether we tell the cell exactly when to die (defaults to true)
     */
    void StartApoptosis(bool setDeathTime=true);

    /**
     * This labels the cell as dead, it does not delete the cell, it remains
     * in the Tissue until AbstractTissue::RemoveDeadCells() is called.
     */
    void Kill();

    /**
     * Returns whether the cell is undergoing apoptosis or not.
     */
    bool HasApoptosisBegun() const;

    /**
     * @return How long until the cell dies (if it is in apoptosis, throws an exception if not)
     */
    double TimeUntilDeath() const;

    /**
     * Return whether the cell is dead or undergoing apoptosis.
     */
    bool IsDead();

    /**
     * Sets a flag to perform special output on this cell only.
     */
    void SetLogged();

    /**
     * @return Whether the cell is being tracked.
     */
    bool IsLogged();

    /**
     * Give the TissueCell an index which it passes to its children.
     *
     * @param ancestorIndex the cell's ancestor index
     */
    void SetAncestor(unsigned ancestorIndex);

    /**
     * @return The ancestor index, inherited from parents or set using the method above,
     * used for monoclonality experiments.
     */
    unsigned GetAncestor() const;

    /**
     * @return The cell identifier.
     */
    unsigned GetCellId() const;

    /**
     * Reset the current max ID.
     */
    static void ResetMaxCellId();
};


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueCell.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const TissueCell * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const CellType cell_type = t->GetCellType();
    const CellMutationState mutation_state = t->GetMutationState();
    const AbstractCellCycleModel* const p_cell_cycle_model = t->GetCellCycleModel();
    ar << cell_type;
    ar << mutation_state;
    ar << p_cell_cycle_model;
}

/**
 * De-serialize constructor parameters and initialise a TissueCell.
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

    // Invoke inplace constructor to initialise instance
    ::new(t)TissueCell(cell_type, mutation_state, p_cell_cycle_model, archiving);
}
}
} // namespace ...

#endif /*TISSUECELL_HPP_*/
