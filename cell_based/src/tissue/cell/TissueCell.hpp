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
#ifndef TISSUECELL_HPP_
#define TISSUECELL_HPP_

#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include "CellProliferativeTypes.hpp"
#include "AbstractCellMutationState.hpp"
#include "AbstractCellCycleModel.hpp"
#include "SimulationTime.hpp"
#include "CellMutationStateRegistry.hpp"
#include "CellPropertyCollection.hpp"


class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).

class TissueCell;

/** Cells shouldn't be copied - it doesn't make sense.  So all access is via this pointer type. */
typedef boost::shared_ptr<TissueCell> TissueCellPtr;

/**
 * TissueCell is the basic container for all the biological information about a cell.
 * It contains the cell cycle model and all other biological properties such as mutation
 * state, cell type, whether it is undergoing apoptosis or not.
 *
 * This class should not store any spatial information - tissue cells are linked to space by the AbstractTissue subclasses.
 */
class TissueCell : boost::noncopyable, public boost::enable_shared_from_this<TissueCell> 
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
        archive & mpMutationState;
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
    /**
     * The cell's property collection.
     */
    CellPropertyCollection mCellPropertyCollection;

    /** The cell's mutation state. */
    boost::shared_ptr<AbstractCellMutationState> mpMutationState;

    /** The cell's cell-cycle model */
    AbstractCellCycleModel* mpCellCycleModel;

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

public:

    /**
     * Create a new tissue cell.
     *
     * @param pMutationState the mutation state of the cell
     * @param pCellCycleModel  the cell cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     * @param archiving  whether this constructor is being called by the archiver - do things slightly differently! (defaults to false)
     * @param cellPropertyCollection the cell property collection (defaults to NULL)
     */
    TissueCell(boost::shared_ptr<AbstractCellMutationState> pMutationState,
               AbstractCellCycleModel* pCellCycleModel,
               bool archiving=false,
               CellPropertyCollection cellPropertyCollection=CellPropertyCollection());

    /**
     * Destructor, which frees the memory allocated for our cell cycle model.
     */
    ~TissueCell();

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
     * Get method for #mpMutationState.
     */
    boost::shared_ptr<AbstractCellMutationState> GetMutationState() const;

    /**
     * Set method for #mpMutationState.
     *
     * @param pMutationState the cell's mutation state
     */
    void SetMutationState(boost::shared_ptr<AbstractCellMutationState> pMutationState);

    /**
     * @return reference to #mCellPropertyCollection.
     */
    CellPropertyCollection& rGetCellPropertyCollection();

    /**
     * @return reference to #mCellPropertyCollection (used in archiving).
     */
    const CellPropertyCollection& rGetCellPropertyCollection() const;

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
    TissueCellPtr Divide();

    /**
     * Make the cell enter apoptosis and sets #mDeathTime using the apoptosis
     * time as defined in the TissueConfig singleton.
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
    double GetTimeUntilDeath() const;

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


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TissueCell)

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
    const boost::shared_ptr<AbstractCellMutationState> p_mutation_state = t->GetMutationState();
    ar << p_mutation_state;

    const AbstractCellCycleModel* const p_cell_cycle_model = t->GetCellCycleModel();
    ar << p_cell_cycle_model;

    const CellPropertyCollection* const p_cell_property_collection = &(t->rGetCellPropertyCollection());
    ar << p_cell_property_collection;
}

/**
 * De-serialize constructor parameters and initialize a TissueCell.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, TissueCell * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    ar >> p_mutation_state;

    AbstractCellCycleModel* p_cell_cycle_model;
    ar >> p_cell_cycle_model;

    bool archiving = true;

    CellPropertyCollection* p_cell_property_collection;
    ar >> p_cell_property_collection;

    // Invoke inplace constructor to initialize instance
    ::new(t)TissueCell(p_mutation_state, p_cell_cycle_model, archiving, *p_cell_property_collection);
}
}
} // namespace ...

#endif /*TISSUECELL_HPP_*/
