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

#ifndef CELL_HPP_
#define CELL_HPP_

#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include "CellProliferativeTypes.hpp"
#include "AbstractCellMutationState.hpp"
#include "CellLabel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "AbstractCellCycleModel.hpp"
#include "SimulationTime.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellPropertyCollection.hpp"

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).

class Cell;

/** Cells shouldn't be copied - it doesn't make sense.  So all access is via this pointer type. */
typedef boost::shared_ptr<Cell> CellPtr;

/**
 * Cell is the basic container for all the biological information about a cell.
 * It contains the cell-cycle model and all other biological properties such as mutation
 * state, cell type, whether it is undergoing apoptosis or not.
 *
 * This class should not store any spatial information - cells are linked to space by the AbstractCellPopulation subclasses.
 */
class Cell : boost::noncopyable, public boost::enable_shared_from_this<Cell>
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
        archive & mpCellCycleModel;
        archive & mUndergoingApoptosis;
        archive & mDeathTime;
        archive & mStartOfApoptosisTime;
        archive & mApoptosisTime;
        archive & mIsDead;
        archive & mIsLogged;
        archive & mAncestor;
        archive & mCellId;
        archive & mMaxCellId;
    }

protected:

    /** The cell's property collection. */
    CellPropertyCollection mCellPropertyCollection;

    /** The cell's cell-cycle model. */
    AbstractCellCycleModel* mpCellCycleModel;

    /** An index which is inherited by all children of this cell. */
    unsigned mAncestor;

    /** An identifier which is unique to this cell. */
    unsigned mCellId;

    /** maximum cell identifier. */
    static unsigned mMaxCellId;

    /** When the cell will/did die. */
    double mDeathTime;

    /** When the cell was commanded to start apoptosis. */
    double mStartOfApoptosisTime;

    /** The time it takes for a cell to fully undergo apoptosis. Has units of hours. */
    double mApoptosisTime;

    /** Whether the cell is currently in apoptosis - don't divide. */
    bool mUndergoingApoptosis;

    /**
     * Whether the cell is dead or not (they exist in the CellPopulation until they are
     * removed by AbstractCellPopulation::RemoveDeadCells().
     */
    bool mIsDead;

    /** Whether the cell is being tracked specially. */
    bool mIsLogged;

public:

    /**
     * Create a new cell.
     *
     * @param pMutationState the mutation state of the cell
     * @param pCellCycleModel  the cell-cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     * @param archiving  whether this constructor is being called by the archiver - do things slightly differently! (defaults to false)
     * @param cellPropertyCollection the cell property collection (defaults to NULL)
     */
    Cell(boost::shared_ptr<AbstractCellProperty> pMutationState,
               AbstractCellCycleModel* pCellCycleModel,
               bool archiving=false,
               CellPropertyCollection cellPropertyCollection=CellPropertyCollection());

    /**
     * Destructor, which frees the memory allocated for our cell-cycle model.
     */
    ~Cell();

    /**
     * Set the birth time of the cell - can be negative so that your cells have an age when a simulation begins
     *
     * @param birthTime  The time the cell was born (in hours)
     */
    void SetBirthTime(double birthTime);

    /**
     * Change the cell-cycle model used. This takes effect immediately.
     *
     * @param pCellCycleModel pointer to the cell-cycle model to use
     */
    void SetCellCycleModel(AbstractCellCycleModel* pCellCycleModel);

    /**
     * Returns a pointer to the Cell's cell-cycle model.
     */
    AbstractCellCycleModel* GetCellCycleModel() const;

    /**
     * Calls Initialise on the cell-cycle model associated with this cell.
     */
    void InitialiseCellCycleModel();

    /**
     * Get the cell's age from its cell-cycle model.
     */
    double GetAge() const;

    /**
     * Get the cell's birth time from its cell-cycle model.
     */
    double GetBirthTime() const;

    /**
     * Get the time at which apoptosis was commanded to start.
     */
    double GetStartOfApoptosisTime() const;

    /**
     * @return mApoptosisTime
     */
    double GetApoptosisTime() const;

    /**
     * Set mApoptosisTime.
     *
     * @param apoptosisTime the new value of mApoptosisTime
     */
    void SetApoptosisTime(double apoptosisTime);

    /**
     * Get the cell's current mutation state.
     */
    boost::shared_ptr<AbstractCellMutationState> GetMutationState() const;

    /**
     * Set the cell's current mutation state.
     *
     * @param pMutationState the cell's new mutation state
     */
    void SetMutationState(boost::shared_ptr<AbstractCellProperty> pMutationState);

    /**
     * @return reference to #mCellPropertyCollection.
     */
    CellPropertyCollection& rGetCellPropertyCollection();

    /**
     * @return reference to #mCellPropertyCollection (used in archiving).
     */
    const CellPropertyCollection& rGetCellPropertyCollection() const;

    /**
     * Add a cell property to the cell. Use this method instead of calling
     *     rGetCellPropertyCollection().AddProperty()
     * directly, to ensure that the cell property keeps correct track of the
     * number of cells with it (if this is done).
     *
     * @param rProperty the property to add
     */
    void AddCellProperty(const boost::shared_ptr<AbstractCellProperty>& rProperty);

    /**
     * Remove a cell property of the given type. Use this method instead of
     * calling
     *     rGetCellPropertyCollection().AddProperty()
     * directly, to ensure that the cell property keeps correct track of the
     * number of cells with it (if this is done).
     */
    template<typename CLASS>
    void RemoveCellProperty()
    {
        bool cell_has_property = false;

        for (std::set<boost::shared_ptr<AbstractCellProperty> >::iterator property_iter = mCellPropertyCollection.Begin();
             property_iter != mCellPropertyCollection.End();
             ++property_iter)
        {
            if ((*property_iter)->IsType<CLASS>())
            {
                cell_has_property = true;
                (*property_iter)->DecrementCellCount();
                break;
            }
        }

        if (cell_has_property)
        {
            mCellPropertyCollection.RemoveProperty<CLASS>();
        }
    }

    /**
     * Test whether the cell property collection contains a property that has the exact type CLASS.
     * Just calls mCellPropertyCollection.HasProperty().
     */
    template<typename CLASS>
    bool HasCellProperty() const
    {
        return mCellPropertyCollection.HasProperty<CLASS>();
    }

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
    CellPtr Divide();

    /**
     * Make the cell enter apoptosis and sets #mDeathTime using the apoptosis
     * time as defined by mApoptosisTime.
     *
     * @param setDeathTime whether we tell the cell exactly when to die (defaults to true)
     */
    void StartApoptosis(bool setDeathTime=true);

    /**
     * This labels the cell as dead, it does not delete the cell, it remains
     * in the CellPopulation until AbstractCellPopulation::RemoveDeadCells() is called.
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
     * Give the Cell an index which it passes to its children.
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
CHASTE_CLASS_EXPORT(Cell)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Cell.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Cell * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const boost::shared_ptr<AbstractCellMutationState> p_mutation_state = t->GetMutationState();
    ar & p_mutation_state;

    const AbstractCellCycleModel* const p_cell_cycle_model = t->GetCellCycleModel();
    ar & p_cell_cycle_model;

    const CellPropertyCollection& r_cell_property_collection = t->rGetCellPropertyCollection();
    ar & r_cell_property_collection;
}

/**
 * De-serialize constructor parameters and initialize a Cell.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Cell * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    ar & p_mutation_state;

    AbstractCellCycleModel* p_cell_cycle_model;
    ar & p_cell_cycle_model;

    bool archiving = true;

    CellPropertyCollection cell_property_collection;
    ar & cell_property_collection;

    // Invoke inplace constructor to initialize instance
    ::new(t)Cell(p_mutation_state, p_cell_cycle_model, archiving, cell_property_collection);
}
}
} // namespace ...

#endif /*CELL_HPP_*/
