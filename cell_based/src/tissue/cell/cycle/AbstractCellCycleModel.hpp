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
#ifndef ABSTRACTCELLCYCLEMODEL_HPP_
#define ABSTRACTCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>

#include "CellProliferativeTypes.hpp"
#include "CellCyclePhases.hpp"
#include "SimulationTime.hpp"
#include "TissueConfig.hpp"
#include "TissueCell.hpp"

class TissueCell; // Circular definition (cells need to know about cycle models and vice-versa)
typedef boost::shared_ptr<TissueCell> TissueCellPtr;

/**
 * The AbstractCellCycleModel contains basic information to all cell cycle models.
 * It handles assignment of birth time, cell cycle phase and a TissueCell.
 * 
 * Cell cycle models are noncopyable since cells are noncopyable.
 */
class AbstractCellCycleModel : boost::noncopyable
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Make sure the SimulationTime and TissueConfig singletons get saved too
        SimulationTime* p_time = SimulationTime::Instance();
        archive & *p_time;
        archive & p_time;
        TissueConfig* p_params = TissueConfig::Instance();
        archive & *p_params;
        archive & p_params;

        // DO NOT archive & mpCell; -- The CellCycleModel is only ever archived from the Cell
        // which knows this and it is handled in the load_construct of TissueCell.
        archive & mBirthTime;
        archive & mCurrentCellCyclePhase;
        archive & mG1Duration;
        archive & mReadyToDivide;
        archive & mDimension;
        archive & mCellProliferativeType;
        archive & mMinimumGapDuration;
    }

protected:

    /** The cell that this model is associated with. */
    TissueCellPtr mpCell;

    /**
     * The time that the cell began to split from its parent
     * (i.e. @b beginning of M phase NOT the end)
     */
    double mBirthTime;

    /** The phase of the cell cycle that this model is in (specified in CellCyclePhases.hpp) */
    CellCyclePhase mCurrentCellCyclePhase;

    /**
     * How long the G1 phase lasts for.
     * Not necessarily a fixed value.
     */
    double mG1Duration;

    /**
     * Whether the cell is currently ready to undergo division.
     */
    bool mReadyToDivide;

    /**
     *  Spatial dimension being used in simulation (defaults to 0, set with SetDimension).
     */
    unsigned mDimension;

    /**
     * The cell type - defined in CellProliferativeTypes.hpp.
     */
    CellProliferativeType mCellProliferativeType;

    /**
     * Minimum possbile duration of either of the gap phases (G1 or G2).
     * Has units of hours.
     * 
     * Used to guarantee a strictly positive duration in cell cycle models that
     * use normal random deviates for G1 or G2 phases.
     */
    double mMinimumGapDuration;

public:

    /**
     * Sets up a new AbstractCellCycleModel, gives it a birth time of the
     * current simulation time (which is overwritten by some subclasses)
     */
    AbstractCellCycleModel();

    /**
     * Base class with virtual methods needs a virtual destructor. The destructor
     * does not delete mpCell. Instead, the cell takes responsibility for deleting
     * the cell cycle model when it is destroyed.
     */
    virtual ~AbstractCellCycleModel();

    /**
     * Gives the cell cycle model a pointer to its host cell.
     *
     * Some cell cycle models pass this pointer to other classes (e.g. WntConcentration),
     * which use this information to determine other information based upon the location
     * of the cell (e.g. the Wnt concentration at this location).
     *
     * @param pCell pointer to the cell
     */
    void SetCell(TissueCellPtr pCell);

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This method will be called precisely once per cell set up in the initial
     * tissue. It is not called on cell division; use ResetForDivision(),
     * CreateCellCycleModel() and InitialiseDaughterCell() for that.
     *
     * By the time this is called, a Tissue will have been set up, so the model
     * can know where its cell is located in space. If relevant to the simulation,
     * the CellwiseData and WntConcentration singletons will also have been initialised.
     */
    virtual void Initialise();

    /**
     * Initialise the new daughter cell's cycle model after a cell division.
     *
     * This is called by TissueCell::Divide once the new cell object
     * has been fully created, to perform any initialisation of the
     * cell cycle which requires access to the cell.
     *
     * Note that much initialisation can be performed using the
     * combination of ResetForDivision() (called on the parent prior to
     * division) and CreateCellCycleModel() (called on the reset
     * parent to create the new cell cycle model object).
     */
    virtual void InitialiseDaughterCell();

    /**
     * @return The cell which plays host to this cell cycle model.
     */
    TissueCellPtr GetCell();

    /**
     * Set the cell's time of birth (usually not required as it should be inside
     * the indivdual cell-cycle-model-constructor, but useful for tests).
     *
     * @param birthTime the simulation time at this cell's birth.
     *
     * (This function is overridden in AbstractOdeBasedCellCycleModel).
     */
    virtual void SetBirthTime(double birthTime);

    /**
     * Set the spatial dimension.
     *
     * @param dimension
     */
    void SetDimension(unsigned dimension);

    /**
     *  Get the dimension this cell cycle model thinks the simulation is in.
     */
    unsigned GetDimension();

    /**
     * @return the time at which the cell was born.
     */
    double GetBirthTime() const;

    /**
     * Returns the cell's age.
     */
    double GetAge();

    /**
     * Determine whether the cell is ready to divide (enter M phase).
     *
     * The intention is that this method is called precisely once at
     * each timestep of the simulation. However this does not appear
     * to always be the case at present, and so it can cope with more
     * unusual usage patterns.
     */
    bool ReadyToDivide();

    /**
     * This method must be implemented by subclasses in order to set the phase
     * the cell cycle model is currently in. It is called from ReadyToDivide()
     * just prior to deciding whether to divide the cell based on how far through
     * the cell cycle it is, i.e. whether it has completed M, G1, S and G2 phases.
     */
    virtual void UpdateCellCyclePhase()=0;

    /**
     * Each cell cycle model must be able to be reset 'after' a cell division.
     *
     * Actually, this method is called from TissueCell::Divide() to
     * reset the cell cycle just before the daughter cell is created.
     * CreateCellCycleModel() can then clone our state to generate a
     * cell cycle model instance for the daughter cell.
     */
    virtual void ResetForDivision();

    /**
     * Builder method to create new instances of the cell cycle model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     *
     * This method is called by TissueCell::Divide() to create a cell
     * cycle model for the daughter cell.  Note that the parent cell
     * cycle model will have had ResetForDivision() called just before
     * CreateCellCycleModel() is called, so performing an exact copy of the
     * parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     */
    virtual AbstractCellCycleModel* CreateCellCycleModel()=0;

    /**
     * @return the current cell cycle phase
     */
    CellCyclePhase GetCurrentCellCyclePhase();

    /**
     * @return the duration of the G1 phase of the cell cycle
     */
    virtual double GetG1Duration();

    /**
     * @return the duration of the S phase of the cell cycle
     */
    virtual double GetSDuration();

    /**
     * @return the duration of the G2 phase of the cell cycle
     */
    virtual double GetG2Duration();

    /**
     * @return the duration of the M phase of the cell cycle
     */
    virtual double GetMDuration();

    /**
     * Return the typical cell cycle duration for a transit cell, in hours.
     * This method is overridden in some subclasses.
     */
    virtual double GetAverageTransitCellCycleTime();

    /**
     * Return the typical cell cycle duration for a stem cell, in hours.
     * This method is overridden in some subclasses.
     */
    virtual double GetAverageStemCellCycleTime();

    /**
     * Whether a cell with this cell cycle model is able to fully (terminally) differentiate.
     */
    virtual bool CanCellTerminallyDifferentiate();

    /**
     * Get method for #mCellProliferativeType.
     */
    CellProliferativeType GetCellProliferativeType() const;

    /**
     * Set method for #mCellProliferativeType.
     *
     * @param cellType the cell's type
     */
    void SetCellProliferativeType(CellProliferativeType cellType);

    /**
     * @return mMinimumGapDuration
     */
    double GetMinimumGapDuration();

    /**
     * Set mMinimumGapDuration.
     * 
     * @param minimumGapDuration the new value of mMinimumGapDuration
     */
    void SetMinimumGapDuration(double minimumGapDuration);
};

CLASS_IS_ABSTRACT(AbstractCellCycleModel)

#endif /*ABSTRACTCELLCYCLEMODEL_HPP_*/
