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
#ifndef ABSTRACTCELLCYCLEMODEL_HPP_
#define ABSTRACTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include <vector>

#include "CellTypes.hpp"
#include "CellCyclePhases.hpp"
#include "SimulationTime.hpp"
#include "CancerParameters.hpp"
#include "TissueCell.hpp"


// Needs to be included last
#include <boost/serialization/export.hpp>

class TissueCell; // Circular definition (cells need to know about cycle models and vice-versa).

/**
 * The AbstractCellCycleModel contains basic information to all cell cycle models.
 * It handles assignment of birth time, cell cycle phase and a TissueCell.
 */
class AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Make sure the simulation and cancer parameters get saved too
        SimulationTime* p_time = SimulationTime::Instance();
        archive & *p_time;
        archive & p_time;
        CancerParameters* p_params = CancerParameters::Instance();
        archive & *p_params;
        archive & p_params;

        // DO NOT archive & mpCell; -- The CellCycleModel is only ever archived from the Cell
        // which knows this and it is handled in the load_construct of TissueCell.
        archive & mBirthTime;
        archive & mCurrentCellCyclePhase;
        archive & mGeneration;
        archive & mG1Duration;
        archive & mReadyToDivide;
    }


protected:
    /** The cell that this model is associated with */
    TissueCell* mpCell;

    /**
     * The time that the cell began to split from its parent
     * (i.e. beginning of M phase NOT the end)
     */
    double mBirthTime;

    /** The phase of the cell cycle that this model is in (specified in CellCyclePhases.hpp) */
    CellCyclePhase mCurrentCellCyclePhase;

    /** The generation of this cell (STEM cells have a generation of 0) */
    unsigned mGeneration;

    /**
     * How long the G1 phase lasts for.
     * Not necessarily a fixed value...
     */
    double mG1Duration;

    /**
     * Whether the cell is currently ready to undergo division.
     */
    bool mReadyToDivide;

public:

    /**
     * Sets up a new AbstractCellCycleModel, gives it a birth time of the
     * current simulation time (which is overwritten by some subclasses)
     */
    AbstractCellCycleModel()
        : mpCell(NULL),
          mBirthTime(SimulationTime::Instance()->GetDimensionalisedTime()),
          mCurrentCellCyclePhase(M_PHASE),
          mGeneration(0),
          mG1Duration(DBL_MAX),
          mReadyToDivide(false)
    {}

    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractCellCycleModel();

    void SetCell(TissueCell* pCell);

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This method will be called precisely once per cell set up in the initial
     * tissue. It is not called on cell division; use ResetModel,
     * CreateDaughterCellCycleModel and InitialiseDaughterCell for that.
     *
     * By the time this is called, a Tissue will have been set up, so the model
     * can know where its cell is located in space.  If relevant to the simulation,
     * the CellwiseData and WntConcentration singletons will also have been initialised.
     */
    virtual void Initialise()
    {}

    /**
     * Initialise the new daughter cell's cycle model after a cell division.
     *
     * This is called by TissueCell::Divide once the new cell object
     * has been fully created, to perform any initialisation of the
     * cell cycle which requires access to the cell.
     *
     * Note that much initialisation can be performed using the
     * combination of ResetModel (called on the parent prior to
     * division) and CreateDaughterCellCycleModel (called on the reset
     * parent to create the new cell cycle model object).
     */
    virtual void InitialiseDaughterCell()
    {}

    TissueCell* GetCell();

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
     * @return the time at which the cell was born.
     */
    double GetBirthTime() const;

    /**
     * Returns the cell's age...
     */
    double GetAge();

    /**
     * Sets the cell's generation...
     */
    void SetGeneration(unsigned generation);

    /**
     * Returns the cell's generation...
     */
    unsigned GetGeneration() const;

    /**
     * Determine whether the cell is ready to divide (enter M phase).
     *
     * The intention is that this method is called precisely once at
     * each timestep of the simulation.  However this does not appear
     * to always be the case at present, and so it can cope with more
     * unusual usage patterns.
     */
    bool ReadyToDivide();

    /**
     * This method must be implemented by subclasses in order to set the phase
     * the cell cycle model is currently in.  It is called from ReadyToDivide()
     * just prior to deciding whether to divide the cell based on how far through
     * the cell cycle it is, i.e. whether it has completed M, G1, S and G2 phases.
     */
    virtual void UpdateCellCyclePhase()=0;

    /**
     * Each cell cycle model must be able to be reset 'after' a cell division.
     *
     * Actually, this method is called from TissueCell::Divide to
     * reset the cell cycle just before the daughter cell is created.
     * CreateDaughterCellCycleModel can then clone our state to generate a
     * cell cycle model instance for the daughter cell.
     */
    virtual void ResetForDivision();

    /**
     * Builder method to create new instances of the cell cycle model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     *
     * This method is called by the copy constructor and operator= of
     * TissueCell to create a copy of the cell cycle model when
     * copying a cell.  It thus just needs to create any instance of
     * the right class, as operator= on the cell cycle model is then
     * called to ensure the model is copied properly.
     *
     * A default implementation is given here which uses
     * CreateDaughterCellCycleModel, in order to reduce coding effort
     * for the refactor.
     */
    virtual AbstractCellCycleModel *CreateCellCycleModel()
    {
        return CreateDaughterCellCycleModel();
    }

    /**
     * Builder method to create new instances of the cell cycle model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     *
     * This method is called by TissueCell.Divide to create a cell
     * cycle model for the daughter cell.  It thus must thus produce a
     * cell cycle model in a suitable state for a newly-born cell
     * spawned from the 'current' cell.  Note that the parent cell
     * cycle model will have had ResetModel called just before
     * CreateDaughterCellCycleModel is called.
     */
    virtual AbstractCellCycleModel *CreateDaughterCellCycleModel()=0;

    /**
     * @return whether the cell cycle model uses beta-catenin levels in cell cycle model, ie Inge models.
     */
    virtual bool UsesBetaCat();

    /**
     * Returns the protein concentrations at the current time. However in most Cell Cycle models this does not exist.
     * We have a "work-around" such that we throw an error if we try and access it for any other cell type.
     * \todo may be better to use dynamic_cast and/or MI (see #697)
     *
     * NB: Will copy the vector - you can't use this to modify the concentrations.
     */
     virtual std::vector<double> GetProteinConcentrations() const;

    /*
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

};


BOOST_IS_ABSTRACT(AbstractCellCycleModel)


#endif /*ABSTRACTCELLCYCLEMODEL_HPP_*/
