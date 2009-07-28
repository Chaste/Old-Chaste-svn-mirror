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
#ifndef ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_
#define ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractCellCycleModel.hpp"
#include "AbstractOdeSystem.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * This class contains all the things common to standard cell cycle
 * ODE models for intracellular protein concentrations (along the lines
 * of Tyson & Novak), such as solving the ODEs until a stopping condition
 * is met.
 */
class AbstractOdeBasedCellCycleModel : public AbstractCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        assert(mpOdeSystem!=NULL);
        archive & mpOdeSystem->rGetStateVariables();
        archive & mLastTime;
        archive & mDivideTime;
        archive & mFinishedRunningOdes;
        archive & mG2PhaseStartTime;
    }

protected:

    /** The system of ODEs for the cell cycle model. */
    AbstractOdeSystem *mpOdeSystem;

    /** The last time the cell cycle ODEs were evaluated.*/
    double mLastTime;

    /** The time at which the cell should divide - Set this to DBL_MAX in constructor.*/
    double mDivideTime;

    /** Whether the cell cycle model is currently in a delay (not solving ODEs).*/
    bool mFinishedRunningOdes;

    /** The start time for the G2 phase */
    double mG2PhaseStartTime;

public:

    /**
     * Creates an AbstractOdeBasedCellCycleModel, calls SetBirthTime on the
     * AbstractCellCycleModel to make sure that can be set 'back in time' for
     * cells which did not divide at the current time.
     *
     * @param lastTime  The birth time of the cell / last time model was evaluated (defaults to the current SimulationTime)
     */
    AbstractOdeBasedCellCycleModel(double lastTime = SimulationTime::Instance()->GetTime());

    /**
     * Copy constructor.
     *
     * This is needed because we store and manage a pointer to an ODE system.
     * Note that this class doesn't actually copy the ODE system, because each
     * subclass will use a different type.  Hence subclasses *must* copy their
     * own ODE system in their copy constructor.
     *
     * @param rOtherModel the cell cycle model being copied.
     */
    AbstractOdeBasedCellCycleModel(const AbstractOdeBasedCellCycleModel& rOtherModel);

    /**
     * This destructor deletes the mpOdeSystem.
     */
    virtual ~AbstractOdeBasedCellCycleModel();

    /**
     * Default UpdateCellCyclePhase() method for an ODE-based cell cycle model.
     * This method calls SolveOdeToTime() for G1 phase and adds time for the other phases.
     *
     * Can be overridden if they should do something more subtle.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * This method must be implemented by each subclass - solves the ODEs to a given time.
     *
     * @param currentTime the current time
     *
     * @return Whether a stopping event occurred.
     */
    virtual bool SolveOdeToTime(double currentTime)=0;

    /**
     * This method must be implemented by each subclass
     *
     * When the ODEs have reached a stopping event it returns the time at which
     * the ODEs stopped running so a delay can be added in for S-G2-M phases if necessary.
     *
     * @return The time at which the ODE reached its stopping event.
     */
    virtual double GetOdeStopTime()=0;

    /**
     * This overrides the AbstractCellCycleModel::SetBirthTime(double birthTime)
     * because an ODE based cell cycle model has more to reset...
     *
     * @param birthTime the simulation time when the cell was born
     */
    void SetBirthTime(double birthTime);

    /**
     * Returns the protein concentrations at the current time (useful for tests)
     *
     * NB: Will copy the vector - you can't use this to modify the concentrations.
     */
    std::vector<double> GetProteinConcentrations() const;

    /**
     * Sets the protein concentrations and time when the model was last evaluated - should only be called by tests
     *
     * @param lastTime the SimulationTime at which the protein concentrations apply
     * @param proteinConcentrations a standard vector of doubles of protein concentrations
     *
     */
    void SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations);

    /**
     * For a naturally cycling model this does not need to be overridden in the
     * subclasses. But most models should override this function and then
     * call AbstractOdeBasedCellCycleModel::ResetForDivision() from inside their version.
     */
    virtual void ResetForDivision();
};

BOOST_IS_ABSTRACT(AbstractOdeBasedCellCycleModel)

#endif /*ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_*/
