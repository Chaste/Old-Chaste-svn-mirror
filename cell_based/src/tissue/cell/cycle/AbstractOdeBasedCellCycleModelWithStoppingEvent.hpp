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
#ifndef ABSTRACTODEBASEDCELLCYCLEMODELWTIHSTOPPINGEVENT_HPP_
#define ABSTRACTODEBASEDCELLCYCLEMODELWTIHSTOPPINGEVENT_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "AbstractOdeSystem.hpp"

/**
 * This class contains all the things common to standard cell cycle
 * ODE models for intracellular protein concentrations (along the lines
 * of Tyson & Novak), such as solving the ODEs until a stopping condition
 * is met.
 */
class AbstractOdeBasedCellCycleModelWithStoppingEvent : public AbstractOdeBasedCellCycleModel
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
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
//        archive & mDivideTime;
//        archive & mFinishedRunningOdes;
//        archive & mG2PhaseStartTime;
    }

protected:

//    /** The time at which the cell should divide - Set this to DBL_MAX in constructor.*/
//    double mDivideTime;
//
//    /** Whether the cell cycle model is currently in a delay (not solving ODEs).*/
//    bool mFinishedRunningOdes;
//
//    /** The start time for the G2 phase */
//    double mG2PhaseStartTime;

public:

    /**
     * Creates an AbstractOdeBasedCellCycleModel, calls SetBirthTime on the
     * AbstractCellCycleModel to make sure that can be set 'back in time' for
     * cells which did not divide at the current time.
     *
     * @param lastTime  The birth time of the cell / last time model was evaluated (defaults to the current SimulationTime)
     * @param pOdeSolver An optional pointer to a cell cycle model ODE solver object (allows the use of different ODE solvers)
     */
    AbstractOdeBasedCellCycleModelWithStoppingEvent(double lastTime = SimulationTime::Instance()->GetTime(),
                                                    boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

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
    AbstractOdeBasedCellCycleModelWithStoppingEvent(const AbstractOdeBasedCellCycleModelWithStoppingEvent& rOtherModel);

    /**
     * This destructor deletes the mpOdeSystem.
     */
    virtual ~AbstractOdeBasedCellCycleModelWithStoppingEvent();

};

CLASS_IS_ABSTRACT(AbstractOdeBasedCellCycleModelWithStoppingEvent)

#endif /*ABSTRACTODEBASEDCELLCYCLEMODELWTIHSTOPPINGEVENT_HPP_*/
