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
#ifndef ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_
#define ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>

#include "AbstractOdeBasedCellCycleModelWithStoppingEvent.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"

/**
 * Oxygen-dependent ODE-based cell cycle model. Published by Alarcon et al.
 * (doi:10.1016/j.jtbi.2004.04.016).
 */
class Alarcon2004OxygenBasedCellCycleModel : public AbstractOdeBasedCellCycleModelWithStoppingEvent
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model and ODE system.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModelWithStoppingEvent>(*this);
    }

public:

    /**
     * Default constructor.
     * 
     * @param pOdeSolver An optional pointer to a cell cycle model ODE solver object (allows the use of different ODE solvers)
     */
    Alarcon2004OxygenBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Constructor used in archiving.
     * 
     * @param unused an unused argument
     */
    Alarcon2004OxygenBasedCellCycleModel(double unused)
    {}

    /**
     * Resets the oxygen-based model to the start of the cell cycle
     * (this model does not cycle naturally). Cells are given a new
     * birth time and cell cycle proteins are reset. Note that the
     * oxygen concentration maintains its current value.
     *
     * Should only be called by the Cell Divide() method.
     */
    virtual void ResetForDivision();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This overridden method sets up a new ODE system.
     */
    void Initialise();

    /**
     * Solve the ODEs up to the current time and return whether a stopping event occurred.
     *
     * @param currentTime the current time
     * @return whether a stopping event occurred
     */
    bool SolveOdeToTime(double currentTime);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Alarcon2004OxygenBasedCellCycleModel)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
    double unused = 0.0;
    ::new(t)Alarcon2004OxygenBasedCellCycleModel(unused);
}
}
} // namespace ...

#endif /*ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_*/
