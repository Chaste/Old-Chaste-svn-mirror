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
#ifndef TYSONNOVAKCELLCYCLEMODEL_HPP_
#define TYSONNOVAKCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"

#include <iostream>

#include "AbstractOdeBasedCellCycleModelWithStoppingEvent.hpp"
#include "TysonNovak2001OdeSystem.hpp"
#include "Exception.hpp"

#include "CellCycleModelOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "HeunIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/**
 *  Tyson-Novak 2001 cell cycle model, taken from the version at  doi:10.1006/jtbi.2001.2293
 *
 *  Note that this is not a model for murine or human colonic-cell cycling, but is
 *  included in chaste as one of the most commonly known ODE based cell cycle models.
 *
 *  Time taken to progress through the cycle is deterministic and given by
 *  an ODE system independent of external factors.
 */
class TysonNovakCellCycleModel : public AbstractOdeBasedCellCycleModelWithStoppingEvent
{
private:

    friend class TestOdeBasedCellCycleModels;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
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
    TysonNovakCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Copy constructor.
     *
     * Also creates a copy of our ODE system.
     *
     * @param rOtherModel the instance being copied.
     */
    TysonNovakCellCycleModel(const TysonNovakCellCycleModel& rOtherModel);

    /**
     * Reset cell cycle model by calling AbstractOdeBasedCellCycleModelWithStoppingEvent::ResetForDivision()
     * and setting initial conditions for protein concentrations.
     */
    void ResetForDivision();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Solve the ODEs up to the current time and return whether a stopping event occurred.
     *
     * @param currentTime the current time
     * @return whether a stopping event occurred
     */
    bool SolveOdeToTime(double currentTime);

    /**
     * Get the duration of the cell's S phase.
     */
    double GetSDuration();

    /**
     * Get the duration of the cell's G2 phase.
     */
    double GetG2Duration();

    /**
     * Get the duration of the cell's M phase.
     */
    double GetMDuration();

    /**
     * If the daughter cell type is stem, change it to transit.
     */
    void InitialiseDaughterCell();

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     */
    double GetAverageStemCellCycleTime();

    /**
     * Overridden CanCellTerminallyDifferentiate() method.
     */
    virtual bool CanCellTerminallyDifferentiate();
};


#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TysonNovakCellCycleModel)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a TysonNovakCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const TysonNovakCellCycleModel * t, const unsigned int file_version)
{
    const boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_ode_solver = t->GetOdeSolver();
    ar & p_ode_solver;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a TysonNovakCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, TysonNovakCellCycleModel * t, const unsigned int file_version)
{
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_ode_solver;
    ar & p_ode_solver;

    ::new(t)TysonNovakCellCycleModel(p_ode_solver);
}
}
} // namespace ...


#ifdef CHASTE_CVODE
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, TysonNovakCellCycleModel, CvodeAdaptor)
#endif //CHASTE_CVODE
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, TysonNovakCellCycleModel, EulerIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, TysonNovakCellCycleModel, HeunIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, TysonNovakCellCycleModel, RungeKutta2IvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, TysonNovakCellCycleModel, RungeKutta4IvpOdeSolver)

#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
