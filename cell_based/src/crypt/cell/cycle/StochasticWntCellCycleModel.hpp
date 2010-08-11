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
#ifndef STOCHASTICWNTCELLCYCLEMODEL_HPP_
#define STOCHASTICWNTCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "WntCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * Wnt-dependent cell cycle model with a stochastic G2 duration.
 *
 * Note that this class uses C++'s default copying semantics, and so
 * doesn't implement a copy constructor or operator=.
 */
class StochasticWntCellCycleModel : public WntCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model and member variables.
     * Used by boost, never directly by chaste code.
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
        archive & boost::serialization::base_object<WntCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;

        archive & mG2Duration;
    }

    /** The duration of the G2 phase, set stochastically. */
    double mG2Duration;

    /**
     * This method overrides that in WntCellCycleModel and
     * introduces the stochastic element of this class.
     *
     * We allow the duration of the G2 phase of the cell cycle to
     * vary as a normal random deviate with a mean of its deterministic
     * duration, a standard deviation of 0.9 hours, and a cutoff to
     * ensure that it is greater than some minimum value.
     *
     * @return the duration of the G2 phases of the cell cycle.
     */
    void SetG2Duration();

public:

    /**
     * The standard constructor called in tests.
     * 
     * @param pOdeSolver An optional pointer to a cell cycle model ODE solver object (allows the use of different ODE solvers)
     */
    StochasticWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * A private constructor for archiving.
     *
     * @param rParentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
     * @param pMutationState the mutation state of the cell (used by ODEs)
     * @param rDimension the spatial dimension
     */
    StochasticWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver,
                                const std::vector<double>& rParentProteinConcentrations,
                                boost::shared_ptr<AbstractCellMutationState> pMutationState,
                                const unsigned& rDimension);

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Set the duration of the G2 phase for the daughter cell.
     */
    void InitialiseDaughterCell();

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This overridden method sets up a new WntCellCycleOdeSystem,
     * sets the cell type according to the current beta catenin level
     * and sets a random G2 duration.
     */
    void Initialise();

    /**
     * Reset cell cycle model by calling AbstractOdeBasedCellCycleModel::ResetForDivision()
     * and setting a new random G2 duration.
     */
    void ResetForDivision();

    /**
     * Get the duration of the G2 phase.
     */
    double GetG2Duration();
};


#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(StochasticWntCellCycleModel)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a StochasticWntCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const StochasticWntCellCycleModel * t, const unsigned int file_version)
{
    const boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_ode_solver = t->GetOdeSolver();
    ar & p_ode_solver;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a StochasticWntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, StochasticWntCellCycleModel * t, const unsigned int file_version)
{
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_ode_solver;
    ar & p_ode_solver;

    /**
     * Invoke inplace constructor to initialise an instance of StochasticWntCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor, provided
     * they are valid parameter values, since the state loaded later from the archive will
     * overwrite their effect in this case.
     */

    std::vector<double> state_vars;
    for (unsigned i=0; i<9; i++)
    {
        state_vars.push_back(0.0);
    }

    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    unsigned dimension = 1;
    ::new(t)StochasticWntCellCycleModel(p_ode_solver, state_vars, p_mutation_state, dimension);
}
}
} // namespace

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(StochasticWntCellCycleModel)

#endif /*STOCHASTICWNTCELLCYCLEMODEL_HPP_*/
