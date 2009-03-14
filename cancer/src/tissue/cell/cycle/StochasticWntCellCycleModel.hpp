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
#ifndef STOCHASTICWNTCELLCYCLEMODEL_HPP_
#define STOCHASTICWNTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "WntCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

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
     * @param archive
     * @param version
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

    /**
     * The standard constructor called in tests.
     */
    StochasticWntCellCycleModel();

     /**
     * A private constructor for daughter cells called by the CreateDaughterCellCycleModel function
     * (which can be called by TissueCell::CommonCopy() and isn't necessarily being born.
     *
     * @param pParentOdeSystem  to copy the state of
     * @param mutationState the mutation state of the cell (used by ODEs)
     * @param birthTime the simulation time when the cell divided (birth time of parent cell)
     * @param lastTime last time the cell cycle model was evaluated
     * @param inSG2MPhase whether the cell is in S-G2-M (not evaluating ODEs and just waiting)
     * @param readyToDivide whether the cell is ready to divide
     * @param divideTime if in the future this is the time at which the cell is going to divide
     * @param g2Duration the duration of the cell's G2 phase
     */
    StochasticWntCellCycleModel(AbstractOdeSystem* pParentOdeSystem,
                                CellMutationState mutationState,
                                double birthTime,
                                double lastTime,
                                bool inSG2MPhase,
                                bool readyToDivide,
                                double divideTime,
                                double g2Duration);

    /**
     * A private constructor for archiving.
     *
     * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
     * @param mutationState the mutation state of the cell (used by ODEs)
     */
    StochasticWntCellCycleModel(std::vector<double> parentProteinConcentrations,
                                CellMutationState mutationState);

    /**
     * Returns a new StochasticWntCellCycleModel, created with the correct
     * initial conditions.
     *
     * This method should be called just after the parent cell cycle model
     * has been reset.
     *
     * @return pointer to the daughter cell cycle model
     */
    AbstractCellCycleModel* CreateDaughterCellCycleModel();

};


// Declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticWntCellCycleModel)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a WntCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const StochasticWntCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a StochasticWntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, StochasticWntCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of StochasticWntCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    std::vector<double> state_vars;
    for (unsigned i=0; i<9; i++)
    {
        state_vars.push_back(0.0);
    }

    CellMutationState mutation_state = HEALTHY;

    ::new(t)StochasticWntCellCycleModel(state_vars, mutation_state);
}
}
} // namespace

#endif /*STOCHASTICWNTCELLCYCLEMODEL_HPP_*/
