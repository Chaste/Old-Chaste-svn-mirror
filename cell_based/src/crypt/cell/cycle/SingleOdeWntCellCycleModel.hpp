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
#ifndef SINGLEODEWNTCELLCYCLEMODEL_HPP_
#define SINGLEODEWNTCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cfloat>

#include "SimpleWntCellCycleModel.hpp"
#include "Mirams2010WntOdeSystem.hpp"

#include "CellCycleModelOdeSolver.hpp"
#ifdef CHASTE_CVODE
#include "CvodeAdaptor.hpp"
#endif //CHASTE_CVODE
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "HeunIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/**
 * Wnt-dependent cell cycle model. Needs to operate with a WntConcentration
 * singleton object.
 *
 * This model has a constant length M phase, runs ODEs to decide when
 * to finish G1 phase then adds time for S and G2 phases. The CellProliferativeType is
 * updated dependent on the concentration of beta-catenin (given by one
 * of the ODEs).
 */
class SingleOdeWntCellCycleModel : public SimpleWntCellCycleModel
{
private:
    friend class boost::serialization::access;
    /**
     * Save the cell cycle model and ODE system to archive.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        assert(mpOdeSystem);
        archive & boost::serialization::base_object<SimpleWntCellCycleModel>(*this);
        /**
         * Reference can be read or written into once mpOdeSystem has been set up
         * mpOdeSystem isn't set up by the first constructor, but is by the second
         * which is now utilised by the load_construct at the bottom of this file.
         *
         * Note mpOdeSystem itself is not archived just the current values of the
         * state variables...
         */
        boost::shared_ptr<AbstractCellMutationState> p_mutation_state = static_cast<Mirams2010WntOdeSystem*>(mpOdeSystem)->GetMutationState();
        archive & p_mutation_state;
        archive & mBetaCateninDivisionThreshold;
        archive & mLastTime;
        archive & mpOdeSystem->rGetStateVariables();
    }
    /**
     * Load the cell cycle model and ODE system from archive.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        /*
         * The ODE system is set up by the archiving constructor, so we can set the mutation state
         * here. This is a horrible hack, but avoids having to regenerate test archives.
         */
        assert(mpOdeSystem);
        archive & boost::serialization::base_object<SimpleWntCellCycleModel>(*this);
        boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
        archive & p_mutation_state;
        static_cast<Mirams2010WntOdeSystem*>(mpOdeSystem)->SetMutationState(p_mutation_state);

        archive & mBetaCateninDivisionThreshold;
        archive & mLastTime;
        archive & mpOdeSystem->rGetStateVariables();
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /**
     * Pointer to the ODE system developed by Mirams et al (2010).
     */
    Mirams2010WntOdeSystem* mpOdeSystem;

    /**
     * The ODE solver.
     * 
     * Subclasses need to set this in their constructor to point to an instance
     * of a suitable class. See for example the CellCycleModelOdeSolver class.
     */
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> mpOdeSolver;

    /**
     * The cell differentiates when the beta-catenin level drops
     * below this value. It is hard coded in
     * Initialise() because there are so many constructors.
     *
     * Set and Get methods are also provided.
     */
    double mBetaCateninDivisionThreshold;

    /** The last time at which the ODEs were solved up to */
    double mLastTime;

    /**
     * Called by ::Initialise() and ::UpdateCellProliferativeType() only.
     * Updates the mpCell::mCellProliferativeType to match mpOdeSystem's
     * beta-catenin levels
     *
     * This carries out the work for ::UpdateCellProliferativeType();
     * But does not check the current time so it can be used by the initialise method.
     */
    void ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();

    /**
     * This method runs the ODEs and updates the beta-catenin level.
     */
    void UpdateBetaCateninLevel();

public:

    /**
     * Default constructor.
     * 
     * @param pOdeSolver An optional pointer to a cell cycle model ODE solver object (allows the use of different ODE solvers)
     */
    SingleOdeWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * A 'private' constructor for archiving.
     * 
     * @param pOdeSolver a pointer to a cell cycle model ODE solver object (allows the use of different ODE solvers)
     * @param rProteinConcs a std::vector of doubles of the protein concentrations (see VanLeeuwen2009WntSwatCellCycleOdeSystem)
     * @param pMutationState the mutation state of the cell (used by ODEs)
     * @param rDimension the spatial dimension
     * @param useTypeDependentG1 whether to make the duration of G1 phase dependent on the cell's proliferative type (defaults to false)
     */
    SingleOdeWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver,
                               std::vector<double>& rProteinConcs,
                               boost::shared_ptr<AbstractCellMutationState> pMutationState,
                               unsigned& rDimension,
                               bool useTypeDependentG1 = false);

    /**
     * Destructor.
     */
    ~SingleOdeWntCellCycleModel();

    /**
     * Copy constructor.
     *
     * This is important to make a copy of the ODE system instead of
     * giving the copied cell cycle model a pointer to the same ODE.
     *
     * @param rOtherModel  the one to copy
     */
    SingleOdeWntCellCycleModel(const SingleOdeWntCellCycleModel& rOtherModel);

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This overridden method sets up a new WntCellCycleOdeSystem,
     * sets the cell type according to the current beta catenin level
     * and sets a random G2 duration.
     */
    void Initialise();

    /**
     * This specialisation updates the beta-catenin level
     */
     void UpdateCellCyclePhase();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Return the total beta-catenin concentration
     */
    double GetBetaCateninConcentration();

    /**
     * Set #mBetaCateninDivisionThreshold.
     *
     * @param betaCateninDivisionThreshold to be set
     */
    void SetBetaCateninDivisionThreshold(double betaCateninDivisionThreshold);

    /**
     * Get #mBetaCateninDivisionThreshold.
     */
    double GetBetaCateninDivisionThreshold();

     /**
      * @return mpOdeSolver (used in archiving).
      */
    const boost::shared_ptr<AbstractCellCycleModelOdeSolver> GetOdeSolver() const;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SingleOdeWntCellCycleModel)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SingleOdeWntCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const SingleOdeWntCellCycleModel * t, const unsigned int file_version)
{
    const boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_ode_solver = t->GetOdeSolver();
    ar & p_ode_solver;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SingleOdeWntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SingleOdeWntCellCycleModel * t, const unsigned int file_version)
{
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_ode_solver;
    ar & p_ode_solver;

    /**
     * Invoke inplace constructor to initialise an instance of SingleOdeWntCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later from the
     * archive will overwrite their effect in this case.
     */

    std::vector<double> state_vars;
    for (unsigned i=0; i<3; i++)
    {
        state_vars.push_back(0.0);
    }

    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    unsigned dimension = 1;
    ::new(t)SingleOdeWntCellCycleModel(p_ode_solver, state_vars, p_mutation_state, dimension);
}
}
} // namespace

#ifdef CHASTE_CVODE
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, SingleOdeWntCellCycleModel, CvodeAdaptor)
#endif //CHASTE_CVODE
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, SingleOdeWntCellCycleModel, BackwardEulerIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, SingleOdeWntCellCycleModel, EulerIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, SingleOdeWntCellCycleModel, HeunIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, SingleOdeWntCellCycleModel, RungeKutta2IvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver)

#endif /*SINGLEODEWNTCELLCYCLEMODEL_HPP_*/
