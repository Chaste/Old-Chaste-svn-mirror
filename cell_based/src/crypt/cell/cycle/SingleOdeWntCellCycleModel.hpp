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

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <cfloat>

#include "SimpleWntCellCycleModel.hpp"
#include "Mirams2010WntOdeSystem.hpp"
#include "CryptCellMutationStates.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "WntConcentration.hpp"

/**
 * Wnt-dependent cell cycle model. Needs to operate with a WntConcentration
 * singleton object.
 *
 * This model has a constant length M phase, runs ODEs to decide when
 * to finish G1 phase then adds time for S and G2 phases. The CellProliferativeType is
 * updated dependent on the concentration of beta-catenin (given by one
 * of the ODEs).
 *
 * Note that this class uses C++'s default copying semantics, and so
 * doesn't implement a copy constructor or operator=.
 */
class SingleOdeWntCellCycleModel : public SimpleWntCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        assert(mpOdeSystem!=NULL);
        archive & boost::serialization::base_object<SimpleWntCellCycleModel>(*this);
        /**
         * Reference can be read or written into once mpOdeSystem has been set up
         * mpOdeSystem isn't set up by the first constructor, but is by the second
         * which is now utilised by the load_construct at the bottom of this file.
         *
         * Note mpOdeSystem itself is not archived just the current values of the
         * state variables...
         */
        archive & static_cast<Mirams2010WntOdeSystem*>(mpOdeSystem)->rGetMutationState();
        archive & mpOdeSystem->rGetStateVariables();
        archive & mBetaCateninDivisionThreshold;
        archive & mLastTime;
    }

    Mirams2010WntOdeSystem* mpOdeSystem;

    /**
     * The cell differentiates when the beta-catenin level drops
     * below this value. It is hard coded in
     * Initialise() because there are so many constructors.
     *
     * Set and Get methods are also provided.
     */
    double mBetaCateninDivisionThreshold;

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

#ifdef CHASTE_CVODE
    /** The ODE solver. */
    static CvodeAdaptor msSolver;
#else
    /** The ODE solver. */
    static RungeKutta4IvpOdeSolver msSolver;
#endif //CHASTE_CVODE

    /** The last time at which the ODEs were solved up to */
    double mLastTime;



public:

    /**
     * Default constructor.
     * @param dimension the spatial dimension
     */
    SingleOdeWntCellCycleModel(unsigned dimension, bool useTypeDependentG1Durations=false)
        : SimpleWntCellCycleModel(dimension, useTypeDependentG1Durations),
          mpOdeSystem(NULL),
          mLastTime(DBL_MAX) // Ensure this is set properly before we try to use it.
    {
#ifdef CHASTE_CVODE
        msSolver.SetMaxSteps(10000);
#endif // CHASTE_CVODE
    };


    /**
	 * Normal constructor
	 */
	SingleOdeWntCellCycleModel(std::vector<double>& rProteinConcs, CryptCellMutationState& rMutationState, unsigned& rDimension, bool useTypeDependentG1 = false);

	/**
	 * Destructor
	 */
	~SingleOdeWntCellCycleModel()
	{
	    delete mpOdeSystem;
	}

	/**
	 * Copy constructor...
	 *
	 * This is important to make a copy of the ODE system instead of
	 * giving the copied cell cycle model a pointer to the same ODE...
	 */
	SingleOdeWntCellCycleModel(const SingleOdeWntCellCycleModel& rOtherModel);

	void Initialise();

	void UpdateCellCyclePhase();

	AbstractCellCycleModel* CreateCellCycleModel();


	double GetBetaCateninConcentration(void)
	{
	    return mpOdeSystem->rGetStateVariables()[0] + mpOdeSystem->rGetStateVariables()[1];
	}

	void SetBetaCateninDivisionThreshold(double value)
	{
	    mBetaCateninDivisionThreshold = value;
	}

	double GetBetaCateninDivisionThreshold(void)
	{
	    return mBetaCateninDivisionThreshold;
	}

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
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SingleOdeWntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SingleOdeWntCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of SingleOdeWntCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    std::vector<double> state_vars;
    for (unsigned i=0; i<1; i++)
    {
        state_vars.push_back(0.0);
    }

    CryptCellMutationState mutation_state = HEALTHY;
    unsigned dimension = UINT_MAX;
    ::new(t)SingleOdeWntCellCycleModel(state_vars, mutation_state, dimension);
}
}
} // namespace

#endif /*SINGLEODEWNTCELLCYCLEMODEL_HPP_*/
