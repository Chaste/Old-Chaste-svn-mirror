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
#ifndef INGEWNTSWATCELLCYCLEMODEL_HPP_
#define INGEWNTSWATCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cfloat>

#include "AbstractOdeSystem.hpp"
#include "AbstractWntOdeBasedCellCycleModel.hpp"
#include "IngeWntSwatCellCycleOdeSystem.hpp"
#include "AbstractCellMutationState.hpp"
#include "Exception.hpp"


/**
 *  Wnt-dependent cell cycle model.
 *
 * Note that this class uses C++'s default copying semantics, and so doesn't implement a copy constructor
 * or operator=.
 */
class IngeWntSwatCellCycleModel : public AbstractWntOdeBasedCellCycleModel
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
    void save(Archive & archive, const unsigned int version) const
    {
        assert(mpOdeSystem);
        archive & boost::serialization::base_object<AbstractWntOdeBasedCellCycleModel>(*this);
        boost::shared_ptr<AbstractCellMutationState> p_mutation_state = static_cast<IngeWntSwatCellCycleOdeSystem*>(mpOdeSystem)->GetMutationState();
        archive & p_mutation_state;
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
    	// The ODE system is set up by the archiving constructor, so we can set the mutation state
    	// here.  This is a horrible hack, but avoids having to regenerate test archives...
    	assert(mpOdeSystem);
        archive & boost::serialization::base_object<AbstractWntOdeBasedCellCycleModel>(*this);
        boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
        archive & p_mutation_state;
        static_cast<IngeWntSwatCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(p_mutation_state);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

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
     * Hypothesis number (1 or 2), concerning the nature of the
     * interactions modelled by the cell cycle ODE system.
     */
    unsigned mHypothesis;

public:

    /**
     * Default constructor.
     */
    IngeWntSwatCellCycleModel();

    /**
     * Copy constructor.
     *
     * Creates an appropriate copy of our ODE system too.
     *
     * @param rOtherModel the instance being copied.
     */
    IngeWntSwatCellCycleModel(const IngeWntSwatCellCycleModel& rOtherModel);

    /**
     * A 'private' constructor for archiving.
     *
     * @param rHypothesis which model hypothesis to use (1 or 2)
     * @param rParentProteinConcentrations a std::vector of doubles of the protein concentrations (see IngeWntSwatCellCycleOdeSystem)
     * @param pMutationState the mutation state of the cell (used by ODEs)
     * @param rDimension the spatial dimension
     */
    IngeWntSwatCellCycleModel(const unsigned& rHypothesis,
                              const std::vector<double>& rParentProteinConcentrations,
                              boost::shared_ptr<AbstractCellMutationState> pMutationState,
                              const unsigned& rDimension);

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * See AbstractCellCycleModel::Initialise()
     *
     * In this case we set up a new ODE system for a daughter cell.
     */
    void Initialise();

    /**
     * Solve the ODE to the current time.
     *
     * @param currentTime the current time
     * @return Whether a stopping event occurred.
     */
    bool SolveOdeToTime(double currentTime);

    /**
     * @return the level of membrane bound beta-catenin. To be used in cell-cell adhesion calculations.
     */
    double GetMembraneBoundBetaCateninLevel();

    /**
     * @return the level of cytoplasmic beta-catenin (including ubiquitinated - awaiting degradation)
     */
    double GetCytoplasmicBetaCateninLevel();

    /**
     * @return the level of nuclear beta-catenin. To be used in transcription
     */
    double GetNuclearBetaCateninLevel();

    /**
     * Set #mHypothesis.
     *
	 * @param hypothesis.
	 */
	void SetHypothesis(unsigned hypothesis);

    /**
     * @return #mHypothesis.
     */
    unsigned GetHypothesis() const; // this function promises not to change the object

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(IngeWntSwatCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a IngeWntSwatCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const IngeWntSwatCellCycleModel * t, const unsigned int file_version)
{
    const unsigned hypothesis = t->GetHypothesis();
    ar & hypothesis;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a IngeWntSwatCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, IngeWntSwatCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of IngeWntSwatCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    std::vector<double> state_vars;
    for (unsigned i=0; i<22; i++)
    {
        state_vars.push_back(0.0);
    }

    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    unsigned dimension = UINT_MAX;

    unsigned hypothesis;
    ar & hypothesis;

    ::new(t)IngeWntSwatCellCycleModel(hypothesis, state_vars, p_mutation_state, dimension);
}
}
} // namespace ...

#endif /*INGEWNTSWATCELLCYCLEMODEL_HPP_*/

