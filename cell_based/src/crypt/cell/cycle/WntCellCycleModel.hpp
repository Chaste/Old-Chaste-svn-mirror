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
#ifndef WNTCELLCYCLEMODEL_HPP_
#define WNTCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cfloat>

#include "AbstractWntOdeBasedCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "AbstractCellMutationState.hpp"


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
class WntCellCycleModel : public AbstractWntOdeBasedCellCycleModel
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
        archive & boost::serialization::base_object<AbstractWntOdeBasedCellCycleModel>(*this);
        boost::shared_ptr<AbstractCellMutationState> p_mutation_state = static_cast<WntCellCycleOdeSystem*>(mpOdeSystem)->GetMutationState();
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
        static_cast<WntCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(p_mutation_state);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /**
     * Update the cell type according to the current beta catenin
     * level as given by the WntCellCycleOdeSystem.
     *
     * This method carries out the work for UpdateCellProliferativeType(), but
     * does not check the current time, so can also be called by
     * Initialise().
     */
    void ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();

public:

    /**
     * Default constructor.
     */
    WntCellCycleModel();

    /**
     * Copy constructor.
     *
     * Also copies our ODE system.
     *
     * @param rOtherModel the instance being copied.
     */
    WntCellCycleModel(const WntCellCycleModel& rOtherModel);

    /**
     * A private constructor for daughter cells called by the CreateDaughterCellCycleModel function.
     *
     * @param pParentOdeSystem  to copy the state of
     * @param pMutationState the mutation state of the cell (used by ODEs)
     * @param birthTime the simulation time when the cell divided (birth time of parent cell)
     * @param lastTime last time the cell cycle model was evaluated
     * @param inSG2MPhase whether the cell is in S-G2-M (not evaluating ODEs and just waiting)
     * @param readyToDivide whether the cell is ready to divide
     * @param divideTime if in the future this is the time at which the cell is going to divide
     * @param dimension the spatial dimension
     */
    WntCellCycleModel(AbstractOdeSystem* pParentOdeSystem,
                      boost::shared_ptr<AbstractCellMutationState> pMutationState,
                      double birthTime,
                      double lastTime,
                      bool inSG2MPhase,
                      bool readyToDivide,
                      double divideTime,
                      unsigned dimension);

    /**
     * A private constructor for archiving.
     *
     * @param rParentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
     * @param pMutationState the mutation state of the cell (used by ODEs)
     * @param rDimension the spatial dimension
     */
    WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                      boost::shared_ptr<AbstractCellMutationState> pMutationState,
                      const unsigned& rDimension);

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This overridden method sets up a new WntCellCycleOdeSystem
     * and sets the cell type according to the current beta catenin
     * level.
     */
    void Initialise();

    /**
     * Solve the ODEs up to the current time and return whether a stopping event occurred.
     *
     * @param currentTime the current time
     * @return whether a stopping event occured
     */
    bool SolveOdeToTime(double currentTime);
};


#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(WntCellCycleModel)

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
    Archive & ar, const WntCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a WntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, WntCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of WntCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    std::vector<double> state_vars;
    for (unsigned i=0; i<9; i++)
    {
        state_vars.push_back(0.0);
    }

    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    unsigned dimension = 1;
    ::new(t)WntCellCycleModel(state_vars, p_mutation_state, dimension);
}
}
} // namespace

#endif /*WNTCELLCYCLEMODEL_HPP_*/
