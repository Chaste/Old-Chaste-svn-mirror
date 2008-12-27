/*

Copyright (C) University of Oxford, 2008

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

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <cfloat>

#include "AbstractWntOdeBasedCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "WntConcentration.hpp"
#include "CellMutationStates.hpp"


// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * Wnt-dependent cell cycle model.
 *
 * Note that this class uses C++'s default copying semantics, and so 
 * doesn't implement a copy constructor or operator=.
 */
class WntCellCycleModel : public AbstractWntOdeBasedCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        assert(mpOdeSystem!=NULL);
        archive & boost::serialization::base_object<AbstractWntOdeBasedCellCycleModel>(*this);
        /**
         * Reference can be read or written into once mpOdeSystem has been set up
         * mpOdeSystem isn't set up by the first constructor, but is by the second 
         * which is now utilised by the load_construct at the bottom of this file.
         */
        archive & static_cast<WntCellCycleOdeSystem*>(mpOdeSystem)->rGetMutationState();
    }

    /**
     * Update the cell type according to the current beta catenin 
     * level as given by the WntCellCycleOdeSystem.
     * 
     * This method carries out the work for UpdateCellType(), but 
     * does not check the current time, so can also be called by 
     * Initialise().
     */
    void ChangeCellTypeDueToCurrentBetaCateninLevel();

public:

    /**
     * Default constructor.
     */
    WntCellCycleModel() {};

    /**
     * A private constructor for daughter cells called by the CreateDaughterCellCycleModel function
     * (which can be called by TissueCell::CommonCopy() and isn't necessarily being born.
     *
     * @param pParentOdeSystem  to copy the state of
     * @param rMutationState the mutation state of the cell (used by ODEs)
     * @param birthTime the simulation time when the cell divided (birth time of parent cell)
     * @param lastTime last time the cell cycle model was evaluated
     * @param inSG2MPhase whether the cell is in S-G2-M (not evaluating ODEs and just waiting)
     * @param readyToDivide whether the cell is ready to divide
     * @param divideTime if in the future this is the time at which the cell is going to divide
     * @param generation the cell's generation
     */
    WntCellCycleModel(AbstractOdeSystem* pParentOdeSystem,
                      const CellMutationState& rMutationState,
                      double birthTime,
                      double lastTime,
                      bool inSG2MPhase,
                      bool readyToDivide,
                      double divideTime,
                      unsigned generation);

    /**
     * A private constructor for archiving.
     *
     * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
     * @param rMutationState the mutation state of the cell (used by ODEs)
     */
    WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                      const CellMutationState& rMutationState);

    /**
     * Returns a new WntCellCycleModel, created with the correct 
     * initial conditions.
     *
     * This method should be called just after the parent cell cycle model 
     * has been reset.
     * 
     * @return pointer to the daughter cell cycle model
     */
    AbstractCellCycleModel* CreateDaughterCellCycleModel();

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


// Declare identifier for the serializer
BOOST_CLASS_EXPORT(WntCellCycleModel)

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

    CellMutationState mutation_state = HEALTHY;
    ::new(t)WntCellCycleModel(state_vars, mutation_state);
}
}
} // namespace

#endif /*WNTCELLCYCLEMODEL_HPP_*/
