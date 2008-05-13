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
#ifndef STOCHASTICWNTCELLCYCLEMODEL_HPP_
#define STOCHASTICWNTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "WntCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

class StochasticWntCellCycleModel : public WntCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<WntCellCycleModel>(*this);
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
        
        archive & mG2Duration;
    }
    
    /**
     * This is a function which overrides that in WntCellCycleModel and 
     * introduces the stochastic element of this class. 
     * We allow the duration of the G2 phase of the cell cycle to 
     * vary with a mean of its deterministic duration and a standard 
     * deviation of 0.9 hours.
     * 
     * @return the duration of the G2 phases of the cell cycle.
     */
    void SetG2Duration();
    
    /// The duration of the G2 phase, set stochastically
    double mG2Duration;
        
public:

    void InitialiseDaughterCell();

    void Initialise();
    
    void ResetForDivision();
    
    double GetG2Duration();

    /**
     * The standard constructor called in tests
     */
    StochasticWntCellCycleModel()
      :  WntCellCycleModel() {};
    
    /**
     * This is needed because a wnt model which is not to be run from the current time is 
     * sometimes needed. Should only be called by the cell itself when it wants to divide.
     */
    StochasticWntCellCycleModel(AbstractOdeSystem* pParentOdeSystem,
                                CellMutationState mutationState,
                                double birthTime, double lastTime,
                                bool inSG2MPhase, bool readyToDivide, double divideTime, unsigned generation, double g2Duration)
      : WntCellCycleModel(pParentOdeSystem, mutationState, birthTime, lastTime, 
                          inSG2MPhase, readyToDivide, divideTime, generation),
        mG2Duration(g2Duration) {};

    /**
     * This is needed because a wnt model which is not to be run from the current time is 
     * sometimes needed. Should only be called by the archiver.
     */
    StochasticWntCellCycleModel(std::vector<double> proteinConcentrations, 
                                CellMutationState mutationState)
      : WntCellCycleModel(proteinConcentrations, mutationState) {};
    
    /**
     * Returns a new StochasticWntCellCycleModel created with the correct initial conditions.
     *
     * Should be called just after the parent cell cycle model has been .Reset().
     *
     */
    AbstractCellCycleModel* CreateDaughterCellCycleModel();
    
};


// declare identifier for the serializer
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
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of my_class

    std::vector<double> state_vars;
    for (unsigned i=0; i<9; i++)
    {
        state_vars.push_back(0.0);
    }   

    CellMutationState mutation_state = HEALTHY;

    ::new(t)StochasticWntCellCycleModel(state_vars, mutation_state);
}
}
} // namespace ...

#endif /*STOCHASTICWNTCELLCYCLEMODEL_HPP_*/
