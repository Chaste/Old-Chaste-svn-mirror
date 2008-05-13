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
#ifndef ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_
#define ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include <cfloat>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CellwiseData.hpp"
#include "Exception.hpp"


// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 *  Oxygen-dependent cell cycle model.
 * 
 * Note that this class uses C++'s default copying semantics, and so 
 * doesn't implement a copy constructor or operator=.
 * 
 * Note also that this model currently only works in 2D, since the 
 * SolveOdeToTime() and GetDivideTime() methods involve instances of 
 * CellwiseData<2>. 
 */
class Alarcon2004OxygenBasedCellCycleModel : public AbstractOdeBasedCellCycleModel
{
    friend class boost::serialization::access;   
    
private:
    static RungeKutta4IvpOdeSolver msSolver;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        assert(mpOdeSystem!=NULL); 
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
        // reference can be read or written into once mpOdeSystem has been set up
        // mpOdeSystem isn't set up by the first constructor, but is by the second
        // which is now utilised by the load_construct at the bottom of this file.
        archive & static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->rGetMutationState(); 
    }
        
public:

    /**
     * Default constructor, variables are set by abstract classes.
     */
    Alarcon2004OxygenBasedCellCycleModel() {};
   

    Alarcon2004OxygenBasedCellCycleModel(AbstractOdeSystem* pParentOdeSystem, 
                                         const CellMutationState& rMutationState, 
                                         double birthTime, 
                                         double lastTime, 
                                         bool inSG2MPhase, 
                                         bool readyToDivide, 
                                         double divideTime, 
                                         unsigned generation);

    Alarcon2004OxygenBasedCellCycleModel(const std::vector<double>& rParentProteinConcentrations, 
                                         const CellMutationState& rMutationState); 
                          
    virtual void ResetForDivision();
    
    AbstractCellCycleModel *CreateDaughterCellCycleModel();
    
    void Initialise();    
    
    bool SolveOdeToTime(double currentTime);
    
    double GetOdeStopTime();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of my_class   
    
    std::vector<double> state_vars;
    for (unsigned i=0; i<6; i++)
    {
        state_vars.push_back(0.0);
    }
    ::new(t)Alarcon2004OxygenBasedCellCycleModel(state_vars, HEALTHY);
}
}
} // namespace ...

#endif /*ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_*/
