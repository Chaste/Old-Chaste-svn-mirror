/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PLANESTIMULUSCELLFACTORY_HPP_
#define PLANESTIMULUSCELLFACTORY_HPP_

#include "LuoRudyIModel1991OdeSystem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LogFile.hpp"

template<unsigned DIM>
class PlaneStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PlaneStimulusCellFactory() : AbstractCardiacCellFactory<DIM>(0.01)//Ode timestep
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(-600, 0.5);
        LOG(1, "Defined a PlaneStimulusCellFactory<"<<DIM<<"> with InitialStimulus(-600, 0.5)\n");
    }
    
    PlaneStimulusCellFactory(double timeStep, double stimulusMagnitude) : AbstractCardiacCellFactory<DIM>(timeStep)
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(stimulusMagnitude, 0.5);
        LOG(1, "Defined a PlaneStimulusCellFactory<"<<DIM<<"> with InitialStimulus("<<stimulusMagnitude<<",0.5)\n");
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (this->mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver,
                                                  this->mTimeStep,
                                                  mpStimulus,
                                                  this->mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver,
                                                  this->mTimeStep,
                                                  this->mpZeroStimulus,
                                                  this->mpZeroStimulus);
        }
    }
    
    ~PlaneStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


#endif /*PLANESTIMULUSCELLFACTORY_HPP_*/
