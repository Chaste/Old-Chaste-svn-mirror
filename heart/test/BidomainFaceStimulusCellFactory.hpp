/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef BIDOMAINFACESTIMULUSCELLFACTORY_HPP_
#define BIDOMAINFACESTIMULUSCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"

class BidomainFaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
    RegularStimulus *mpRegStimulus;
    
public:
    //Pdetime step is (by default) 0.01
    //Odetime step set below to 0.001 (10:1)
    BidomainFaceStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.001)
    {
        mpRegStimulus = new RegularStimulus(-900.0*1000, 0.5, 100.0, 0.0);//Same as above, but every 100ms
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpRegStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
        }
    }
    
    ~BidomainFaceStimulusCellFactory(void)
    {
        delete mpRegStimulus;
    }
};

#endif /*BIDOMAINFACESTIMULUSCELLFACTORY_HPP_*/
