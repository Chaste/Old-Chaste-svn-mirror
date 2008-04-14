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

#ifndef STIMULUSBOUNDARYCONDITION_HPP_
#define STIMULUSBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"
#include "AbstractStimulusFunction.hpp"
#include "PdeSimulationTime.hpp"

/**
 * Boundary condition defined by an AbstractStimlus object.
 */
template<unsigned SPACE_DIM>
class StimulusBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{
private:
    AbstractStimulusFunction* mpStimulus;
    
public:
    /**
     * Create a new boundary condition object.
     * 
     * @param pStimulus Stimulus object defining the parameters of the boundary condition
     */
    StimulusBoundaryCondition(AbstractStimulusFunction* pStimulus)
    {
        mpStimulus = pStimulus;
    }
    
    /**
     * @param x The point at which this boundary condition is to be evaluated.
     * @return The constant value given in the constructor.
     */
    double GetValue( const ChastePoint<SPACE_DIM>& ) const
    {
        return mpStimulus->GetStimulus(PdeSimulationTime::GetTime());
    }
};

#endif /*STIMULUSBOUNDARYCONDITION_HPP_*/
