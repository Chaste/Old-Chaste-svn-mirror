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

#ifndef MULTISTIMULUS_HPP_
#define MULTISTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class provides a stimulus function which is the 
 * sum of an arbitrary number of stimuli.
 * 
 * After creation it behaves like a ZeroStimulus until
 * any number of stimuli are added.
 */

class MultiStimulus : public AbstractStimulusFunction
{
private:
    std::vector<AbstractStimulusFunction*> mStimuli;
        
public:
    /**
     * Combine a stimulus with the existing ones.
     * 
     * @param pStimulus pointer to the stimulus to be added.
     */   
     void AddStimulus(AbstractStimulusFunction* pStimulus);

    /**
     * Get the magnitude of the multiple stimuli at time 'time'
     *
     * @return  Magnitude of stimulus at time 'time'.
     */
     double GetStimulus(double time);
};

#endif /*MULTISTIMULUS_HPP_*/
