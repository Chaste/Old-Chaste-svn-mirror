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

#ifndef SUMSTIMULUS_HPP_
#define SUMSTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"

/**
 * Provides a stimulus function which is the sum of two such functions.
 */
class SumStimulus : public AbstractStimulusFunction
{
private:
    AbstractStimulusFunction* mpStimulus1;
    AbstractStimulusFunction* mpStimulus2;
    
public:
    SumStimulus(AbstractStimulusFunction *mStimulus1, AbstractStimulusFunction *mStimulus2);
    double GetStimulus(double time);
    
};
#endif /*SUMSTIMULUS_HPP_*/
