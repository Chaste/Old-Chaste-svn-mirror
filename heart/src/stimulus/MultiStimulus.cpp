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

#include "MultiStimulus.hpp"

void MultiStimulus::AddStimulus(AbstractStimulusFunction* pStimulus)
{
    mStimuli.push_back(pStimulus);
}

double MultiStimulus::GetStimulus(double time)
{
    double total_stimulus = 0.0;
    
    for (unsigned current_stimulus = 0; current_stimulus < mStimuli.size(); ++current_stimulus)
    {
        total_stimulus += mStimuli[current_stimulus]->GetStimulus(time);
    }
    
    return total_stimulus;
}
