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

#include "SumStimulus.hpp"

SumStimulus::SumStimulus(AbstractStimulusFunction* pStimulus1, AbstractStimulusFunction* pStimulus2)
  : mpStimulus1(pStimulus1), mpStimulus2(pStimulus2)
{
}

double SumStimulus::GetStimulus(double time)
{
    return mpStimulus1->GetStimulus(time)+mpStimulus2->GetStimulus(time);
}
