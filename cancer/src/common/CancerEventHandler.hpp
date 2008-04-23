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

#ifndef CANCEREVENTHANDLER_HPP_
#define CANCEREVENTHANDLER_HPP_

#include "GenericEventHandler.hpp"

typedef enum CancerEventType_
{
    SETUP=0,
    DEATH,
    BIRTH,
    REMESH,
    TESSELLATION,
    VELOCITY,
    POSITION,
    OUTPUT,
    CANCER_EVERYTHING
} CancerEventType;


class CancerEventNames
{
    public:
    const static char* EVENT_NAME[9];
};

typedef GenericEventHandler<9, CancerEventNames::EVENT_NAME> CancerEventHandler;

#endif /*CANCEREVENTHANDLER_HPP_*/
