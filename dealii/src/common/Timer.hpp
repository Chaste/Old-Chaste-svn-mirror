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

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <time.h>
#include <iostream>
#include <string>
#include "LogFile.hpp"

class Timer
{
private:
    static time_t StartTime;

public:
    static void Reset()
    {
        StartTime = std::clock();
    }
    static void Print(std::string message)
    {
        std::cout << message << " time is "  
                  << (std::clock() - StartTime)/(CLOCKS_PER_SEC) 
                  << "s\n" << std::flush;
        LOG(1,"    "<<message << " time is "<<(std::clock()-StartTime)/(CLOCKS_PER_SEC)<<"s");                  
    }
    static void PrintAndReset(std::string message)
    {
        Print(message);
        Reset();
    }
};

time_t Timer::StartTime;

#endif /*TIMER_HPP_*/
