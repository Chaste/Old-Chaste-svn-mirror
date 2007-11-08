#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <time.h>
#include <iostream>
#include <string>

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
    }
    static void PrintAndReset(std::string message)
    {
        Print(message);
        Reset();
    }
};

time_t Timer::StartTime;

#endif /*TIMER_HPP_*/
