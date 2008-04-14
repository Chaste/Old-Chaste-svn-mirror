#ifndef CANCEREVENTHANDLER_HPP_
#define CANCEREVENTHANDLER_HPP_

#include <cassert>
#include <time.h>
#include <iostream>

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


class CancerEventHandler
{
public:
    
    const static unsigned NUM_EVENTS = 9;
    
    const static char* EVENT_NAME[NUM_EVENTS];
    
    static double mCpuTime[NUM_EVENTS];
    
    static void BeginEvent(CancerEventType event);
    
    static void EndEvent(CancerEventType event);
    
    static void Report();
    
    static void Headings();
};

#endif /*CANCEREVENTHANDLER_HPP_*/
