#ifndef CANCEREVENTHANDLER_HPP_
#define CANCEREVENTHANDLER_HPP_

#include <cassert>
//#include <petsc.h>
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
    const static unsigned NUM_EVENTS=9;
    const static char* EVENT_NAME[NUM_EVENTS];
    //static PetscEvent mPetscEvent[NUM_EVENTS];//Unused?
    static double mCpuTime[NUM_EVENTS];
    
    static void BeginEvent(CancerEventType event)
    {
        mCpuTime[event]-= clock()/(CLOCKS_PER_SEC/1000.0); 
        //std::cout << "Begining " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }
    
    static void EndEvent(CancerEventType event)
    {
        mCpuTime[event]+= clock()/(CLOCKS_PER_SEC/1000.0);
        //std::cout << "Ending " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }
    
    
    static void Report()
    {
        // times are in milliseconds
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            printf("%2.1e\t", mCpuTime[event]);
            mCpuTime[event]=0.0;
        }
        std::cout << "(milliseconds) \n";
    }
    
    static void Headings()
    {
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            printf("%6s\t", EVENT_NAME[event]);
        } 
        std::cout << "\n";
    }
};

#endif /*CANCEREVENTHANDLER_HPP_*/
