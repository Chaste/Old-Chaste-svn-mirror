#ifndef ABSTRACTPDE_HPP_
#define ABSTRACTPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include <petscvec.h>

class AbstractPde
{
public :
    /**
    * PrepareForAssembleSystem is a virtual method.
    * It's called by the AssembleSystem method of the assembler before any other
    * useful work happens.  The idea is that a *coupled system* will want to 
    * solve all the ODE systems before the PDE is solved.  A *parallel* coupled
    * system will want to solve the ODE systems and distribute the answers 
    * before anything else happens.
    */
    virtual void PrepareForAssembleSystem(Vec currentSolution, double /*time - not used here*/)
    {
    }

    virtual ~AbstractPde()
    {
    }
};

#endif /*ABSTRACTPDE_HPP_*/
