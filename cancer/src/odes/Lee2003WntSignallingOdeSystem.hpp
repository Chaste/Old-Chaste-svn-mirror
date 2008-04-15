#ifndef _LEE2003ODESYSTEM_HPP_
#define _LEE2003ODESYSTEM_HPP_

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

/**
 * Represents the Lee et al. (2003) Wnt Signalling pathway model equations
 *
 * The variables are
 *   1. X2 Dsh_active
 *   2. X3 APC/axin/GSK3
 *   3. X4 APC/axin/GSK3
 *   4. X9 beta-cat/APC/axin/GSK3
 *   5. X10 beta-cat
 *   6. X11 beta-cat
 *   7. X12 axin
 *   8. WntLevel
 */
class Lee2003WntSignallingOdeSystem: public AbstractOdeSystem
{
private:

    /**
     * Parameters for the Lee et al. (2003) Model
     */ 
    double mDsh0;
    double mAPC0;
    double mTCF0;
    double mGSK0;
    double mK7;
    double mK8;
    double mK16;
    double mK17;
    double mk1;
    double mk2;
    double mk3;
    double mk4;
    double mk5;
    double mk6;
    double mk_6;
    double mk9;
    double mk10;
    double mk11;
    double mk13;
    double mk15;
    double mv12;
    double mv14;

public:
    
    /**
     * Constructor.
     *
     * @param WntLevel is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
     * 
     */
    Lee2003WntSignallingOdeSystem(double WntStimulus=0.0);
    
    /**
     * Destructor.
     */
    ~Lee2003WntSignallingOdeSystem();
    
    /**
     * Initialise parameter values.
     */
    void Init();
    
    /**
     * Compute the RHS of the Lee et al. (2003) system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *  
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using Alarcons et al. (2004) system of equations).
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

};

#endif //_LEE2003ODESYSTEM_HPP_
