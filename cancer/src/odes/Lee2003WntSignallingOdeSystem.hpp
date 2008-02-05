#ifndef _LEE2003ODESYSTEM_HPP_
#define _LEE2003ODESYSTEM_HPP_

#include <vector>
#include <cmath>

#include "AbstractOdeSystem.hpp"

/**
 * Represents the Lee et al. (2003) Wnt Signalling pathway model equations
 *
 * The variables are
 * % 1. X2 Dsh_active
 * % 2. X3 APC* /axin* /GSK3
 * % 3. X4 APC/axin/GSK3
 * % 4. X9 beta-cat* /APC* /axin* /GSK3
 * % 5. X10 beta-cat*
 * % 6. X11 beta-cat
 * % 7. X12 axin
 * % 8. WntLevel
 *
 */
class Lee2003WntSignallingOdeSystem: public AbstractOdeSystem
{
private:
    // Constants for the Lee et al. (2003) Model
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
    // Constructor
    Lee2003WntSignallingOdeSystem(double WntStimulus = 0.0);
    
    // Destructor
    ~Lee2003WntSignallingOdeSystem();
    
    void Init(); //Sets up the parameter values
    
    // Compute the RHS of the WntCellCycle system of ODEs
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

};
#endif //_LEE2003ODESYSTEM_HPP_
