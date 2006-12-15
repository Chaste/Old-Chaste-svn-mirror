#ifndef _WNTCELLCYCLEODESYSTEM_HPP_
#define _WNTCELLCYCLEODESYSTEM_HPP_

#include <vector>
#include <cmath>
#include "AbstractOdeSystem.hpp"

/**
 * Represents the Mirams et al. (2007-with any luck) system of ODEs.
 * Based on Swat et al. (2004) and a simple Wnt model.
 * 
 * The variables are
 * 
 1. r = pRb
 2. e = E2F1
 3. i = CycD (inactive)
 4. j = CycD (active)
 5. p = pRb-p
 6. c = destruction complex (Active)
 7. b = Beta-Catenin
 8. WntLevel 
 * 
 */
class WntCellCycleOdeSystem: public AbstractOdeSystem
{
private:
	// Constants for the Swat et al. (2004) Model
    double mk2d;
	double mk3d;
	double mk34d;
	double mk43d;
	double mk23d;
	double mad;
	double mJ11d;
	double mJ12d;
	double mJ13d;
	double mJ61d;
	double mJ62d;
	double mJ63d;
	double mKm1d;
	double mkpd;
	double mphi_r;
	double mphi_i;
	double mphi_j;
	double mphi_p;
	double ma2d;
	double ma3d;
	double ma4d;
	double ma5d;
	double mk16d;
	double mk61d;
	double mPhiE2F1;
    
public:
    // Constructor
    WntCellCycleOdeSystem();
    
    // Destructor
    ~WntCellCycleOdeSystem();
    
    void Init(); //Sets up the parameter values
        
    // Compute the RHS of the T&N system of ODEs
    std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
    
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {
         return (fabs(rY[1]-1.0) < 1.0e-2 && EvaluateYDerivatives(time, rY)[1] > 0.0);
    }
    
};
#endif //_WNTCELLCYCLEODESYSTEM_HPP_
