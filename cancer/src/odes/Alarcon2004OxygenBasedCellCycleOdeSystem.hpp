#ifndef _ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
#define _ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_

#include <vector>
#include <cmath>
#include "AbstractOdeSystem.hpp"

/**
 * Represents the Alarcon et al. (2004) system of ODEs (see ticket #461).
 *
 * The variables are
 *
 0. x = Cdh1-APC complexes
 1. y = cyclin-CDK
 2. z = p27
 3. m = mass
 4. u = RBNP
 5. P = oxygen concentration
 *
 */
      
class Alarcon2004OxygenBasedCellCycleOdeSystem : public AbstractOdeSystem
{
private:
    // Constants for the Alarcon et al. (2004) model
    double ma1;
    double ma2;
    double ma3;
    double ma4;
    double mb3;
    double mb4;
    double mc1;
    double mc2;
    double md1;
    double md2;
    double mJ3;
    double mJ4;
    double mEta;
    double mMstar;
    double mB;
    double mxThreshold;
    double myThreshold;
    bool mIsCancerCell;
 
        
public:
    // Constructor
    Alarcon2004OxygenBasedCellCycleOdeSystem(double oxygenConcentration, const bool &rIsCancerCell);  
    
    // Destructor
    ~Alarcon2004OxygenBasedCellCycleOdeSystem();
    
    void Init(); // sets up the parameter values
    
    void SetIsCancerCell(const bool &rIsCancerCell);
    
//    /** 
//     * Called by the archive function on the Alarcon et al. (2004) model.
//     * @return IsCancerCell whether the cell is a 'cancer cell' (in the sense of the Alarcon et al. (2004) model)
//     */
//    bool& rGetIsCancerCell()
//    {
//        return mIsCancerCell;
//    }
    
    // Compute the RHS of the Alarcon et al. (2004) system of ODEs
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {        
        return (rY[0] < mxThreshold && rY[1] > myThreshold);
    }
    
};
#endif //_ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
