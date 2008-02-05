#ifndef _ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
#define _ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_

#include <vector>
#include <cmath>

#include "AbstractOdeSystem.hpp"
#include "CellMutationStates.hpp"

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
    CellMutationState mMutationState;
 
        
public:
    // Constructor
    Alarcon2004OxygenBasedCellCycleOdeSystem(double oxygenConcentration, const CellMutationState& rMutationState);  
    
    // Destructor
    ~Alarcon2004OxygenBasedCellCycleOdeSystem();
    
    void Init(); // sets up the parameter values
    
    void SetMutationState(const CellMutationState &rMutationState);
    
    /** 
     * Called by the archive function on the cell cycle model.
     * @return mMutationState the mutation state of the cell defined by 
     * CellMutationStates.hpp
     */
    CellMutationState& rGetMutationState()
    {
        return mMutationState;
    }
    
    // Compute the RHS of the Alarcon et al. (2004) system of ODEs
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {        
        return (rY[0] < mxThreshold && rY[1] > myThreshold);
    }
    
};
#endif //_ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
