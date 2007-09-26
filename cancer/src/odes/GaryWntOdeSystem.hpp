#ifndef _GARYWNTODESYSTEM_HPP_
#define _GARYWNTODESYSTEM_HPP_

#include <vector>
#include <cmath>
#include "AbstractOdeSystem.hpp"
#include "CryptCellMutationStates.hpp"

/**
 * Represents the Mirams et al. (2007-with any luck) system of ODEs.
 * Based on Swat et al. (2004) and a simple Wnt model.
 *
 * The variables are
 *
 1. c = destruction complex (Active)
 2. b1 = Beta-Catenin (from 1st allele)
 3. b2 = Beta-Catenin (from 1st allele)
 4. WntLevel
 5. mutation state (0/1/2/3)
 *
 */
class GaryWntOdeSystem : public AbstractOdeSystem
{
private:
    double ma2d;
    double ma3d;
    double ma4d;
    double ma5d;
    double mPhiE2F1;  // TODO: remove
    CryptCellMutationState mMutationState;
    
public:
    // Constructor
    GaryWntOdeSystem(double WntStimulus = 0.0, const CryptCellMutationState& rMutationState = HEALTHY);
    
    // Destructor
    ~GaryWntOdeSystem();
    
    void Init(); //Sets up the parameter values
    
    void SetMutationState(const CryptCellMutationState &rMutationState);
    
    /** 
     * Called by the archive function on the wnt cell cycle model.
     * @return mMutationState the mutation state of the cell defined by 
     * CryptCellMutationStates.hpp
     */
    CryptCellMutationState& rGetMutationState()
    {
        return mMutationState;
    }
    
    // Compute the RHS of the WntCellCycle system of ODEs
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {
        std::vector<double> dy(rY.size());
        EvaluateYDerivatives(time, rY, dy);
        return (fabs(rY[1]-1.0) < 1.0e-2 && dy[1] > 0.0);
    }
    
};
#endif //_GARYWNTODESYSTEM_HPP_
