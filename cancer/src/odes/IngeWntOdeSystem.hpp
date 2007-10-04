#ifndef _INGEWNTODESYSTEM_HPP_
#define _INGEWNTODESYSTEM_HPP_

#include <vector>
#include <cmath>
#include "AbstractOdeSystem.hpp"
#include "CellMutationStates.hpp"

/**
 * Represents the van Leeuwen et al. (2007) system of ODEs.
 * doi:10.1016/j.jtbi.2007.01.019
 * 
 * Taking Hypothesis One Equations by default.
 * 
 */
class IngeWntOdeSystem : public AbstractOdeSystem
{
private:
    // Constants 
    double mSa;
    double mSca;
    double mSc;
    double mSct;
    double mSd;
    double mSt;
    double mSx;
    double mSy;
    double mDa;
    double mDca;
    double mDc;
    double mDct;
    double mDd;
    double mDdx;
    double mDt;
    double mDu;
    double mDx;
    double mDy;
    double mKc;
    double mKd;
    double mKt;
    double mPc;
    double mPu;
    double mXiD;
    double mXiDx;
    double mXiX;
    double mXiC;
    CellMutationState mMutationState;
    
public:
    // Constructor
    IngeWntOdeSystem(double WntStimulus = 0.0, const CellMutationState& rMutationState = HEALTHY);
    
    // Destructor
    ~IngeWntOdeSystem();
    
    void Init(); //Sets up the parameter values
    
    void SetMutationState(const CellMutationState &rMutationState);
    
    CellMutationState& rGetMutationState();

    
    // Compute the RHS of the WntCellCycle system of ODEs
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    void SetUseHypothesisTwo(bool hypothesisTwo);

};
#endif //_INGEWNTODESYSTEM_HPP_
