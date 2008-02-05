#ifndef _INGEWNTSWATCELLCYCLEODESYSTEM_HPP_
#define _INGEWNTSWATCELLCYCLEODESYSTEM_HPP_

#include <vector>
#include <cmath>

#include "AbstractOdeSystem.hpp"
#include "CellMutationStates.hpp"

/**
 * Represents the van Leeuwen et al. (2007) system of ODEs
 * coupled to the Swat et al. Cell cycle model equations.
 * doi:10.1016/j.jtbi.2007.01.019
 *
 * The variables are
 *
    % 0. r = pRb
    % 1. e = E2F1 (This is the S-phase indicator)
    % 2. i = CycD (inactive)
    % 3. j = CycD (active)
    % 4. p = pRb-p
    % 5. D = APC destruction complex
    % 6. X = Axin
    % 7. Cu = Beta Cat marked for ubiquitination
    % 8. Co = Open form Beta Cat
    % 9. Cc = Closed form Beta Cat
    % 10. Mo = Open form Mutant Beta Cat
    % 11. Mc = Closed form Mutant Beta Cat
    % 12. A = Free Adhesion molecules
    % 13. Ca = BetaCat/Adhesion
    % 14. Ma = Mutant BetaCat/Adhesion
    % 15. T = free TCF
    % 16. Cot = Open BetaCat/TCF
    % 17. Cct = Closed BetaCat/TCF
    % 18. Mot = Open Mutant BetaCat/TCF
    % 19. Mct = Closed Mutant BetaCat/TCF
    % 20. Y = Wnt Target protein
    % 21. Wnt level
 *
 */
class IngeWntSwatCellCycleOdeSystem : public AbstractOdeSystem
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
    double mk16d;
    double mk61d;
    double mPhiE2F1;    
    
    // Inge's Constants 
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
    unsigned mHypothesis;
        
public:
    // Constructor
    IngeWntSwatCellCycleOdeSystem(unsigned hypothesis, double WntStimulus = 0.0, const CellMutationState& rMutationState = HEALTHY);
    
    // Destructor
    ~IngeWntSwatCellCycleOdeSystem();
    
    void Init(); //Sets up the parameter values
    
    void SetMutationState(const CellMutationState &rMutationState);
    
    /** 
     * Called by the archive function on the wnt cell cycle model.
     * @return mMutationState the mutation state of the cell defined by 
     * CryptCellMutationStates.hpp
     */
    CellMutationState& rGetMutationState()
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
#endif //_INGEWNTSWATCELLCYCLEODESYSTEM_HPP_
