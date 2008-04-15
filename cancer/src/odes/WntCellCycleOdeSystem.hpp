#ifndef _WNTCELLCYCLEODESYSTEM_HPP_
#define _WNTCELLCYCLEODESYSTEM_HPP_

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"
#include "CellMutationStates.hpp"

/**
 * Represents the Mirams et al. system of ODEs, based on Swat et al. (2004) 
 * and a simple Wnt model.
 *
 * The variables are
 *
 * 0. r = pRb
 * 1. e = E2F1
 * 2. i = CycD (inactive)
 * 3. j = CycD (active)
 * 4. p = pRb-p
 * 5. c = destruction complex (Active)
 * 6. b1 = Beta-Catenin (from 1st allele)
 * 7. b2 = Beta-Catenin (from 1st allele)
 * 8. WntLevel 
 */
class WntCellCycleOdeSystem : public AbstractOdeSystem
{
private:

    // Parameters for the Swat et al. (2004) Model
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
        
    CellMutationState mMutationState;
        
public:
        
    /**
     * Constructor.
     *
     * @param WntStimulus is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
     * @param rMutationState affects the ODE system and is given by CellMutationStates.hpp
     */
    WntCellCycleOdeSystem(double WntStimulus=0.0, const CellMutationState& rMutationState=HEALTHY);
    
    /**
     * Destructor.
     */
    ~WntCellCycleOdeSystem();
    
    /**
     * Initialise parameter values.
     */
    void Init();
    
    /**
     * Set the mutation state of the cell.
     * 
     * This should be called by the relevant cell cycle model before any solving
     * of the ODE system (as it is used to evaluate the Y derivatives).
     * 
     * @param rMutationState the mutation state.
     */
    void SetMutationState(const CellMutationState &rMutationState);
    
    /** 
     * Called by the archive function on the wnt cell cycle model.
     * 
     * @return mMutationState the mutation state of the cell defined by 
     * CellMutationStates.hpp
     */
    CellMutationState& rGetMutationState();
    
    /**
     * Compute the RHS of the WntCellCycle system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *  
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using Alarcons et al. (2004) system of equations).
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    /**
     * This also contains a calculation of dY[1], copied from EvaluateYDerivatives.
     * Ensure they do not get out of sync!
     */
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY);
    
};

#endif //_WNTCELLCYCLEODESYSTEM_HPP_
