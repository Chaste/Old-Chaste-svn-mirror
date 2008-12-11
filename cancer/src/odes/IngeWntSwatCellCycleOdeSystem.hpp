/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef _INGEWNTSWATCELLCYCLEODESYSTEM_HPP_
#define _INGEWNTSWATCELLCYCLEODESYSTEM_HPP_

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"
#include "CellMutationStates.hpp"

/**
 * Represents the van Leeuwen et al. (2007) system of ODEs
 * coupled to the Swat et al. cell cycle model equations.
 * doi:10.1016/j.jtbi.2007.01.019
 *
 * The variables are
 *
 *   0. r = pRb
 *   1. e = E2F1 (This is the S-phase indicator)
 *   2. i = CycD (inactive)
 *   3. j = CycD (active)
 *   4. p = pRb-p
 *   5. D = APC destruction complex
 *   6. X = Axin
 *   7. Cu = Beta Cat marked for ubiquitination
 *   8. Co = Open form Beta Cat
 *   9. Cc = Closed form Beta Cat
 *   10. Mo = Open form Mutant Beta Cat
 *   11. Mc = Closed form Mutant Beta Cat
 *   12. A = Free Adhesion molecules
 *   13. Ca = BetaCat/Adhesion
 *   14. Ma = Mutant BetaCat/Adhesion
 *   15. T = free TCF
 *   16. Cot = Open BetaCat/TCF
 *   17. Cct = Closed BetaCat/TCF
 *   18. Mot = Open Mutant BetaCat/TCF
 *   19. Mct = Closed Mutant BetaCat/TCF
 *   20. Y = Wnt Target protein
 *   21. Wnt level
 */
class IngeWntSwatCellCycleOdeSystem : public AbstractOdeSystem
{
private:

    /**
     * Parameters for the Swat et al. (2004) model
     */
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

    /**
     * Parameters for the Van Leeuwen et al. (2007) model
     */
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

    /**
     * Constructor.
     *
     * @param hypothesis takes the value 1 or 2 and affects the ODE system.
     * @param wntLevel is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
     * @param rMutationState affects the ODE system and is given by CryptCellMutationStates.hpp
     */
    IngeWntSwatCellCycleOdeSystem(unsigned hypothesis, double WntStimulus = 0.0, const CellMutationState& rMutationState = HEALTHY);

    /**
     * Destructor
     */
    ~IngeWntSwatCellCycleOdeSystem();

    /**
     * Initialise parameter values
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
     * Called by the archive function on the Wnt cell cycle model.
     *
     * @return mMutationState the mutation state of the cell defined by
     * CryptCellMutationStates.hpp
     */
    CellMutationState& rGetMutationState();

    /**
     * Compute the RHS of the system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using van Leeuwen et al. (2007) system of equations)
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    /**
     * Calculate whether the conditions for the cell cycle to finish have been met.
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS.
     */
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY);
    
    /**
     * When using CVODE this function is called instead of CalculateStoppingEvent.
     * It allows the point at which rY[1] reaches 1 to be found to greater precision.
     */
    double CalculateRootFunction(double time, const std::vector<double> &rY);

};

#endif //_INGEWNTSWATCELLCYCLEODESYSTEM_HPP_
