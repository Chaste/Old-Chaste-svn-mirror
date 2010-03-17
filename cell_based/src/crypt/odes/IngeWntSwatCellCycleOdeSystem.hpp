/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef INGEWNTSWATCELLCYCLEODESYSTEM_HPP_
#define INGEWNTSWATCELLCYCLEODESYSTEM_HPP_

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"
#include "AbstractCellMutationState.hpp"
#include "UblasCustomFunctions.hpp"

/**
 * Represents the van Leeuwen et al. (2007) system of ODEs
 * [doi:10.1016/j.jtbi.2007.01.019]
 * coupled to the Swat et al. cell cycle model equations.
 * [doi:10.1093/bioinformatics/bth110]
 *
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

    /** Dimensional parameter k_2. */
    double mk2d;
    /** Dimensional parameter k_3. */
    double mk3d;
    /** Dimensional parameter k_34. */
    double mk34d;
    /** Dimensional parameter k_2. */
    double mk43d;
    /** Dimensional parameter k_23. */
    double mk23d;
    /** Dimensional parameter a. */
    double mad;
    /** Dimensional parameter J_11. */
    double mJ11d;
    /** Dimensional parameter J_12. */
    double mJ12d;
    /** Dimensional parameter J_13. */
    double mJ13d;
    /** Dimensional parameter J_13. */
    double mJ61d;
    /** Dimensional parameter J_62. */
    double mJ62d;
    /** Dimensional parameter J_63. */
    double mJ63d;
    /** Dimensional parameter K_m1. */
    double mKm1d;
    /** Dimensional parameter k_p. */
    double mkpd;
    /** Dimensionless parameter phi_r. */
    double mphi_r;
    /** Dimensionless parameter phi_i. */
    double mphi_i;
    /** Dimensionless parameter phi_j. */
    double mphi_j;
    /** Dimensionless parameter phi_p. */
    double mphi_p;
    /** Dimensional parameter k_16. */
    double mk16d;
    /** Dimensional parameter k_61. */
    double mk61d;
    /** Dimensionless parameter phi_E2F1. */
    double mPhiE2F1;

    /**
     * Parameters for the Van Leeuwen et al. (2007) model
     */

    /** Dimensionless parameter s_A. */
    double mSa;
    /** Dimensionless parameter s_CA. */
    double mSca;
    /** Dimensionless parameter s_C. */
    double mSc;
    /** Dimensionless parameter s_CT. */
    double mSct;
    /** Dimensionless parameter s_D. */
    double mSd;
    /** Dimensionless parameter s_T. */
    double mSt;
    /** Dimensionless parameter s_X. */
    double mSx;
    /** Dimensionless parameter s_Y. */
    double mSy;
    /** Dimensionless parameter d_A. */
    double mDa;
    /** Dimensionless parameter d_CA. */
    double mDca;
    /** Dimensionless parameter d_C. */
    double mDc;
    /** Dimensionless parameter d_CT. */
    double mDct;
    /** Dimensionless parameter d_D. */
    double mDd;
    /** Dimensionless parameter d_Dx. */
    double mDdx;
    /** Dimensionless parameter d_T. */
    double mDt;
    /** Dimensionless parameter d_U. */
    double mDu;
    /** Dimensionless parameter d_X. */
    double mDx;
    /** Dimensionless parameter d_Y. */
    double mDy;
    /** Dimensionless parameter K_c. */
    double mKc;
    /** Dimensionless parameter K_D. */
    double mKd;
    /** Dimensionless parameter K_T. */
    double mKt;
    /** Dimensionless parameter p_c. */
    double mPc;
    /** Dimensionless parameter p_u. */
    double mPu;
    /** Dimensionless parameter xi_D. */
    double mXiD;
    /** Dimensionless parameter xi_Dx. */
    double mXiDx;
    /** Dimensionless parameter xi_X. */
    double mXiX;
    /** Dimensionless parameter xi_C. */
    double mXiC;

    /**
     * The mutation state of the cell
     */
    boost::shared_ptr<AbstractCellMutationState> mpMutationState;

    /**
     * The hypothesis we are using
     *  = 1u for Van Leeuwen Hypothesis I
     *  = 2u for Van Leeuwen Hypothesis II
     */
    unsigned mHypothesis;

public:

    /**
     * Constructor.
     *
     * @param hypothesis takes the value 1 or 2 and affects the ODE system.
     * @param wntLevel is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
     * @param pMutationState cell mutation; some affect the ODE system
     */
    IngeWntSwatCellCycleOdeSystem(unsigned hypothesis,
    		                      double wntLevel = 0.0,
    		                      boost::shared_ptr<AbstractCellMutationState> pMutationState=boost::shared_ptr<AbstractCellMutationState>());

    /**
     * Destructor.
     */
    ~IngeWntSwatCellCycleOdeSystem();

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
     * @param pMutationState the mutation state.
     */
    void SetMutationState(boost::shared_ptr<AbstractCellMutationState> pMutationState);

    /**
     * Called by the archive function on the Wnt cell cycle model.
     *
     * @return #mpMutationState the mutation state of the cell.
     */
    boost::shared_ptr<AbstractCellMutationState> GetMutationState();

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
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    /**
     * Calculate whether the conditions for the cell cycle to finish have been met.
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS
     *
     * @return whether or not stopping conditions have been met
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    /**
     * When using CVODE this function is called instead of CalculateStoppingEvent.
     * It allows the point at which rY[1] reaches 1 to be found to greater precision.
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS
     *
     * @return How close we are to the root of the stopping condition
     */
    double CalculateRootFunction(double time, const std::vector<double>& rY);

};

#endif /*INGEWNTSWATCELLCYCLEODESYSTEM_HPP_*/
