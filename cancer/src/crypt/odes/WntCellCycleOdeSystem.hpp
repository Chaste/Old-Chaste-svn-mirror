/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef WNTCELLCYCLEODESYSTEM_HPP_
#define WNTCELLCYCLEODESYSTEM_HPP_

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"
#include "CryptCellMutationStates.hpp"

/**
 * Represents the Mirams et al. system of ODEs, based on Swat et al. (2004)
 * [doi:10.1093/bioinformatics/bth110]
 * and a simple Wnt model (unpublished)
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
    /** Dimensional parameter a_2. */
    double ma2d;
    /** Dimensional parameter a_3. */
    double ma3d;
    /** Dimensional parameter a_4. */
    double ma4d;
    /** Dimensional parameter a_5. */
    double ma5d;
    /** Dimensional parameter k_16. */
    double mk16d;
    /** Dimensional parameter k_61. */
    double mk61d;
    /** Dimensionless parameter phi_E2F1. */
    double mPhiE2F1;

    /** The mutation state of the cell - Wnt pathway behaviour (and hence cell cycle time) changes depending on this */
    CryptCellMutationState mMutationState;

public:

    /**
     * Constructor.
     *
     * @param WntStimulus is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
     * @param rMutationState affects the ODE system and is given by CryptCellMutationStates.hpp
     */
    WntCellCycleOdeSystem(double WntStimulus=0.0, const CryptCellMutationState& rMutationState=HEALTHY);

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
    void SetMutationState(const CryptCellMutationState& rMutationState);

    /**
     * Called by the archive function on the Wnt cell cycle model.
     *
     * @return mMutationState the mutation state of the cell defined by
     * CryptCellMutationStates.hpp
     */
    CryptCellMutationState& rGetMutationState();

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
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    /**
     * This also contains a calculation of dY[1], copied from EvaluateYDerivatives.
     * Ensure they do not get out of sync!
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS
     *
     * @return whether we have reached the stopping event
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    /**
     * When using CVODE this function is called instead of CalculateStoppingEvent.
     * It allows the point at which rY[1] reaches 1 to be found to greater precision.
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS
     *
     * @return function value - giving CVODE an estimate of how close we are to the root.
     */
    double CalculateRootFunction(double time, const std::vector<double>& rY);

};

#endif /*WNTCELLCYCLEODESYSTEM_HPP_*/
