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
#ifndef LEE2003WNTSIGNALLINGODESYSTEM_HPP_
#define LEE2003WNTSIGNALLINGODESYSTEM_HPP_

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

/**
 * Represents the Lee et al. (2003) Wnt Signalling pathway model equations.
 * [doi:10.1371/journal.pbio.0000010]
 *
 * The variables are
 *   1. X2 Dsh_active
 *   2. X3 APC/axin/GSK3
 *   3. X4 APC/axin/GSK3
 *   4. X9 beta-cat/APC/axin/GSK3
 *   5. X10 beta-cat
 *   6. X11 beta-cat
 *   7. X12 axin
 *   8. WntLevel
 */
class Lee2003WntSignallingOdeSystem: public AbstractOdeSystem
{
private:

    /**
     * Parameters for the Lee et al. (2003) Model
     */

    /** Dimensional parameter Dsh_0. */
    double mDsh0;
    /** Dimensional parameter APC_0. */
    double mAPC0;
    /** Dimensional parameter TCF_0. */
    double mTCF0;
    /** Dimensional parameter GSK_0. */
    double mGSK0;
    /** Dimensional parameter K_7. */
    double mK7;
    /** Dimensional parameter K_8. */
    double mK8;
    /** Dimensional parameter K_16. */
    double mK16;
    /** Dimensional parameter K_17. */
    double mK17;
    /** Dimensional parameter k_1. */
    double mk1;
    /** Dimensional parameter k_2. */
    double mk2;
    /** Dimensional parameter k_3. */
    double mk3;
    /** Dimensional parameter k_4. */
    double mk4;
    /** Dimensional parameter k_5. */
    double mk5;
    /** Dimensional parameter k_6. */
    double mk6;
    /** Dimensional parameter k_-6. */
    double mk_6;
    /** Dimensional parameter k_9. */
    double mk9;
    /** Dimensional parameter k_10. */
    double mk10;
    /** Dimensional parameter k_11. */
    double mk11;
    /** Dimensional parameter k_13. */
    double mk13;
    /** Dimensional parameter k_15. */
    double mk15;
    /** Dimensional parameter v_12. */
    double mv12;
    /** Dimensional parameter v_14. */
    double mv14;

public:

    /**
     * Constructor.
     *
     * @param wntStimulus is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
     *
     */
    Lee2003WntSignallingOdeSystem(double wntStimulus=0.0);

    /**
     * Destructor.
     */
    ~Lee2003WntSignallingOdeSystem();

    /**
     * Initialise parameter values.
     */
    void Init();

    /**
     * Compute the RHS of the Lee et al. (2003) system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using Alarcons et al. (2004) system of equations).
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

};

#endif /*LEE2003WNTSIGNALLINGODESYSTEM_HPP_*/
