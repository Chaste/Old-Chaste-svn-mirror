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
#ifndef _TYSONNOVAK2001ODESYSTEM_HPP_
#define _TYSONNOVAK2001ODESYSTEM_HPP_

#include <cmath>

#include "AbstractOdeSystemWithAnalyticJacobian.hpp"

/**
 * Represents the Tyson & Novak 2001 system of ODEs.
 */
class TysonNovak2001OdeSystem : public AbstractOdeSystemWithAnalyticJacobian
{
private:

    // Parameters for the Tyson Novak 2001 model
    double mK1;
    double mK2d;
    double mK2dd;
    double mK2ddd;
    double mCycB_threshold;
    double mK3d;
    double mK3dd;
    double mK4d;
    double mK4;
    double mJ3;
    double mJ4;
    double mK5d;
    double mK5dd;
    double mK6;
    double mJ5;
    double mN;
    double mK7;
    double mK8;
    double mJ7;
    double mJ8;
    double mMad;
    double mK9;
    double mK10;
    double mK11;
    double mK12d;
    double mK12dd;
    double mK12ddd;
    double mKeq;
    double mK13;
    double mK14;
    double mK15d;
    double mK15dd;
    double mK16d;
    double mK16dd;
    double mJ15;
    double mJ16;
    double mMu;
    double mMstar;
    
public:

    /**
     * Constructor.
     */
    TysonNovak2001OdeSystem();
    
    /**
     *  Destructor.
     */
    ~TysonNovak2001OdeSystem();
    
    /**
     * Initialise parameter values.
     */
    void Init();
    
    /**
     * Compute the RHS of the Alarcon et al. (2004) system of ODEs.
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
     * Calculate whether the conditions for the cell cycle to finish have been met.
     * 
     * @param time at which to calculate whether the stopping event has occured
     * @param rY value of the solution vector used to evaluate the RHS.
     */
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY);
    
    /**
     * Compute the Jacobian of the ODE system.
     * 
     * @param solutionGuess initial guess for the solution vector.
     * @param jacobian the Jacobian of the ODE system.
     * @param time at which to calculate the Jacobian.
     * @param timeStep used to calculate the Jacobian.
     */
    virtual void AnalyticJacobian(const std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep);
};

#endif //_TYSONNOVAK2001ODESYSTEM_HPP_
