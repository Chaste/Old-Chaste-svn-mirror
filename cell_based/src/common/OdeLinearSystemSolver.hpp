/*

Copyright (C) University of Oxford, 2005-2011

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
#ifndef ODELINEARSYSTEMSOLVER_HPP_
#define ODELINEARSYSTEMSOLVER_HPP_

#include "LinearSystem.hpp"

// todo doxygen
class OdeLinearSystemSolver
{
private:
        
    /** Timestep for solver */
    double mTimeStep;
    
    /** The LHS matrix and the force vector. Will call solve on this. */
    LinearSystem mLinearSystem;
    
    /** Initial conditions vector */
    Vec mInitialConditionsVector;
    
    /** Force vector (F in Mdr/dt = F) */
    Vec mForceVector; 
    
public:

    /** Constructor
     *  @param systemSize size of the ODE system
     *  @param timeStep timestep
     */
    OdeLinearSystemSolver(unsigned systemSize, double timeStep);
    
    /** Destructor */
    ~OdeLinearSystemSolver();
       
    /** Get the timestep for the solver */       
    double GetTimeStep();
    
    /** Get the matrix */
    Mat& rGetLhsMatrix();
    
    /** Get the force vector */
    Vec& rGetForceVector();
    
    /** Get the vector of initial conditions */
    Vec& rGetInitialConditionVector();
    
    /** Solve method */
    Vec SolveOneTimeStep();
    
};

#endif /*ODELINEARSYSTEMSOLVER_HPP_*/


