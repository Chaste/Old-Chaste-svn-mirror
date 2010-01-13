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

#ifndef _RUNGEKUTTA4IVPODESOLVER_HPP_
#define _RUNGEKUTTA4IVPODESOLVER_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractOneStepIvpOdeSolver.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * A concrete one step ODE solver class that employs the Runge Kutta
 * 4th order solver (RK4).
 */
class RungeKutta4IvpOdeSolver : public AbstractOneStepIvpOdeSolver
{
private:
    
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the abstract IVP Solver, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractOneStepIvpOdeSolver>(*this);
    }
    
protected:

    /**
     * Calculate the solution to the ODE system at the next timestep.
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rNextYValues  the state at the next timestep
     */
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& rCurrentYValues,
                             std::vector<double>& rNextYValues);

private:

    std::vector<double> k1;  /**< Working memory: expression k1 in the RK4 method. */
    std::vector<double> k2;  /**< Working memory: expression k2 in the RK4 method. */
    std::vector<double> k3;  /**< Working memory: expression k3 in the RK4 method. */
    std::vector<double> k4;  /**< Working memory: expression k4 in the RK4 method. */
    std::vector<double> yki; /**< Working memory: expression yki in the RK4 method. */

};

#include "TemplatedExport.hpp"
CHASTE_CLASS_EXPORT(RungeKutta4IvpOdeSolver);

#endif //_RUNGEKUTTA4IVPODESOLVER_HPP_
