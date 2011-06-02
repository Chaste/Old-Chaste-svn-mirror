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

#ifndef _ABSTRACTCVODESYSTEM_HPP_
#define _ABSTRACTCVODESYSTEM_HPP_

#include <vector>
#include <string>
#include <algorithm>


#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>
#include "ClassIsAbstract.hpp"

#include "AbstractParameterisedSystem.hpp"
#include "Exception.hpp"

// CVODE headers
#include <nvector/nvector_serial.h>

/**
 * Abstract OdeSystem class for Cvode systems (NVector instead of std::vector)
 *
 * Sets up variables and functions for a general CVODE system.
 *
 * ODE systems are specified primarily by the EvaluateYDerivatives() method,
 * which calculates the right-hand side of the system.
 *
 * Instances can store their state internally in the mStateVariables vector
 * in our base class AbstractParameterisedSystem (see also
 * GetNumberOfStateVariables(), SetStateVariables() and rGetStateVariables()),
 * although this is not essential - the vector may be empty, in which case
 * AbstractIvpOdeSolver::SolveAndUpdateStateVariable may not be used to
 * solve the system.
 *
 * CVODE systems may also have a vector of parameters, which can be accessed
 * through the GetParameter() and SetParameter() methods of our base class.
 *
 * Information about what the parameters and state variables represent is
 * provided by a subclass of AbstractOdeSystemInformation.  Various wrapper
 * methods (e.g. rGetStateVariableNames()) are provided in our base class to
 * access this information.
 *
 * Also, subclasses may define a condition at which ODE solvers should stop
 * prematurely. For this class CVODE solvers are being used, so
 * CalculateRootFunction() should be used to detect the stopping time.
 */
class AbstractCvodeSystem : public AbstractParameterisedSystem<N_Vector>
{

protected:

    /** Whether to use an analytic Jacobian. */
    bool mUseAnalyticJacobian;

public:

    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system
     */
    AbstractCvodeSystem(unsigned numberOfStateVariables);

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractCvodeSystem();

    /**
     * Method to evaluate the derivatives of the system.
     *
     * @param time  the current time
     * @param y  the current values of the state variables
     * @param ydot  storage for the derivatives of the system; will be filled in on return
     */
    virtual void EvaluateYDerivatives(realtype time,
                                      N_Vector y,
                                      N_Vector ydot)=0;

//    /**
//     * An alternative approach to stopping events; currently only useful with CVODE.
//     * CVODE can search for roots (zeros) of this function while solving the ODE system,
//     * and home in on them to find sign transitions to high precision.
//     *
//     * The default implementation here fakes a root function using CalculateStoppingEvent.
//     *
//     * @param time  the current time
//     * @param rY  the current values of the state variables
//     */
//    virtual double CalculateRootFunction(double time, const std::vector<double>& rY);
//
//    /**
//     * Get whether an analytic Jacobian is used.
//     *
//     * @return mUseAnalyticJacobian
//     */
//    bool GetUseAnalyticJacobian();

//    /**
//     * \todo move to AbstractParameterisedSystem? (1540)
//     *
//     * @return const reference to the state variables in the ODE system (used in archiving).
//     */
//    const std::vector<double>& rGetConstStateVariables() const;

};

//CLASS_IS_ABSTRACT(AbstractCvodeSystem)

#endif //_ABSTRACTCVODESYSTEM_HPP_
