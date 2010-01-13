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

#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_


#include <cassert>
#include <climits> // Work around a boost bug - see #1024.

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include "ClassIsAbstract.hpp"
#include <boost/shared_ptr.hpp>

#include "Exception.hpp"
#include "AbstractOdeSystemInformation.hpp"

/**
 * Abstract OdeSystem class.
 *
 * Sets up variables and functions for a general ODE system.
 *
 * ODE systems are specified primarily by the EvaluateYDerivatives method,
 * which calculates the right-hand side of the system.
 *
 * Instances can store their state internally in the #mStateVariables vector
 * (see GetNumberOfStateVariables, SetStateVariables and rGetStateVariables),
 * although this is not essential - the vector may be empty, although in this
 * case AbstractIvpOdeSolver::SolveAndUpdateStateVariable may not be used to
 * solve the system.
 *
 * ODE systems may also have a vector of parameters, which can be accessed
 * through the GetParameter and SetParameter methods.
 *
 * Information about what the parameters and state variables represent is
 * provided by a subclass of AbstractOdeSystemInformation.  Various wrapper
 * methods (e.g. rGetVariableNames) are provided in this class to access
 * this information.
 *
 * There are two more advanced facilities available for subclass authors.
 * An analytic form for the Jacobian matrix of the system may be provided,
 * in which case you must subclass AbstractOdeSystemWithAnalyticJacobian.
 * The GetUseAnalyticJacobian method will test whether this is the case.
 *
 * Also, subclasses may define a condition at which ODE solvers should stop
 * prematurely.  For the Chaste solvers this is done by overriding
 * CalculateStoppingEvent; if the more advanced CVODE solvers are being used
 * then implement CalculateRootFunction instead to detect the stopping time
 * more accurately.
 */
class AbstractOdeSystem
{
    friend class TestAbstractOdeSystem;

    /** Needed for serialization. */
    friend class boost::serialization::access;
private:
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mNumberOfStateVariables;
        archive & mUseAnalyticJacobian;
        archive & mStateVariables;
        archive & mParameters;

//        archive &mpSystemInfo;
        /// \todo archive mpSystemInfo
    }

protected:

    /** The number of state variables in the ODE system. */
    unsigned mNumberOfStateVariables;

    /** Vector containing the current values of the state variables. */
    std::vector<double> mStateVariables;

    /** Vector containing parameters. */
    std::vector<double> mParameters;

    /**
     * Information about the concrete ODE system class.
     *
     * Subclasses need to set this in their constructor to point to an instance
     * of a suitable class.  See for example the OdeSystemInformation class.
     */
    boost::shared_ptr<AbstractOdeSystemInformation> mpSystemInfo;

    /** Whether to use an analytic Jacobian. */
    bool mUseAnalyticJacobian;

public:

    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system (defaults to 0)
     */
    AbstractOdeSystem(unsigned numberOfStateVariables = 0);

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractOdeSystem();

    /**
     * Method to evaluate the derivatives of the system.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @param rDY  storage for the derivatives of the system; will be filled in on return
     */
    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                                      std::vector<double>& rDY)=0;

    /**
     * Get the number of state variables in the ODE system.
     *
     * @return mNumberOfStateVariables
     */
    unsigned GetNumberOfStateVariables() const;


    /**
     * Get the number of parameters.
     */
    unsigned GetNumberOfParameters() const;

    /**
     * Get the value of a given parameter.
     *
     * @param index the index of the parameter
     */
    double GetParameter(unsigned index) const;

    /**
     * Set the value of a given parameter.
     *
     * @param index the index of the parameter
     * @param value the value
     */
    void SetParameter(unsigned index, double value);

    /**
     * Get the names of the parameters in the ODE system.
     */
    const std::vector<std::string>& rGetParameterNames() const;

    /**
     * Get the units of the parameters in the ODE system.
     */
    const std::vector<std::string>& rGetParameterUnits() const;


    /**
     * Set the initial conditions for the ODE system.
     *
     * @param rInitialConditions  vector containing initial values for the state variables
     */
    void SetInitialConditions(const std::vector<double>& rInitialConditions);

    /**
     * Set the initial condition one state variable.
     *
     * @param index  the index of the state variable in the system
     * @param initialCondition  the initial value for the state variable
     */
    void SetInitialConditionsComponent(unsigned index, double initialCondition);

    /**
     * Get the initial conditions for the ODE system.
     */
    std::vector<double> GetInitialConditions() const;

    /**
     * Set the values of the state variables in the ODE system.
     *
     * @param rStateVariables vector containing values for the state variables
     */
    void SetStateVariables(const std::vector<double>& rStateVariables);

    /**
	 * Set the value of a single state variable in the ODE system.
	 *
	 * @param stateVariable index of the state variable to be set
	 * @param newValue new value of the state variable
	 */
	void SetStateVariable(unsigned stateVariable, double newValue);

    /**
     * Get the values of the state variables in the ODE system.
     */
    std::vector<double>& rGetStateVariables();

    /**
     * Get the names of the state variables in the ODE system.
     */
    const std::vector<std::string>& rGetVariableNames() const;

    /**
     * Get the units of the state variables in the ODE system.
     */
    const std::vector<std::string>& rGetVariableUnits() const;

    /**
     *  CalculateStoppingEvent() - can be overloaded if the ODE is to be solved
     *  only until a particular event (for example, only until the y value becomes
     *  negative.
     *
     *  After each timestep the solver will call this method on the ODE to see if
     *  it should stop there. By default, false is returned here.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     */
    virtual bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    /**
     * An alternative approach to stopping events; currently only useful with CVODE.
     * CVODE can search for roots (zeros) of this function while solving the ODE system,
     * and home in on them to find sign transitions to high precision.
     *
     * The default implementation here fakes a root function using CalculateStoppingEvent.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     */
    virtual double CalculateRootFunction(double time, const std::vector<double>& rY);

    /**
     * Get whether an analytic Jacobian is used.
     *
     * @return mUseAnalyticJacobian
     */
    bool GetUseAnalyticJacobian();

    /**
     * This method is used to establish a state variable's position within
     * the vector of state variables of an ODE system. This number can
     * then be used with the methods GetStateVariableValueByNumber and
     * GetStateVariableUnitsByNumber.
     *
     * @param name  the name of a state variable.
     *
     * @return the state variable's position within the vector of state variables
     *         associated with the ODE system.
     */
    unsigned GetStateVariableNumberByName(const std::string name);

    /**
     * Get the value of a state variable given its index in the ODE system.
     *
     * @param varNumber  a state variable's position within the vector of
     *                   state variables associated with the ODE system.
     *
     * @return the current value of the state variable.
     */
    double GetStateVariableValueByNumber(unsigned varNumber) const;

    /**
     * Get the units of a state variable given its index in the ODE system.
     *
     * @param varNumber  a state variable's position within the vector of
     *                   state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetStateVariableUnitsByNumber(unsigned varNumber) const;

    /**
     * Get the object which provides information about this ODE system.
     */
    boost::shared_ptr<const AbstractOdeSystemInformation> GetSystemInformation() const;

protected:

    /**
     * Used to include extra debugging information in exception messages.
     * For example,
     *      EXCEPTION(DumpState("Gating variable out of range"));
     *
     * @param rMessage  the exception message
     * @param Y  the values of the state variables (optional input argument)
     */
    std::string DumpState(const std::string& rMessage,
                          std::vector<double> Y = std::vector<double>());
};

CLASS_IS_ABSTRACT(AbstractOdeSystem)


#endif //_ABSTRACTODESYSTEM_HPP_
