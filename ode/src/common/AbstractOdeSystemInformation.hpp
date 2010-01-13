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


#ifndef _ABSTRACTODESYSTEMINFORMATION_HPP_
#define _ABSTRACTODESYSTEMINFORMATION_HPP_

#include <vector>
#include <string>
#include <stdexcept>

/**
 * An abstract class which provides access to information about a particular
 * ODE system *class* (as opposed to an instance).
 *
 * The information available includes:
 *  - names of state variables
 *  - units of state variables
 *  - suggested initial conditions
 *
 * This class requires a subclass defining the Initialise method in order to set
 * up the information.  Developers may do this by defining their own subclass, but
 * the most convenient method is likely to be to use the OdeSystemInformation
 * class, which is a templated singleton subclass of this AbstractOdeSystemInformation
 * class.  See its documentation for details of how to use it.
 */
class AbstractOdeSystemInformation
{
protected:

    /** State variable names */
    std::vector<std::string> mVariableNames;

    /** State variable units */
    std::vector<std::string> mVariableUnits;

    /** Parameter names */
    std::vector<std::string> mParameterNames;

    /** Parameter units */
    std::vector<std::string> mParameterUnits;

    /** Suggested initial conditions */
    std::vector<double> mInitialConditions;

    /** Whether a 'real' Initialise method has been called */
    bool mInitialised;

    /**
     * Initialise the ODE system information.
     *
     * This must be provided by subclasses.
     */
    virtual void Initialise()=0;

public:

    /**
     * Constructor.
     */
    AbstractOdeSystemInformation();

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractOdeSystemInformation();

    /**
     * Set the suggested initial conditions to use.
     *
     * @param rInitialConditions  vector containing initial values for the state variables
     */
    void SetInitialConditions(const std::vector<double>& rInitialConditions);

    /**
     * Set a single component of the suggested initial conditions to use.
     *
     * @param index  the index of the state variable in the system
     * @param initialCondition  the initial value for the state variable
     */
    void SetInitialConditionsComponent(unsigned index, double initialCondition) throw(std::out_of_range);

    /**
     * Get a copy of the suggested initial conditions.
     */
    std::vector<double> GetInitialConditions() const;

    /**
     * Get the variable names vector.
     */
    const std::vector<std::string>& rGetVariableNames() const;

    /**
     * Get the variable units vector.
     */
    const std::vector<std::string>& rGetVariableUnits() const;

    /**
     * This method is used to establish a state varible's position within
     * the vector of state variables of an ODE system. This number can
     * then be used with the method GetStateVariableUnitsByNumber.
     *
     * @param rName  the name of a state variable
     * @return the state variable's position within the vector of state
     *         variables associated with the ODE system.
     */
    unsigned GetStateVariableNumberByName(const std::string& rName) const;

    /**
     * Get the units of a state variable given its index in the ODE system.
     *
     * @param varNumber  a state variable's position within the vector of
     *                   state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetStateVariableUnitsByNumber(unsigned varNumber) const;
    
    /**
     * Get the vector of parameter names.
     */
    const std::vector<std::string>& rGetParameterNames() const;

    /**
     * Get the vector of parameter units.
     */
    const std::vector<std::string>& rGetParameterUnits() const;
};


#endif /*_ABSTRACTODESYSTEMINFORMATION_HPP_*/
