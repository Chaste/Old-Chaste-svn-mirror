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

#ifndef ABSTRACTPARAMETERISEDSYSTEM_HPP_
#define ABSTRACTPARAMETERISEDSYSTEM_HPP_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "AbstractOdeSystemInformation.hpp"

/**
 * This class contains the state variable and parameter vectors for an ODE system,
 * along with methods to access these.  It also holds the AbstractOdeSystemInformation
 * pointer, and methods to access this object to provide information about the ODE
 * system, such as state variable/parameter names and units.
 * 
 * Its main purpose is to be a common base class for both AbstractOdeSystem and
 * AbstractCvodeCell, which require similar functionality but use different vector
 * types.
 */
template<typename VECTOR>
class AbstractParameterisedSystem
{
protected:
    /** The number of state variables in the system. */
    unsigned mNumberOfStateVariables;

    /** Vector containing the current values of the state variables. */
    VECTOR mStateVariables;

    /** Vector containing parameter values. */
    VECTOR mParameters;

    /**
     * Information about the concrete ODE system class.
     *
     * Subclasses @b need to set this in their constructor to point to an instance
     * of a suitable class.  See for example the OdeSystemInformation class.
     */
    boost::shared_ptr<AbstractOdeSystemInformation> mpSystemInfo;

public:
    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system
     */
    AbstractParameterisedSystem(unsigned numberOfStateVariables);
    
    /**
     * Virtual destructor.
     */
    virtual ~AbstractParameterisedSystem();

    /**
     * Get the object which provides information about this ODE system.
     */
    boost::shared_ptr<const AbstractOdeSystemInformation> GetSystemInformation() const;
    
    /**
     * Get the name of this system.
     */
    std::string GetSystemName() const;

    //
    // State variable methods
    //

    /**
     * Get the number of state variables in the ODE system.
     *
     * @return mNumberOfStateVariables
     */
    unsigned GetNumberOfStateVariables() const;

    /**
     * Get the value of a given state variable.
     *
     * @param index the index of the state variable
     */
    double GetStateVariable(unsigned index) const;

    /**
     * Set the value of a single state variable in the ODE system.
     *
     * @param index index of the state variable to be set
     * @param newValue new value of the state variable
     */
    void SetStateVariable(unsigned index, double newValue);

    /**
     * Get the names of the state variables in the ODE system.
     */
    const std::vector<std::string>& rGetStateVariableNames() const;

    /**
     * Get the units of the state variables in the ODE system.
     */
    const std::vector<std::string>& rGetStateVariableUnits() const;

    /**
     * This method is used to establish a state variable's position within
     * the vector of state variables of an ODE system.  This number can
     * then be used with the methods GetStateVariable and GetStateVariableUnits.
     *
     * @param rName  the name of a state variable.
     *
     * @return the state variable's position within the vector of state variables
     *         associated with the ODE system.
     */
    unsigned GetStateVariableIndex(const std::string& rName) const;

    /**
     * Get the units of a state variable given its index in the ODE system.
     *
     * @param index  a state variable's position within the vector of
     *               state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetStateVariableUnits(unsigned index) const;

    //
    // Parameter methods
    //

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
     * This method is used to establish a parameter's position within
     * the vector of parameters of an ODE system. This number can
     * then be used with the methods GetParameterUnits and GetParameter.
     *
     * @param rName  the name of a parameter
     * @return the parameter's position within the vector of parameters
     *         associated with the ODE system.
     */
    unsigned GetParameterIndex(const std::string& rName) const;

    /**
     * Get the units of a parameter given its index in the ODE system.
     *
     * @param index  a state variable's position within the vector of
     *               state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetParameterUnits(unsigned index) const;

    //
    // "Any variable" methods
    //

    /**
     * Get the value of a variable, whether a state variable, parameter,
     * or derived quantity.
     * 
     * Note that if the variable is a derived quantity, this method will compute
     * all derived quantities, so may not be very efficient.
     *
     * @param index the index of the variable, as given by GetAnyVariableIndex.
     * @param time  the current simulation time, possibly needed if the variable
     *     is a derived quantity
     */
    double GetAnyVariable(unsigned index, double time=0.0);

    /**
     * Get the index of a variable, whether a state variable, parameter,
     * or derived quantity, with the given name.
     * The returned index is suitable for use with GetAnyVariableUnits
     * and GetAnyVariable.
     *
     * @param rName  the name of a variable
     */
    unsigned GetAnyVariableIndex(const std::string& rName) const;

    /**
     * Get the units of a variable, whether a state variable, parameter, or
     * derived quantity, given its index as returned by GetAnyVariableIndex.
     *
     * @param index  an index from GetAnyVariableIndex.
     * @return the units of the variable.
     */
    std::string GetAnyVariableUnits(unsigned index) const;

    //
    // Derived quantity methods
    //
    
    /**
     * Get the number of derived quantities.
     */
    unsigned GetNumberOfDerivedQuantities() const;
    
    /**
     * Compute the derived quantities from the given system state.
     * Uses the current values for the parameters.
     * 
     * @param time  the time at which to compute the derived quantities
     * @param rState  values for the state variables
     */
    virtual VECTOR ComputeDerivedQuantities(double time,
                                            const VECTOR& rState);
    
    /**
     * Compute the derived quantities based on the current system state.
     * 
     * @param time  the time at which to compute the derived quantities
     */
    VECTOR ComputeDerivedQuantitiesFromCurrentState(double time);

    /**
     * Get the vector of derived quantity names.
     */
    const std::vector<std::string>& rGetDerivedQuantityNames() const;

    /**
     * Get the vector of derived quantity units.
     */
    const std::vector<std::string>& rGetDerivedQuantityUnits() const;
    
    /**
     * Get the index of a derived quantity, given its name.
     *
     * @param rName  the name of a derived quantity.
     */
    unsigned GetDerivedQuantityIndex(const std::string& rName) const;

    /**
     * Get the units of a derived quantity.
     *
     * @param index  an index from GetDerivedQuantityIndex.
     * @return the units of the variable.
     */
    std::string GetDerivedQuantityUnits(unsigned index) const;
};

#endif /*ABSTRACTPARAMETERISEDSYSTEM_HPP_*/
