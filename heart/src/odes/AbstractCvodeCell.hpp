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

#ifdef CHASTE_CVODE
#ifndef _ABSTRACTCVODECELL_HPP_
#define _ABSTRACTCVODECELL_HPP_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

// Chaste headers
#include "OdeSolution.hpp"
#include "AbstractOdeSystemInformation.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractIvpOdeSolver.hpp"

// CVODE headers
#include <nvector/nvector_serial.h>

/**
 * A cardiac cell that is designed to be simulated using CVODE.
 * It uses CVODE's vector type natively.
 *
 * Functionality is similar to that provided by AbstractCardiacCell and AbstractOdeSystem,
 * but not identical.  It also includes a direct interface to the CVODE solver, via the
 * Solve methods, since the CvodeAdaptor class doesn't work for us.
 *
 * Assumes that it will be solving stiff systems, so uses BDF/Newton.
 *
 * Note that a call to Solve will initialise the CVODE solver, and free its
 * working memory when done.  There is thus a non-trivial overhead involved.
 *
 * \todo Add an option to just initialise once, and assume subsequent Solve
 *   calls are continuing from where we left off.
 *
 * \todo Integrate this better into the main hierarchy?
 */
class AbstractCvodeCell
{
protected:
    /** The number of state variables. */
    unsigned mNumberOfStateVariables;
    /** The index of the transmembrane potential within the state variable vector */
    unsigned mVoltageIndex;
    /** The state variables. */
    N_Vector mStateVariables;

    /**
     * Information about the concrete ODE system class.
     *
     * Subclasses need to set this in their constructor to point to an instance
     * of a suitable class.  See for example the OdeSystemInformation class.
     */
    boost::shared_ptr<AbstractOdeSystemInformation> mpSystemInfo;

    /** The intracellular stimulus current. */
    boost::shared_ptr<AbstractStimulusFunction> mpIntracellularStimulus;

    /** Relative tolerance for solver. */
    double mRelTol;
    /** Absolute tolerance for solver. */
    double mAbsTol;

    /** CVODE's internal data. */
    void* mpCvodeMem;
    /**
     * The maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int mMaxSteps;

    /** The size of the previous timestep. */
    double mLastInternalStepSize;

    /**
     * Whether to clamp the voltage by setting dV/dt to zero.
     */
    bool mSetVoltageDerivativeToZero;

    /**
     * Used to include extra debugging information in exception messages, e.g.
     *  EXCEPTION(DumpState("Gating variable out of range"));
     */
    std::string DumpState(const std::string& message,
                          N_Vector Y = NULL);

    /**
     * Can be called by concrete subclass constructors to initialise the state
     * variables.
     *
     * Mainly here to reduce differences between this class and AbstractCardiacCell.
     */
    void Init();

    /**
     * This class is used by GetStateVariables and SetStateVariablesUsingACopyOfThisVector()
     */
    N_Vector CopyVector(N_Vector originalVec);

public:
    /**
     * Create a new cardiac cell.
     *
     * @param numberOfStateVariables  the size of the ODE system modelling this cell
     * @param voltageIndex  the index of the transmembrane potential within the vector of state variables
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    AbstractCvodeCell(boost::shared_ptr<AbstractIvpOdeSolver> /* unused; should be empty */,
                      unsigned numberOfStateVariables,
                      unsigned voltageIndex,
                      boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Free the state variables, if they have been set.
     */
    virtual ~AbstractCvodeCell();

    /** Get the index of V within the state variables. */
    unsigned GetVoltageIndex();

    /**
     * Get the current value of the transmembrane potential, as given
     * in our state variable vector.
     */
    double GetVoltage();

    /**
     * Set the transmembrane potential
     * @param voltage  new value
     */
    void SetVoltage(double voltage);

    /**
     * Set the intracellular stimulus.
     * Shorthand for SetIntracellularStimulusFunction.
     * @param pStimulus  new stimulus function
     */
    void SetStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus);

    /**
     * Get the value of the intracellular stimulus.
     * Shorthand for GetIntracellularStimulus.
     * @param time  the time at which to evaluate the stimulus
     */
    double GetStimulus(double time);

    /** Get the number of state variables */
    unsigned GetNumberOfStateVariables();

    /** Get the names of this cell's state variables. */
    const std::vector<std::string>& rGetStateVariableNames() const;

    /** Get the units of this cell's state variables. */
    const std::vector<std::string>& rGetStateVariableUnits() const;

    /**
     * This method is used to establish a state variable's position within
     * the vector of state variables of an ODE system. This number can
     * then be used with the methods GetStateVariable and
     * GetStateVariableUnits.
     *
     * @param name  the name of a state variable.
     *
     * @return the state variable's position within the vector of state variables
     *         associated with the ODE system.
     */
    unsigned GetStateVariableIndex(const std::string name) const;

    /**
     * Get the value of a state variable given its index in the ODE system.
     *
     * @param varNumber  a state variable's position within the vector of
     *                   state variables associated with the ODE system.
     *
     * @return the current value of the state variable.
     */
    double GetStateVariable(unsigned varNumber) const;

    /**
     * Get the units of a state variable given its index in the ODE system.
     *
     * @param varNumber  a state variable's position within the vector of
     *                   state variables associated with the ODE system.
     * @return the units of the state variable.
     */
    std::string GetStateVariableUnits(unsigned varNumber) const;

    /**
     * Get the object which provides information about this ODE system.
     */
    boost::shared_ptr<const AbstractOdeSystemInformation> GetSystemInformation() const;

    /**
     * Get the initial conditions for the cell.
     *
     * Creates and returns a fresh N_Vector, which must be destroyed
     * by the caller when finished with.
     */
    N_Vector GetInitialConditions();

    /**
     * Assign a vector to be used for this cell's state.
     *
     * The cell takes responsibility for freeing this vector when it
     * is destroyed.  If the cell already has state, it will be freed.
     *
     * @param stateVars  new state variables vector
     */
    void SetStateVariables(N_Vector stateVars);

    /**
     * Assign a vector to be copied for this cell's state.
     *
     * Caller retains responsibility for freeing the vector.
     *
     * @param stateVars  new state variables vector
     */
    void SetStateVariablesUsingACopyOfThisVector(N_Vector stateVars);

    /**
     * Takes a copy of the state variable vector.
     * Doesn't really return a vector (N_Vector is a pointer type)
     * but named like this to match AbstractCardiacCell.
     * Caller takes responsibility for freeing the vector.
     */
    N_Vector GetStateVariables();

    /**
     * Get the state variable vector.
     * Doesn't really return a reference (N_Vector is a pointer type)
     * but named like this to match AbstractCardiacCell.
     * This cell will retain responsibility for freeing the vector.
     */
    N_Vector rGetStateVariables();

    /**
     * RHS evaluation function, to be provided by subclasses.
     *
     * @param t  the time at which to evaluate the RHS
     * @param y  the values of the state variables at time t
     * @param ydot  to be filled in with the derivatives of the state variables
     */
    virtual void EvaluateRhs(realtype t,
                             N_Vector y,
                             N_Vector ydot)=0;

    /**
     * Empty method which can be over-ridden in concrete cell class which should
     * go through the current state vector and go range checking on the values
     * (eg check that concentrations are positive and gating variables are between
     * zero and one). This method is called in the Solve methods, at the end of each
     * sampling timestep, and the end of a simulation.
     */
    virtual void VerifyStateVariables()
    {
    }

    //
    // Solver methods
    //

    /**
     * Set whether to clamp the voltage by setting its derivative to zero.
     * @param clamp whether to clamp
     */
    void SetVoltageDerivativeToZero(bool clamp=true);

    /**
     * Simulate the cell, returning a sampling of the state variables.
     *
     * Uses the current values of the state variables at initial conditions.
     * If the state variables have not been set (either by a prior solve, or
     * a call to SetStateVariables) the initial conditions (given by
     * GetInitialConditions) will be used.
     *
     * The final values of the state variables will also be stored in this object.
     *
     * @param tStart  start time of simulation
     * @param tEnd  end time of simulation
     * @param maxDt  maximum time step to be taken by the adaptive solver
     *   (set this appropriately to avoid missing a stimulus)
     * @param tSamp  sampling interval at which to store results
     */
    OdeSolution Solve(realtype tStart,
                      realtype tEnd,
                      realtype maxDt,
                      realtype tSamp);

    /**
     * Simulate the cell, updating its internal state variables.
     *
     * Uses the current values of the state variables at initial conditions.
     * If the state variables have not been set (either by a prior solve, or
     * a call to SetStateVariables) the initial conditions (given by
     * GetInitialConditions) will be used.
     *
     * @param tStart  start time of simulation
     * @param tEnd  end time of simulation
     * @param maxDt  maximum time step to be taken by the adaptive solver
     *   (set this appropriately to avoid missing a stimulus)
     */
    void Solve(realtype tStart,
               realtype tEnd,
               realtype maxDt);

    /**
     * Change the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.  Default is 500 (set by CVODE).
     */
    void SetMaxSteps(long int numSteps);

    /**
     * Get the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int GetMaxSteps();

    /**
     * Set relative and absolute tolerances; both scalars.
     * If no parameters are given, tolerances will be reset to default values.
     *
     * @param relTol  the relative tolerance for the solver (defaults to 1e-4)
     * @param absTol  the absolute tolerance for the solver (defaults to 1e-6)
     */
    void SetTolerances(double relTol=1e-4, double absTol=1e-6);

    /**
     * Get the relative tolerance.
     */
    double GetRelativeTolerance();

    /**
     * Get the absolute tolerance.
     */
    double GetAbsoluteTolerance();

    /**
     * Get the last step size used internally by CVODE in the last Solve call.
     */
    double GetLastStepSize();

private:

    /**
     * Set up the CVODE data structures needed to solve the given system.
     *
     * @param initialConditions  initial conditions
     * @param tStart  start time of simulation
     * @param maxDt  maximum time step to take
     */
    void SetupCvode(N_Vector initialConditions,
                    realtype tStart,
                    realtype maxDt);

    /** Free CVODE memory after a solve. */
    void FreeCvodeMemory();

    /**
     * Report an error from CVODE.
     *
     * @param flag  CVODE error code
     * @param msg  Our description of the error
     */
    void CvodeError(int flag, const char * msg);

    /**
     * Copy an N_Vector into a std::vector<double>
     *
     * @param v  vector to copy
     */
    std::vector<double> MakeStdVec(N_Vector v);

};


#endif // _ABSTRACTCVODECELL_HPP_
#endif // CHASTE_CVODE
