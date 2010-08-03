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


#ifndef ABSTRACTBACKWARDEULERCARDIACCELL_HPP_
#define ABSTRACTBACKWARDEULERCARDIACCELL_HPP_

#include <cassert>
#include <cmath>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "ClassIsAbstract.hpp"

#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"

/**
 * This is the base class for cardiac cells solved using a (decoupled) backward
 * Euler approach.
 *
 * The basic approach to solving such models is:
 *  \li Update the transmembrane potential, either from solving an external PDE,
 *      or using a forward Euler step.
 *  \li Update any gating variables (or similar) using a backward euler step.
 *      Suitable ODEs can be written in the form  \f$du/dt = g(V) + h(V)*u\f$.  The update
 *      expression is then \f$u_n = ( u_{n-1} + g(V_n)*dt ) / ( 1 - h(V_n)*dt )\f$.
 *  \li Update the remaining state variables using Newton's method to solve the
 *      nonlinear system \f$U_n - U_{n-1} = dt*F(U_n, V_n)\f$.
 *      The template parameter to the class specifies the size of this nonlinear system.
 */
template<unsigned SIZE>
class AbstractBackwardEulerCardiacCell : public AbstractCardiacCell
{
    private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
    }
public:

    /**
     * Standard constructor for a cell.
     *
     * @param numberOfStateVariables  the size of the ODE system
     * @param voltageIndex  the index of the variable representing the transmembrane
     *     potential within the state variable vector
     * @param pIntracellularStimulus  the intracellular stimulus function
     *
     * Some notes for future reference:
     *  \li We may want to remove the timestep from this class, and instead pass it to
     *      the Compute* methods, especially if variable timestepping is to be used.
     *  \li It's a pity that inheriting from AbstractCardiacCell forces us to store a
     *      null pointer (for the unused ODE solver) in every instance.  We may want
     *      to revisit this design decision at a later date.
     */
    AbstractBackwardEulerCardiacCell(
        unsigned numberOfStateVariables,
        unsigned voltageIndex,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /** Virtual destructor */
    virtual ~AbstractBackwardEulerCardiacCell();

    /**
     * Compute the residual of the nonlinear system portion of the cell model.
     *
     * @param time  the current time
     * @param rCurrentGuess  the current guess for \f$U_n\f$
     * @param rResidual  to be filled in with the residual vector
     */
    virtual void ComputeResidual(double time, const double rCurrentGuess[SIZE], double rResidual[SIZE])=0;

    /**
     * Compute the Jacobian matrix for the nonlinear system portion of the cell model.
     *
     * @param time  the current time
     * @param rCurrentGuess  the current guess for \f$U_n\f$
     * @param rJacobian  to be filled in with the Jacobian matrix
     */
    virtual void ComputeJacobian(double time, const double rCurrentGuess[SIZE], double rJacobian[SIZE][SIZE])=0;

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.  Uses a forward Euler step to update the transmembrane
     * potential at each timestep.
     *
     * The length of the time interval must be a multiple of the timestep.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     * @param tSamp  sampling interval for returned results (defaults to #mDt)
     * @return  the values of each state variable, at intervals of tSamp.
     */
    virtual OdeSolution Compute(double tStart, double tEnd, double tSamp=0.0);

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.  The transmembrane potential is kept fixed throughout.
     *
     * The length of the time interval must be a multiple of the timestep.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void ComputeExceptVoltage(double tStart, double tEnd);

private:
#define COVERAGE_IGNORE
    /**
     * This function should never be called - the cell class incorporates its own solver.
     *
     * @param time
     * @param rY
     * @param rDY
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
        NEVER_REACHED;
    }
#undef COVERAGE_IGNORE

protected:
    /**
     * Compute the values of all state variables, except the voltage, using backward Euler,
     * for one timestep from tStart.
     *
     * \note This method must be provided by subclasses.
     *
     * @param tStart  start of this timestep
     */
    virtual void ComputeOneStepExceptVoltage(double tStart)=0;

    /**
     * Perform a forward Euler step to update the transmembrane potential.
     *
     * \note This method must be provided by subclasses.
     *
     * @param time  start of this timestep
     */
    virtual void UpdateTransmembranePotential(double time)=0;
};


/*
 * NOTE: Explicit instantiation is not used for this class, because the SIZE
 * template parameter could take arbitrary values.
 */


template <unsigned SIZE>
AbstractBackwardEulerCardiacCell<SIZE>::AbstractBackwardEulerCardiacCell(
    unsigned numberOfStateVariables,
    unsigned voltageIndex,
    boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver>(),
                              numberOfStateVariables,
                              voltageIndex,
                              pIntracellularStimulus)
{}

template <unsigned SIZE>
AbstractBackwardEulerCardiacCell<SIZE>::~AbstractBackwardEulerCardiacCell()
{}

template <unsigned SIZE>
OdeSolution AbstractBackwardEulerCardiacCell<SIZE>::Compute(double tStart, double tEnd, double tSamp)
{
    // In this method, we iterate over timesteps, doing the following for each:
    //   - update V using a forward Euler step
    //   - call ComputeExceptVoltage(t) to update the remaining state variables
    //     using backward Euler
    
    // Check length of time interval
    if (tSamp < mDt)
    {
        tSamp = mDt;
    }
    double _n_steps = (tEnd - tStart) / tSamp;
    const unsigned n_steps = (unsigned) floor(_n_steps+0.5);
    assert(fabs(tStart+n_steps*tSamp - tEnd) < 1e-12);
    const unsigned n_small_steps = (unsigned) floor(tSamp/mDt+0.5);
    assert(fabs(mDt*n_small_steps - tSamp) < 1e-12);

    // Initialise solution store
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(n_steps);
    solutions.rGetSolutions().push_back(rGetStateVariables());
    solutions.rGetTimes().push_back(tStart);
    solutions.SetOdeSystemInformation(this->mpSystemInfo);

    // Loop over time
    double curr_time = tStart;
    for (unsigned i=0; i<n_steps; i++)
    {
        for (unsigned j=0; j<n_small_steps; j++)
        {
            curr_time = tStart + i*tSamp + j*mDt;

            // Compute next value of V
            UpdateTransmembranePotential(curr_time);
    
            // Compute other state variables
            ComputeOneStepExceptVoltage(curr_time);

            // check gating variables are still in range
            VerifyStateVariables();
        }

        // Update solutions
        solutions.rGetSolutions().push_back(rGetStateVariables());
        solutions.rGetTimes().push_back(curr_time+mDt);
    }

    return solutions;
}

template <unsigned SIZE>
void AbstractBackwardEulerCardiacCell<SIZE>::ComputeExceptVoltage(double tStart, double tEnd)
{
    // This method iterates over timesteps, calling ComputeExceptVoltage(t) at
    // each one, to update all state variables except for V, using backward Euler.
    // Check length of time interval
    double _n_steps = (tEnd - tStart) / mDt;
    unsigned n_steps = (unsigned) floor(_n_steps+0.5);
    assert(fabs(tStart+n_steps*mDt - tEnd) < 1e-12);

    // Loop over time
    double curr_time;
    for (unsigned i=0; i<n_steps; i++)
    {
        curr_time = tStart + i*mDt;

        // Compute other state variables
        ComputeOneStepExceptVoltage(curr_time);

        // check gating variables are still in range
        VerifyStateVariables();
    }
}

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractBackwardEulerCardiacCell)

#endif /*ABSTRACTBACKWARDEULERCARDIACCELL_HPP_*/
