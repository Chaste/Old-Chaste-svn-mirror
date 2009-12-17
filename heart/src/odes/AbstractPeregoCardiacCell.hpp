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


#ifndef _ABSTRACTPEREGOCARDIACCELL_HPP_
#define _ABSTRACTPEREGOCARDIACCELL_HPP_

#include <boost/serialization/access.hpp>
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCardiacCell.hpp"

#include <cassert>
#include <cmath>

/**
 * This is the base class for cardiac cells solved using a Perego Veneziani predictor corrector scheme
 *
 */
class AbstractPeregoCardiacCell : public AbstractCardiacCell
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
     */
    AbstractPeregoCardiacCell(
        unsigned numberOfStateVariables,
        unsigned voltageIndex,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /** Virtual destructor */
    virtual ~AbstractPeregoCardiacCell();
    
    /**
     * Computes the predictor step of the scheme
     *
     * @param solutionAtPreviousTime is the solution of the state variables at the previous time step
     * @param rPredictedSolution is returned vector with the predicted values of the gates and other state variables (this is assumed to be initialised to the correct size)
     * @param currentTime is the current time
     */
    void EvaluatePredictedGates(std::vector<double> solutionAtPreviousTime, std::vector<double>& rPredictedSolution, double currentTime);
    
    /**
     * Computes some parameters needed by the Perego Veneziani algorithm.
     * Implemented in the child class. 
     *
     * @param stateVariablesAtPrevousTime is the solution of the state variables at the previous time step
     * @param currentTime is the current time
     */
    virtual void ComputeSystemParameters(std::vector<double> stateVariablesAtPrevousTime, double currentTime)=0;

private:
    /**
     * This function should never be called - the cell class incorporates its own solver.
     *
     * @param time
     * @param rY
     * @param rDY
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);


protected:

    std::vector<unsigned> mGatingVariableIndices; /**< Indices of those variables associated with gates (not concentrations or voltages etc.) */
    double mc0bar;/**< weights for Adams-Bashforth integration*/
    double mc1bar;/**< weights for Adams-Bashforth integration*/
    std::vector<double> ma_current;/**< current value for system variables, first part*/
    std::vector<double> ma_previous;/**< value for system variables, first part, from previous time step*/
    std::vector<double> mb_current;/**< current value for system variables, second part*/
    std::vector<double> mb_previous;/**< value for system variables, second part, from previous time step*/
};



#endif // _ABSTRACTPEREGOCARDIACCELL_HPP_
