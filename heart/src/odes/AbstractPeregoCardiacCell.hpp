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
 * Reference:
 * Perego M, Veneziani A. An efficient generalisation of the Rush-Larsen method for solving electrophysiology membrane equations
 * Technical Report TR-2009-005.
 * www.mathcs.emory.edu/technical-reports/techrep-00153.pdf
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
    void EvaluatePredictedValues(std::vector<double> solutionAtPreviousTime, std::vector<double>& rPredictedSolution, double currentTime);
    
    
    void EvaluateCorrectedValues(std::vector<double> predictedSolution, std::vector<double>& rCorrectedSolution, double currentTime);
    
    /**
     * Computes some parameters needed by the Perego Veneziani algorithm.
     * Implemented in the child class. 
     *
     * @param stateVariablesAtPrevousTime is the solution of the state variables at the previous time step
     * @param currentTime is the current time
     */
    virtual void ComputeSystemParameters(std::vector<double> stateVariablesAtPrevousTime, double currentTime)=0;
    
    /**
     * Overloaded Compute method.
     * 
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    OdeSolution Compute(double tStart, double tEnd);

    /**
     * Overloaded ComputeExceptVoltage method.
     * 
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt, but does not update the voltage.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    void ComputeExceptVoltage(double tStart, double tEnd)
    {
        NEVER_REACHED; // not tested in tissue yet
    }
    
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
    std::vector<double> mSolutionAtPreviousTimeStep;
    //Nomenclature of variables is based on the paper (see class documentation for reference) 
    
    double mc0bar;/**< weights for Adams-Bashforth integration*/
    double mc1bar;/**< weights for Adams-Bashforth integration*/
    
    double mc0;/**< weights for Adams-Moulton integration*/
    double mc1;/**< weights for Adams-Moulton integration*/
    double mcMinus1;/**< weights for Adams-Moulton integration*/
    
    std::vector<double> ma_current;/**< current value for system variables, first part*/
    std::vector<double> ma_predicted;/**< predicted value for system variables, first part*/
    std::vector<double> ma_previous;/**< value for system variables, first part, from previous time step*/
    std::vector<double> mb_current;/**< current value for system variables, second part*/
    std::vector<double> mb_predicted;/**< predicted value for system variables, second part*/
    std::vector<double> mb_previous;/**< value for system variables, second part, from previous time step*/
    
    bool mIsTheCorrectorStep; /**< An helper boolean to flag whether we are in the corrector step*/
    bool mIsTheFirstStep; /**< A helper boolean to indicate whether we have taken any timesteps yet*/

};



#endif // _ABSTRACTPEREGOCARDIACCELL_HPP_
