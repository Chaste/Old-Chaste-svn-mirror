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

#include "ChasteSerialization.hpp"
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
    friend class TestPeregoCellModels;
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
    
    /**
     * This function should never be called - the cell class incorporates its own solver.
     *
     * @param time
     * @param rY
     * @param rDY
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
     /**
     * Computes the predictor step of the scheme
     *
     * @param rSolutionAtPreviousTime is the solution of the state variables at the previous time step
     * @param rPredictedSolution is returned vector with the predicted values of the gates and other state variables (this is assumed to be initialised to the correct size)
     * @param currentTime is the current time
     */
    void EvaluatePredictedValues(const std::vector<double>& rSolutionAtPreviousTime, std::vector<double>& rPredictedSolution, double currentTime);
    
    /**
     * Computes the corrector step of the scheme
     *
     * @param rPredictedSolution is input vector with the predicted values of the gates and other state variables 
     * @param rCorrectedSolution is the corrected solution of the state variables at the previous time step (this is assumed to be initialised to the correct size)
     * @param currentTime is the current time
     */
    void EvaluateCorrectedValues(const std::vector<double>& rPredictedSolution, std::vector<double>& rCorrectedSolution, double currentTime);
    
    /**
     * Computes the error between the predicted solution and the corrected solution
     *
     * @param rErrors is the outpue errors (assumed to be initialised to the correct size)
     * @param rPredictedSolution is the predicted values of the gates and other state variables 
     * @param rCorrectedSolution is the corrected solution of the state variables at the previous time step
     * @param currentTime is the current time
     */
    void EvaluateErrors(std::vector<double>& rErrors, const std::vector<double>& rPredictedSolution, const std::vector<double>& rCorrectedSolution, double currentTime);
    
   
  protected:
    
    /** 
     * 
     * The local time step
     * Initially set to the mDt in the AbstractOdeSystem class.
     * It is then adjusted in the adaptive part of the algorithm 
     * and used to calculate the gates at the following time step
     * 
     * */
    double mLocalTimeStep;
    
    std::vector<unsigned> mGatingVariableIndices; /**< Indices of those variables associated with gates (not concentrations or voltages etc.) */
    std::vector<double> mSolutionAtPreviousTimeStep; /**< Cache of previous solution */
    //Nomenclature of variables is based on the paper (see class documentation for reference) 
    
    std::vector<double> mCorrectedSolution;/**< Variable to store the corrected solution. It's needed by the compute method that will use it as "previous" solution */
    
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
    
    std::vector<double> ma_error;/**< value for system variables needed for error evaluation, first part*/
    std::vector<double> mb_error;/**< value for system variables needed for error evaluation, second part*/
    
    bool mIsTheCorrectorStep; /**< An helper boolean to flag whether we are in the corrector step*/
    bool mIsTheFirstStep; /**< A helper boolean to indicate whether we have taken any timesteps yet*/
    bool mIsTheErrorEvaluationStep; /**< A helper boolean to flag whether we're in the error evaluation step*/
    
    double mThetaP; /**< A numerical solver parameter which changes as dt does*/
    double mThetaC; /**< Another numerical solver parameter which changes as dt does*/
    
    double mNewDt;  /**< Variable to store the new timestep size until it is copied into mLocalTimestep*/
    double mNewDtFromEndOfPreviousPdeStep; /**< Variable to store the timestep that would've been used at the end of the last PDE timestep if it hadn't caused overshoot of endTime. */ 
    
    std::vector<double> mWeightedErrorTolerances; /**< Vector of tolerances to error in each of the system variables, will be weighted by a small tolerance factor*/
    
    bool mUseAdaptTimestep; /**< For testing purposes, so we can test the algorithm without adaptivity. To be removed eventually. */
    
    unsigned mNumberOfStateVariables; /**< stores the number of state variables (for memory allocation)*/ 
    /**
     * Return true if the error in any variable exceeds tolerances
     * 
     * @param rErrors vector of the errors that will be checked against the tolerances
     * 
     */
    bool IsThereTooMuchError(std::vector<double>& rErrors);
     
    bool mIsThereTooMuchError; /**< To hold the return value of IsThereTooMuchError as it is required multiple times.*/
    
    /** 
     * Change the timestep size for the next step based upon the a posteriori errors
     * 
     * @param rErrors a vector of errors on which the error analysis to determine whether to adapt or not is performed
     * 
     */ 
    void AdaptTimestep(std::vector<double>& rErrors); 
      
    /**
     * Computes some parameters needed by the Perego Veneziani algorithm.
     * Implemented in the child class. 
     *
     * @param stateVariablesAtPrevousTime is the solution of the state variables at the previous time step
     * @param currentTime is the current time
     */
    virtual void ComputeSystemParameters(const std::vector<double>& stateVariablesAtPrevousTime, double currentTime)=0;
    
    /**
     * Adjusts the system parameters after a change in the timestep size. This depends on the ratio
     * of the old timestep to the new.
     * 
     * @param oldDt is the old timestep
     * @param newDt is the new timestep to be set
     */ 
    void ChangeTimestepAndRecomputeParameters(double oldDt, double newDt);    
         
public:

    /**
     * Standard constructor for a cell.
     *
     * @param numberOfStateVariables  the size of the ODE system
     * @param voltageIndex  the index of the variable representing the transmembrane
     *     potential within the state variable vector
     * @param pIntracellularStimulus  the intracellular stimulus function
     * @param useAdaptTimestep For testing purposes, so we can test the algorithm without adaptivity. To be removed eventually.
     */
    AbstractPeregoCardiacCell(
        unsigned numberOfStateVariables,
        unsigned voltageIndex,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus,
        bool useAdaptTimestep=true);

    /** Virtual destructor */
    virtual ~AbstractPeregoCardiacCell();
   
    
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
        double saved_voltage = GetVoltage();
        //PRINT_VARIABLE(tStart);
//        mSetVoltageDerivativeToZero = true;
        // Recover the ODE timestep size from the end of the last PDE timestep
        mNewDt = mNewDtFromEndOfPreviousPdeStep;
        Compute(tStart, tEnd);
//        mSetVoltageDerivativeToZero = false;
    
        SetVoltage(saved_voltage);
        
      //  PRINT_VECTOR(this->mStateVariables);
    
        VerifyStateVariables();        
//        NEVER_REACHED;
        // not tested in tissue yet
    }
    
    /**
     * Set method to switch on and off the time adaptivity
     * 
     * @param flag is true if you want adaptivity to be on, false otherwise
     */
    void  SetAdaptivityFlag (bool flag);
    
};



#endif // _ABSTRACTPEREGOCARDIACCELL_HPP_
