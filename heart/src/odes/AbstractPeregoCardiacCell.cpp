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

#include "AbstractPeregoCardiacCell.hpp"
#include "TimeStepper.hpp"
#include "Debug.hpp"

AbstractPeregoCardiacCell::AbstractPeregoCardiacCell(unsigned numberOfStateVariables,
                                unsigned voltageIndex,
                                boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver>(),
                              numberOfStateVariables,
                              voltageIndex,
                              pIntracellularStimulus)
{
    mIsTheCorrectorStep = false;
    mIsTheFirstStep = true;
    mSolutionAtPreviousTimeStep.resize(numberOfStateVariables);
}


AbstractPeregoCardiacCell::~AbstractPeregoCardiacCell()
{}


void  AbstractPeregoCardiacCell::EvaluatePredictedValues(const std::vector<double>& rSolutionAtPreviousTime, std::vector<double>& rPredictedSolution, double currentTime)
{
    //Get the weigths from the end of the corrector step at the previous time
    //hardcoded for the moment. To be changed when error analysis is implemented
    mc0bar = 3.0/2.0;
    mc1bar = -0.5;
    
    mSolutionAtPreviousTimeStep = rSolutionAtPreviousTime;
    //Compute parameters (done in the child class, will modify the member variables here)
    ComputeSystemParameters(rSolutionAtPreviousTime, currentTime);

    bool it_is_a_gating_variable =false;
    //loop over the state variables
    for (unsigned i = 0; i<rSolutionAtPreviousTime.size();i++)
    {
        //reset the flag
        it_is_a_gating_variable=false;
        //apply the predictor scheme to the gating variables only
        for (unsigned j =0;j<mGatingVariableIndices.size();j++)
        {
            if (i == mGatingVariableIndices[j])
            {
                it_is_a_gating_variable = true;
            }
        }

        if (it_is_a_gating_variable)
        {
            double gate_derivative = (exp((mc0bar*ma_current[i]+mc1bar*ma_previous[i])*this->mDt)-1) / 
                                     ((mc0bar*ma_current[i] + mc1bar * ma_previous[i])*this->mDt)*((mc0bar*ma_current[i]+mc1bar*ma_previous[i])*rSolutionAtPreviousTime[i]+mc0bar*mb_current[i]+mc1bar*mb_previous[i]);
            rPredictedSolution[i]=rSolutionAtPreviousTime[i] + this->mDt*gate_derivative;
        }
        // All other variables are updated by a weighted version of the forward Euler:
        else
        {
            double variable_derivative = mc0bar * ma_current[i] + mc1bar * ma_previous[i];            
            
            rPredictedSolution[i] = rSolutionAtPreviousTime[i] + this->mDt * variable_derivative;
        }
    }

}

void  AbstractPeregoCardiacCell::EvaluateCorrectedValues(const std::vector<double>& rPredictedSolution, std::vector<double>& rCorrectedSolution, double currentTime)
{
    mIsTheCorrectorStep = true;
    
    //Get the weigths from the end of the corrector step at the previous time
    //hardcoded for the moment. To be changed when error analysis is implemented
    mc0 = 8.0/12.0;
    mc1 = -1.0/12.0;
    mcMinus1 =  5.0/12.0;

    //Compute parameters (done in the child class, will modify the member variables here)
    ComputeSystemParameters(rPredictedSolution, currentTime);

    bool it_is_a_gating_variable =false;
    //loop over the state variables
    for (unsigned i = 0; i<mSolutionAtPreviousTimeStep.size();i++)
    {
        //reset the flag
        it_is_a_gating_variable=false;
        //apply the corrector scheme to the gating variables only
        for (unsigned j =0;j<mGatingVariableIndices.size();j++)
        {
            if (i == mGatingVariableIndices[j])
            {
                it_is_a_gating_variable = true;
            }
        }

        if (it_is_a_gating_variable)
        {
            double gate_derivative = (exp((mcMinus1*ma_predicted[i] +  mc0*ma_current[i]+mc1*ma_previous[i])*this->mDt)-1)/((mcMinus1*ma_predicted[i]+mc0*ma_current[i] + mc1 * ma_previous[i])*this->mDt)*((mcMinus1*ma_predicted[i]+mc0*ma_current[i]+mc1*ma_previous[i])*mSolutionAtPreviousTimeStep[i]+ mcMinus1*mb_predicted[i] + mc0*mb_current[i]+mc1*mb_previous[i]);
            rCorrectedSolution[i]=mSolutionAtPreviousTimeStep[i] + this->mDt*gate_derivative;
        }
        // All other variables are updated by a weighted version of the forward Euler:
        else
        {
            double variable_derivative = mcMinus1*ma_predicted[i] + mc0 * ma_current[i] + mc1 * ma_previous[i];
            rCorrectedSolution[i]=mSolutionAtPreviousTimeStep[i] + this->mDt * variable_derivative;
        }
    }
    //corrector step is over
    mIsTheCorrectorStep = false;
}

void  AbstractPeregoCardiacCell::EvaluateErrors(std::vector<double>& rErrors, const std::vector<double>& rPredictedSolution, const std::vector<double>& rCorrectedSolution, double currentTime)
{
    mIsTheErrorEvaluationStep = true;
    
    //Compute parameters (done in the child class, will modify the member variables here)
    ComputeSystemParameters(rCorrectedSolution, currentTime);
    mIsTheErrorEvaluationStep = false;
    
    bool it_is_a_gating_variable = false;
    
    // These are fixed for now, but will change with dt when the adaptive timestep size is implemented
    mThetaP = -0.5;
    mThetaC = 1.0/3.0; // Set this to 2.0/3.0 in order to see errors (for testing purposes).   
    
    double error_weight_factor = (mThetaC - 1.0/3.0) / (mThetaP - mThetaC);   
    
    //loop over the state variables
    for (unsigned i = 0; i<rPredictedSolution.size();i++)
    {
        //reset the flag
        it_is_a_gating_variable=false;
        //apply the predictor scheme to the gating variables only
        for (unsigned j = 0; j<mGatingVariableIndices.size(); j++)
        {
            if (i == mGatingVariableIndices[j])
            {
                it_is_a_gating_variable = true;
            }
        }

        if (it_is_a_gating_variable)
        {

            double gate_error =  fabs(error_weight_factor * ( rCorrectedSolution[i] - rPredictedSolution[i]) + (1.0/12.0) * (ma_error[i] * mb_previous[i] - ma_previous[i] * mb_error[i]) * this->mDt * this->mDt);
            rErrors[i] = gate_error;
        }
        // All other variables are updated by a weighted version of the forward Euler:
        else
        {
            double variable_error = fabs( error_weight_factor * ( rCorrectedSolution[i] - rPredictedSolution[i] ) );            
            rErrors[i] = variable_error;

            
        }
        
    }
   
    
    
}

void AbstractPeregoCardiacCell::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    NEVER_REACHED;
}

OdeSolution AbstractPeregoCardiacCell::Compute(double startTime, double endTime)
{
    TimeStepper stepper(startTime, endTime, this->mDt);

    // setup solutions if output is required
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(this->mStateVariables);
    solutions.rGetTimes().push_back(startTime);
    
    std::vector<double> previous_yvalues = mStateVariables;
    std::vector<double> predicted_values = mStateVariables;
    std::vector<double> corrected_values = mStateVariables;
    std::vector<double> errors = mStateVariables;

 
    // Solve the ODE system
    while ( !stepper.IsTimeAtEnd() )
    {
        // Initialise the previous state variables to be a copy of the current state variables
        previous_yvalues = corrected_values;
        
        //predict the next value
        EvaluatePredictedValues(previous_yvalues, predicted_values, stepper.GetTime());
        EvaluateCorrectedValues(predicted_values, corrected_values, stepper.GetTime());
        
        // Evaluate the error in the variables using their predictor and corrector values
                
        EvaluateErrors(errors, predicted_values, corrected_values, stepper.GetTime());
        

        
        
        stepper.AdvanceOneTimeStep();
        
        // Record the state variables after the timestep is taken:
        solutions.rGetSolutions().push_back(corrected_values);
        solutions.rGetTimes().push_back(stepper.GetTime());        
    }
    return solutions;
}

