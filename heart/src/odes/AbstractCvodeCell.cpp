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

#ifdef CHASTE_CVODE

#include <sstream>
#include <iostream>
#include <cmath>

#include "AbstractCvodeCell.hpp"
#include "CvodeAdaptor.hpp" // For CvodeErrorHandler
#include "TimeStepper.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"
#include "VectorHelperFunctions.hpp"

// CVODE headers
#include <cvode/cvode.h>
#include <sundials/sundials_nvector.h>
#include <cvode/cvode_dense.h>

/**
 * Callback function provided to CVODE to allow it to 'call' C++ member functions
 * (in particular, AbstractCvodeCell::EvaluateYDerivatives).
 *
 * @param t  current time
 * @param y  state variable vector
 * @param ydot  derivatives vector to be filled in
 * @param pData  pointer to the cell being simulated
 */
int AbstractCvodeCellRhsAdaptor(realtype t, N_Vector y, N_Vector ydot, void *pData)
{
    assert(pData != NULL);
    AbstractCvodeCell* pCell = (AbstractCvodeCell*) pData;
    try
    {
        pCell->EvaluateYDerivatives(t, y, ydot);
    }
    catch (const Exception &e)
    {
        std::cerr << "CVODE RHS Exception: " << e.GetMessage()
                  << std::endl << std::flush;
        return -1;
    }
    return 0;
}

AbstractCvodeCell::AbstractCvodeCell(boost::shared_ptr<AbstractIvpOdeSolver> /* unused */,
                                     unsigned numberOfStateVariables,
                                     unsigned voltageIndex,
                                     boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCellInterface(boost::shared_ptr<AbstractIvpOdeSolver>(), voltageIndex, pIntracellularStimulus),
      AbstractCvodeSystem(numberOfStateVariables),
      mpCvodeMem(NULL),
      mMaxSteps(0),
      mMaxDt(DOUBLE_UNSET)
{
    SetTolerances();
}


AbstractCvodeCell::~AbstractCvodeCell()
{
}

double AbstractCvodeCell::GetVoltage()
{
    assert(mStateVariables);
    return GetAnyVariable(mVoltageIndex);
}

void AbstractCvodeCell::SetVoltage(double voltage)
{
    assert(mStateVariables);
    SetAnyVariable(mVoltageIndex, voltage);
}

void AbstractCvodeCell::SetVoltageDerivativeToZero(bool clamp)
{
    mSetVoltageDerivativeToZero = clamp;
}


void AbstractCvodeCell::SetTimestep(double maxDt)
{
    mMaxDt = maxDt;
}


void AbstractCvodeCell::SolveAndUpdateState(double tStart, double tEnd)
{
    if (mMaxDt == DOUBLE_UNSET)
    {
        SetTimestep(HeartConfig::Instance()->GetPrintingTimeStep());
    }
    Solve(tStart, tEnd, mMaxDt);
}

OdeSolution AbstractCvodeCell::Compute(double tStart, double tEnd, double tSamp)
{
    if (tSamp == 0.0)
    {
        tSamp = HeartConfig::Instance()->GetPrintingTimeStep();
    }
    if (mMaxDt == DOUBLE_UNSET)
    {
        SetTimestep(tSamp);
    }
    return Solve(tStart, tEnd, mMaxDt, tSamp);
}


void AbstractCvodeCell::ComputeExceptVoltage(double tStart, double tEnd)
{
    EXCEPTION("This method is not yet implemented for CVODE cells.");
}


OdeSolution AbstractCvodeCell::Solve(realtype tStart,
                                     realtype tEnd,
                                     realtype maxDt,
                                     realtype tSamp)
{
    assert(tEnd >= tStart);
    assert(tSamp > 0.0);

    SetupCvode(mStateVariables, tStart, maxDt);

    TimeStepper stepper(tStart, tEnd, tSamp);

    // Set up ODE solution
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(MakeStdVec(mStateVariables));
    solutions.rGetTimes().push_back(tStart);
    solutions.SetOdeSystemInformation(mpSystemInfo);

    // Main time sampling loop
    while (!stepper.IsTimeAtEnd())
    {
        double cvode_stopped_at;
        int ierr = CVode(mpCvodeMem, stepper.GetNextTime(), mStateVariables,
                         &cvode_stopped_at, CV_NORMAL);
        if (ierr<0)
        {
            FreeCvodeMemory();
            CvodeError(ierr, "CVODE failed to solve system");
        }
        // Not root finding, so should have reached requested time
        assert(fabs(cvode_stopped_at - stepper.GetNextTime()) < DBL_EPSILON);

        VerifyStateVariables();

        // Store solution
        solutions.rGetSolutions().push_back(MakeStdVec(mStateVariables));
        solutions.rGetTimes().push_back(cvode_stopped_at);
        stepper.AdvanceOneTimeStep();
    }

    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTotalTimeStepsTaken());

    int ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning
    FreeCvodeMemory();

    return solutions;
}

void AbstractCvodeCell::Solve(realtype tStart,
                              realtype tEnd,
                              realtype maxDt)
{
    assert(tEnd >= tStart);

    SetupCvode(mStateVariables, tStart, maxDt);

    double cvode_stopped_at;
    int ierr = CVode(mpCvodeMem, tEnd, mStateVariables, &cvode_stopped_at, CV_NORMAL);
    if (ierr<0)
    {
        FreeCvodeMemory();
        CvodeError(ierr, "CVODE failed to solve system");
    }
    // Not root finding, so should have reached requested time
    assert(fabs(cvode_stopped_at - tEnd) < DBL_EPSILON);

    ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning
    FreeCvodeMemory();

    VerifyStateVariables();
}


void AbstractCvodeCell::SetMaxSteps(long int numSteps)
{
    mMaxSteps = numSteps;
}

long int AbstractCvodeCell::GetMaxSteps()
{
    return mMaxSteps;
}

void AbstractCvodeCell::SetTolerances(double relTol, double absTol)
{
    mRelTol = relTol;
    mAbsTol = absTol;
}

double AbstractCvodeCell::GetRelativeTolerance()
{
    return mRelTol;
}

double AbstractCvodeCell::GetAbsoluteTolerance()
{
    return mAbsTol;
}

double AbstractCvodeCell::GetLastStepSize()
{
    return mLastInternalStepSize;
}


void AbstractCvodeCell::SetupCvode(N_Vector initialConditions,
                                   realtype tStart,
                                   realtype maxDt)
{
    assert((unsigned)NV_LENGTH_S(initialConditions) == GetNumberOfStateVariables());
    assert(maxDt >= 0.0);

    mpCvodeMem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (mpCvodeMem == NULL) EXCEPTION("Failed to SetupCvode CVODE");
    // Set error handler
    CVodeSetErrHandlerFn(mpCvodeMem, CvodeErrorHandler, NULL);
    // Set the user data & setup CVODE
#if CHASTE_SUNDIALS_VERSION >= 20400
    CVodeSetUserData(mpCvodeMem, (void*)(this));
    CVodeInit(mpCvodeMem, AbstractCvodeCellRhsAdaptor, tStart, initialConditions);
    CVodeSStolerances(mpCvodeMem, mRelTol, mAbsTol);
#else
    CVodeSetFdata(mpCvodeMem, (void*)(this));
    CVodeMalloc(mpCvodeMem, AbstractCvodeCellRhsAdaptor, tStart, initialConditions,
                CV_SS, mRelTol, &mAbsTol);
#endif
    CVodeSetMaxStep(mpCvodeMem, maxDt);
    // Change max steps if wanted
    if (mMaxSteps > 0)
    {
        CVodeSetMaxNumSteps(mpCvodeMem, mMaxSteps);
    }
    // Attach a linear solver for Newton iteration
    CVDense(mpCvodeMem, NV_LENGTH_S(initialConditions));
}


void AbstractCvodeCell::FreeCvodeMemory()
{
    CVodeFree(&mpCvodeMem);
}


void AbstractCvodeCell::CvodeError(int flag, const char * msg)
{

    std::stringstream err;
    char* p_flag_name = CVodeGetReturnFlagName(flag);
    err << msg << ": " << p_flag_name;
    free(p_flag_name);
    std::cerr << err.str() << std::endl << std::flush;
    EXCEPTION(err.str());
}



#endif // CHASTE_CVODE
