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
#include <cassert>

#include "AbstractCvodeSystem.hpp"
#include "Exception.hpp"
#include "VectorHelperFunctions.hpp"
#include "TimeStepper.hpp"
#include "CvodeAdaptor.hpp" // For CvodeErrorHandler

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
int AbstractCvodeSystemRhsAdaptor(realtype t, N_Vector y, N_Vector ydot, void *pData)
{
    assert(pData != NULL);
    AbstractCvodeSystem* p_ode_system = (AbstractCvodeSystem*) pData;
    try
    {
        p_ode_system->EvaluateYDerivatives(t, y, ydot);
    }
    catch (const Exception &e)
    {
        std::cerr << "CVODE RHS Exception: " << e.GetMessage()
                  << std::endl << std::flush;
        return -1;
    }
    return 0;
}

AbstractCvodeSystem::AbstractCvodeSystem(unsigned numberOfStateVariables)
    : AbstractParameterisedSystem<N_Vector>(numberOfStateVariables),
      mUseAnalyticJacobian(false),
      mpCvodeMem(NULL),
      mMaxSteps(0)
{
    SetTolerances();
}

void AbstractCvodeSystem::Init()
{
    DeleteVector(mStateVariables);
    mStateVariables = GetInitialConditions();
    DeleteVector(mParameters);
    mParameters = N_VNew_Serial(rGetParameterNames().size());
    for (int i=0; i<NV_LENGTH_S(mParameters); i++)
    {
        NV_Ith_S(mParameters, i) = 0.0;
    }
}

AbstractCvodeSystem::~AbstractCvodeSystem()
{
    DeleteVector(mStateVariables);
    DeleteVector(mParameters);
}

//
//double AbstractCvodeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
//{
//    bool stop = CalculateStoppingEvent(time, rY);
//    return stop ? 0.0 : 1.0;
//}
//
//bool AbstractCvodeSystem::GetUseAnalyticJacobian()
//{
//    return mUseAnalyticJacobian;
//}



OdeSolution AbstractCvodeSystem::Solve(realtype tStart,
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

void AbstractCvodeSystem::Solve(realtype tStart,
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


void AbstractCvodeSystem::SetMaxSteps(long int numSteps)
{
    mMaxSteps = numSteps;
}

long int AbstractCvodeSystem::GetMaxSteps()
{
    return mMaxSteps;
}

void AbstractCvodeSystem::SetTolerances(double relTol, double absTol)
{
    mRelTol = relTol;
    mAbsTol = absTol;
}

double AbstractCvodeSystem::GetRelativeTolerance()
{
    return mRelTol;
}

double AbstractCvodeSystem::GetAbsoluteTolerance()
{
    return mAbsTol;
}

double AbstractCvodeSystem::GetLastStepSize()
{
    return mLastInternalStepSize;
}


void AbstractCvodeSystem::SetupCvode(N_Vector initialConditions,
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
    CVodeInit(mpCvodeMem, AbstractCvodeSystemRhsAdaptor, tStart, initialConditions);
    CVodeSStolerances(mpCvodeMem, mRelTol, mAbsTol);
#else
    CVodeSetFdata(mpCvodeMem, (void*)(this));
    CVodeMalloc(mpCvodeMem, AbstractCvodeSystemRhsAdaptor, tStart, initialConditions,
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


void AbstractCvodeSystem::FreeCvodeMemory()
{
    CVodeFree(&mpCvodeMem);
}


void AbstractCvodeSystem::CvodeError(int flag, const char * msg)
{

    std::stringstream err;
    char* p_flag_name = CVodeGetReturnFlagName(flag);
    err << msg << ": " << p_flag_name;
    free(p_flag_name);
    std::cerr << err.str() << std::endl << std::flush;
    EXCEPTION(err.str());
}


#endif // CHASTE_CVODE
