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

#include "CvodeAdaptor.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "TimeStepper.hpp"

#include <iostream>

// CVODE headers
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <cvode/cvode_dense.h>

/**
 * A helper function to copy an N_Vector into a std::vector.
 */
void CopyToStdVector(N_Vector src, std::vector<realtype>& rDest)
{
    // Check for no-op
    if (NV_DATA_S(src) == &(rDest[0])) return;
    // Set dest size
    rDest.resize(NV_LENGTH_S(src));
    // Copy data
    for (long i=0; i<NV_LENGTH_S(src); i++)
    {
        rDest[i] = NV_Ith_S(src, i);
    }
}

int CvodeRhsAdaptor(realtype t, N_Vector y, N_Vector ydot, void* pData)
{
    assert(pData != NULL);
    CvodeData* p_data = (CvodeData*) pData;
    // Get y, ydot into std::vector<>s
    static std::vector<realtype> ydot_vec;
    CopyToStdVector(y, *p_data->pY);
    CopyToStdVector(ydot, ydot_vec);
    // Call our function
    try
    {
        p_data->pSystem->EvaluateYDerivatives(t, *(p_data->pY), ydot_vec);
    }
    catch (Exception &e)
    {
        std::cerr << "CVODE RHS Exception: " << e.GetMessage() << std::endl << std::flush;
        return -1;
    }
    // Copy derivative back
    for (long i=0; i<NV_LENGTH_S(ydot); i++)
        NV_Ith_S(ydot, i) = ydot_vec[i];
    return 0;
}

int CvodeRootAdaptor(realtype t, N_Vector y, realtype* pGOut, void* pData)
{
    assert(pData != NULL);
    CvodeData* p_data = (CvodeData*) pData;
    // Get y into a std::vector
    CopyToStdVector(y, *p_data->pY);
    // Call our function
    try
    {
        *pGOut = p_data->pSystem->CalculateRootFunction(t, *p_data->pY);
    }
    catch (Exception &e)
    {
        std::cerr << "CVODE Root Exception: " << e.GetMessage() << std::endl << std::flush;
        return -1;
    }
    return 0;
}

// int CvodeDenseJacobianAdaptor(long int numberOfStateVariables, DenseMat J,
//                               realtype t, N_Vector y, N_Vector fy,
//                               void* pData,
//                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
// {
//     AbstractOdeSystemWithAnalyticJacobian* pSystem
//         = (AbstractOdeSystemWithAnalyticJacobian*) pData;
//     // Get the current time step
//     double dt;
//     CVodeGetCurrentStep(CvodeMem, &dt);
//     // Get std::vector<> for y and double** for J
//     std::vector<realtype>& y_vec = *NV_DATA_STL(y);
//     double** ppJ = J->data; // organised column-wise: J_{i,j} = ppJ[j][i]
//     // Call our function
//     try
//     {
//         pSystem->AnalyticJacobian(y_vec, ppJ, t, dt);
//     }
//     catch (Exception &e)
//     {
//         std::cerr << "CVODE J Exception: " << e.GetMessage() << std::endl << std::flush;
//         return -1;
//     }
//     // Update J (if needed)
//     return 0;
// }


void CvodeErrorHandler(int errorCode, const char *module, const char *function,
                       char *message, void* pData)
{
    std::stringstream err;
    err << "CVODE Error " << errorCode << " in module " << module
        << " function " << function << ": " << message;
    std::cerr << "*" << err.str() << std::endl << std::flush;
    EXCEPTION(err.str());
}



void CvodeAdaptor::SetupCvode(AbstractOdeSystem* pOdeSystem,
                              std::vector<double>& rInitialY,
                              double startTime, double maxStep)
{
    assert(rInitialY.size() == pOdeSystem->GetNumberOfStateVariables());
    assert(maxStep > 0.0);

    mInitialValues = N_VMake_Serial(rInitialY.size(), &(rInitialY[0]));
    assert(NV_DATA_S(mInitialValues) == &(rInitialY[0]));
    assert(!NV_OWN_DATA_S(mInitialValues));
//    std::cout << " Initial values: "; N_VPrint_Serial(mInitialValues);
//    std::cout << " Rtol: " << mRelTol << ", Atol: " << mAbsTol << std::endl;
//    std::cout << " Start: " << startTime << " max dt=" << maxStep << std::endl << std::flush;

    mpCvodeMem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (mpCvodeMem == NULL) EXCEPTION("Failed to SetupCvode CVODE");
    // Set error handler
    CVodeSetErrHandlerFn(mpCvodeMem, CvodeErrorHandler, NULL);
    // Set the user data
    mData.pSystem = pOdeSystem;
    mData.pY = &rInitialY;
    CVodeSetFdata(mpCvodeMem, (void*)(&mData));
    // Setup CVODE
    CVodeMalloc(mpCvodeMem, CvodeRhsAdaptor, startTime, mInitialValues,
                CV_SS, mRelTol, &mAbsTol);
    CVodeSetMaxStep(mpCvodeMem, maxStep);
    // Set the rootfinder function if wanted
    if (mCheckForRoots)
    {
        CVodeRootInit(mpCvodeMem, 1, CvodeRootAdaptor, (void*)(&mData));
    }
    // Change max steps if wanted
    if (mMaxSteps > 0)
    {
        CVodeSetMaxNumSteps(mpCvodeMem, mMaxSteps);
    }
    // Attach a linear solver for Newton iteration
    CVDense(mpCvodeMem, rInitialY.size());
}

void CvodeAdaptor::FreeCvodeMemory()
{
//    assert(!NV_OWN_DATA_STL(mInitialValues));
//    std::vector<double>* pVec = NV_DATA_STL(mInitialValues);
//    double val = (*pVec)[0];
    N_VDestroy_Serial(mInitialValues); mInitialValues = NULL;
//    std::cout << "  a: " << val << ", b: " << (*pVec)[0] << std::endl;

    CVodeFree(&mpCvodeMem);
}


#define COVERAGE_IGNORE
void CvodeAdaptor::CvodeError(int flag, const char * msg)
{
    std::stringstream err;
    err << msg << ": " << CVodeGetReturnFlagName(flag);
    std::cerr << err.str() << std::endl << std::flush;
    EXCEPTION(err.str());
}
#undef COVERAGE_IGNORE


OdeSolution CvodeAdaptor::Solve(AbstractOdeSystem* pOdeSystem,
                                std::vector<double>& rYValues,
                                double startTime,
                                double endTime,
                                double maxStep,
                                double timeSampling)
{
    assert(endTime > startTime);
    assert(timeSampling > 0.0);

    mStoppingEventOccurred = false;
    if (mCheckForRoots && pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true)
    {
        EXCEPTION("(Solve with sampling) Stopping event is true for initial condition");
    }

    SetupCvode(pOdeSystem, rYValues, startTime, maxStep);

    TimeStepper stepper(startTime, endTime, timeSampling);
    N_Vector yout = mInitialValues;

    // Set up ODE solution
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    solutions.SetOdeSystemInformation(pOdeSystem->GetSystemInformation());

    // Main time sampling loop
    while (!stepper.IsTimeAtEnd() && !mStoppingEventOccurred)
    {
//        std::cout << "Solving to time " << stepper.GetNextTime() << std::endl << std::flush;
        double tend;
        int ierr;
        try
        {
            ierr = CVode(mpCvodeMem, stepper.GetNextTime(), yout,
                         &tend, CV_NORMAL);
        }
        catch (...)
        {
            FreeCvodeMemory();
            throw;
        }
        if (ierr<0) CvodeError(ierr, "CVODE failed to solve system");
        // Store solution
        solutions.rGetSolutions().push_back(rYValues);
        solutions.rGetTimes().push_back(tend);
        if (ierr == CV_ROOT_RETURN)
        {
            // Stopping event occurred
            mStoppingEventOccurred = true;
            mStoppingTime = tend;
        }
        stepper.AdvanceOneTimeStep();
    }

    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTotalTimeStepsTaken());

//    std::cout << " Solved to " << stepper.GetTime() << " in " << stepper.GetTotalTimeStepsTaken() << " samples.\n" << std::flush;
    int ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning
    FreeCvodeMemory();

    return solutions;
}


void CvodeAdaptor::Solve(AbstractOdeSystem* pOdeSystem,
                         std::vector<double>& rYValues,
                         double startTime,
                         double endTime,
                         double maxStep)
{
    assert(endTime > startTime);

    mStoppingEventOccurred = false;
    if (mCheckForRoots && pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true)
    {
        EXCEPTION("(Solve) Stopping event is true for initial condition");
    }

    SetupCvode(pOdeSystem, rYValues, startTime, maxStep);

    N_Vector yout = mInitialValues;
    double tend;
    int ierr;
    try
    {
        ierr = CVode(mpCvodeMem, endTime, yout, &tend, CV_NORMAL);
    }
    catch (...)
    {
        FreeCvodeMemory();
        throw;
    }
    if (ierr<0) CvodeError(ierr, "CVODE failed to solve system");
    if (ierr == CV_ROOT_RETURN)
    {
        // Stopping event occurred
        mStoppingEventOccurred = true;
        mStoppingTime = tend;
//        std::cout << "CVODE Stopped at t = " << tend << std::endl;
    }
    assert(NV_DATA_S(yout) == &(rYValues[0]));
    assert(!NV_OWN_DATA_S(yout));

//    long int steps;
//    CVodeGetNumSteps(mpCvodeMem, &steps);
//    std::cout << " Solved to " << endTime << " in " << steps << " steps.\n";

    ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
//    if (mStoppingEventOccurred)
//    {
//        std::cout << "Last internal dt was " << mLastInternalStepSize << std::endl;
//    }
    assert(ierr == CV_SUCCESS);
    ierr=ierr; // avoid unused var warning
    FreeCvodeMemory();
}

CvodeAdaptor::CvodeAdaptor(double relTol, double absTol)
    : AbstractIvpOdeSolver(),
      mpCvodeMem(NULL), mInitialValues(NULL),
      mRelTol(relTol), mAbsTol(absTol),
      mLastInternalStepSize(-0.0),
      mMaxSteps(0),
      mCheckForRoots(false)
{
}

void CvodeAdaptor::SetTolerances(double relTol, double absTol)
{
    mRelTol = relTol;
    mAbsTol = absTol;
}

double CvodeAdaptor::GetRelativeTolerance()
{
    return mRelTol;
}

double CvodeAdaptor::GetAbsoluteTolerance()
{
    return mAbsTol;
}

double CvodeAdaptor::GetLastStepSize()
{
    return mLastInternalStepSize;
}

void CvodeAdaptor::CheckForStoppingEvents()
{
    mCheckForRoots = true;
}

void CvodeAdaptor::SetMaxSteps(long int numSteps)
{
    mMaxSteps = numSteps;
}

long int CvodeAdaptor::GetMaxSteps()
{
    return mMaxSteps;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CvodeAdaptor);

#endif // CHASTE_CVODE
