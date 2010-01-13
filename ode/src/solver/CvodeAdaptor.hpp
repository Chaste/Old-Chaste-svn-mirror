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
#ifndef _CVODEADAPTOR_HPP_
#define _CVODEADAPTOR_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractIvpOdeSolver.hpp"
#include "OdeSolution.hpp"

// CVODE headers
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <cvode/cvode_dense.h>

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * CVODE right-hand-side function adaptor.
 *
 * The CVODE solvers require the RHS of the ODE system to be defined
 * by a function of this type.  We use pData to access the Chaste ODE system,
 * and call EvaluateYDerivatives appropriately.
 *
 * Note that this requires copying the state variable and derivatives vectors,
 * and thus introduces a slight overhead.
 *
 * @param t  the current time
 * @param y  the current state variable values
 * @param ydot  to be filled in with the derivatives
 * @param pData  a pointer to the AbstractOdeSystem to evaluate
 * @return 0 on success, -1 for an unrecoverable error
 */
int CvodeRhsAdaptor(realtype t, N_Vector y, N_Vector ydot, void* pData);

/**
 * CVODE root-finder function adaptor.
 *
 * Adapt the Chaste AbstractOdeSystem::CalculateStoppingEvent method for use by CVODE.
 *
 * This function computes a vector-valued function g(t, y) such that the roots of the
 * components g_i(t, y) are to be found during the integration.
 *
 * Unfortunately, AbstractOdeSystem::CalculateStoppingEvent returns a boolean value,
 * so we have to cheat in the definition of g.
 *
 * Note that this function requires copying the state variable vector, and thus
 * introduces a slight overhead.
 *
 * @param t  the current time
 * @param y  the current state variable values
 * @param pGOut  pointer to array to be filled in with the g_i(t, y) values
 * @param pData  a pointer to the AbstractOdeSystem to use
 * @return 0 on success, negative on error
 */
int CvodeRootAdaptor(realtype t, N_Vector y, realtype* pGOut, void* pData);

// /**
//  * Jacobian computation adaptor function.
//  *
//  * If solving an AbstractOdeSystemWithAnalyticJacobian, this function
//  * can be used to allow CVODE to compute the Jacobian analytically.
//  *
//  * Note to self: can test using pSystem->GetUseAnalyticJacobian().
//  */
// int CvodeDenseJacobianAdaptor(long int numberOfStateVariables, DenseMat J,
//                               realtype t, N_Vector y, N_Vector fy,
//                               void* pData,
//                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/**
 * CVODE error handling function.
 *
 * Throw an Exception to report errors, rather than the CVODE approach of magic
 * return codes.
 */
void CvodeErrorHandler(int errorCode, const char *module, const char *function,
                       char *message, void* pData);


typedef struct CvodeData_ {
    std::vector<realtype>* pY;
    AbstractOdeSystem* pSystem;
} CvodeData;

/**
 * The CVODE adaptor ODE solver class.
 *
 * Assumes that it will be solving stiff systems, so uses BDF/Newton.
 *
 * The timeStep parameters of the abstract class are here used to specify
 * *maximum* steps, since the solver is adaptive.
 *
 * Note that a call to Solve will initialise the CVODE solver, and free its
 * working memory when done.  There is thus a non-trivial overhead involved.
 *
 * \todo Add an option to just initialise once, and assume subsequent Solve
 *   calls are continuing from where we left off.
 */
class CvodeAdaptor : public AbstractIvpOdeSolver
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class - all member variables instantiated on construction or temporary.
        archive & boost::serialization::base_object<AbstractIvpOdeSolver>(*this);
        archive & mRelTol;
        archive & mAbsTol;
        // All other member variables given values on each call.
    }

    /** Pointer to the CVODE memory block. */
    void* mpCvodeMem;

    /** Initial conditions for the ODE solver. */
    N_Vector mInitialValues;

    /** The CVODE data structure. */
    CvodeData mData;

    /** Relative tolerance for the ODE solver. */
    double mRelTol;

    /** Absolute tolerance for the ODE solver. */
    double mAbsTol;

    /** The size of the previous timestep. */
    double mLastInternalStepSize;

    /**
     * The maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int mMaxSteps;

    /** Whether to check for stopping events. */
    bool mCheckForRoots;

protected:

    /**
     * Set up the CVODE data structures needed to solve the given system.
     */
    void SetupCvode(AbstractOdeSystem* pOdeSystem,
                    std::vector<double>& rInitialY,
                    double startTime, double maxStep);

    /**
     * Free CVODE memory after solving a system of ODEs.
     */
    void FreeCvodeMemory();

    /**
     * Report an error from CVODE.
     *
     * This will (probably) never be called, since we supply an error handler function
     * which throws an exception.
     */
    void CvodeError(int flag, const char * msg);

public:

    /**
     * Default constructor.
     * Can optionally set relative and absolute tolerances.
     *
     * @param relTol the relative tolerance for the solver (defaults to 1e-4)
     * @param absTol the absolute tolerance for the solver (defaults to 1e-6)
     */
    CvodeAdaptor(double relTol=1e-4, double absTol=1e-6);

    /**
     * Set relative and absolute tolerances; both scalars.
     * If no parameters are given, tolerances will be reset to default values.
     *
     * @param relTol the relative tolerance for the solver (defaults to 1e-4)
     * @param absTol the absolute tolerance for the solver (defaults to 1e-6)
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
     * Get the last step size used internally by CVODE in the last Solve call
     */
    double GetLastStepSize();

    /**
     * Solve the given ODE system, returning the solution at sampling intervals.
     *
     * @param pOdeSystem  the ODE system to solve
     * @param rYValues  the initial state variable values
     *   (note: this vector will also be used as working memory)
     * @param startTime  the time to start solving at
     * @param endTime  the time to solve to
     * @param maxStep  the maximum time step to be taken by the adaptive solver
     * @param timeSampling  the interval at which to sample the solution
     * @return  the solution
     */
    OdeSolution Solve(AbstractOdeSystem* pOdeSystem,
                      std::vector<double>& rYValues,
                      double startTime,
                      double endTime,
                      double maxStep,
                      double timeSampling);

    /**
     * Solve the given ODE system, storing the final result in rYValues.
     *
     * @param pOdeSystem  the ODE system to solve
     * @param rYValues  the initial state variable values; will be filled in with
     *   the final values on return
     * @param startTime  the time to start solving at
     * @param endTime  the time to solve to
     * @param maxStep  the maximum time step to be taken by the adaptive solver
     */
    void Solve(AbstractOdeSystem* pOdeSystem,
               std::vector<double>& rYValues,
               double startTime,
               double endTime,
               double maxStep);

    /**
     * Make the solver check for stopping events using CVODE's rootfinding functionality.
     *
     * By default we do not check.
     */
    void CheckForStoppingEvents();

    /**
     * Change the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.  Default is 500.
     */
    void SetMaxSteps(long int numSteps);

    /**
     * Get the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int GetMaxSteps();

};

#include "TemplatedExport.hpp"
CHASTE_CLASS_EXPORT(CvodeAdaptor);

#endif // _CVODEADAPTOR_HPP_
#endif // CHASTE_CVODE
