
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

#ifndef ABSTRACTDYNAMICLINEARPDESOLVER_HPP_
#define ABSTRACTDYNAMICLINEARPDESOLVER_HPP_

#include "TimeStepper.hpp"
#include "AbstractLinearPdeSolver.hpp"
#include "PdeSimulationTime.hpp"
#include "AbstractTimeAdaptivityController.hpp"
#include "Hdf5DataWriter.hpp"
#include "Hdf5ToVtkConverter.hpp"

/**
 * Abstract class for dynamic linear PDE solves.
 * This class defines the Solve() method. The concrete class should implement
 * the SetupLinearSystem() method (defined in AbstractLinearPdeSolver), based
 * on the PDE being solved and the numerical method.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractDynamicLinearPdeSolver : public AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
    friend class TestSimpleLinearParabolicSolver;
protected:

    /** Simulation start time. */
    double mTstart;

    /** Simulation end time. */
    double mTend;

    /** Whether SetTimes has been called with suitable parameters. */
    bool mTimesSet;

    /** The initial condition vector. */
    Vec mInitialCondition;

   /** Whether the matrix has been assembled for the current time step. */
    bool mMatrixIsAssembled;

    /**
     * Whether the matrix is constant in time (if so the system need not be assembled at each time step).
     * Defaults to false.
     */
    bool mMatrixIsConstant;

    /**
     * The timestep to use. This is either the last timestep passed in in SetTimeStep,
     * or the last timestep suggested by the time adaptivity controller.
     */
    double mIdealTimeStep;

    /** The last actual timestep used. */
    double mLastWorkingTimeStep;

    /** A controller which determines what timestep to use (defaults to NULL). */
    AbstractTimeAdaptivityController* mpTimeAdaptivityController;

    /**
     * Flag to say if we need to output to Meshalyzer.
     * Defaults to false in the constructor.
     */
    bool mOutputToMeshalyzer;

    /**
     * Flag to say if we need to output to VTK.
     * Defaults to false in the constructor.
     */
    bool mOutputToVtk;

    /**
     * Flag to say if we need to output to VTK parallel (.pvtu).
     * Defaults to false in the constructor.
     */
    bool mOutputToParallelVtk;

    /** Output directory (a subfolder of tmp/[USERNAME]/testoutput). */
    std::string mOutputDirectory;

    /** Filename prefix for HDF5 and other files. */
    std::string mFilenamePrefix;

    /**
     * The ratio of the number of actual timesteps to the number
     * of timesteps at which results are output to HDF5 and other files.
     * Defaults to 1 in the constructor.
     */
    unsigned mPrintingTimestepMultiple;

    /** The object used to write results to HDF5 file. */
    Hdf5DataWriter* mpHdf5Writer;

    /** List of variable column IDs as written to HDF5 file. */
    std::vector<int> mVariableColumnIds;

    /**
     * Create and initialise the HDF5 writer.
     * Called by Solve() if results are to be output.
     */
    void InitialiseHdf5Writer();

    /**
     * Write one timestep of output data to HDF5 file.
     *
     * @param time  the time
     * @param solution  the solution vector to write
     */
    void WriteOneStep(double time, Vec solution);

public:

    /**
     * Constructor.
     *
     * @param pMesh the mesh
     */
    AbstractDynamicLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * Set the times to solve between.
     *
     * @param tStart the start time
     * @param tEnd the end time
     */
    void SetTimes(double tStart, double tEnd);

    /**
     * Set (or reset) the timestep to use.
     *
     * @param dt timestep
     */
    void SetTimeStep(double dt);

    /**
     * Set the initial condition.
     *
     * @param initialCondition the initial condition
     */
    void SetInitialCondition(Vec initialCondition);

    /** Dynamic solve method. */
    Vec Solve();

    /** Tell the solver to assemble the matrix again next timestep. */
    void SetMatrixIsNotAssembled();

    /**
     * Set a controller class which alters the dt used.
     *
     * @param pTimeAdaptivityController the controller
     */
    void SetTimeAdaptivityController(AbstractTimeAdaptivityController* pTimeAdaptivityController);

    /**
     * @param output whether to output to Meshalyzer files
     */
    void SetOutputToMeshalyzer(bool output);

    /**
     * @param output whether to output to VTK (.vtu) file
     */
    void SetOutputToVtk(bool output);

    /**
     * @param output whether to output to parallel VTK (.pvtu) file
     */
    void SetOutputToParallelVtk(bool output);

    /**
     * @param outputDirectory the output directory
     * @param prefix the filename prefix
     */
    void SetOutputDirectoryAndPrefix(std::string outputDirectory, std::string prefix);

    /**
     * @param multiple the ratio of the number of actual timesteps to the number
     * of timesteps at which results are output to HDF5 and other files.
     */
    void SetPrintingTimestepMultiple(unsigned multiple);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseHdf5Writer()
{
	// Check that everything is set up correctly
	if ((mOutputDirectory=="") || (mFilenamePrefix==""))
	{
        EXCEPTION("Output directory or filename prefix has not been set");
	}

	// Create writer
	mpHdf5Writer = new Hdf5DataWriter(*(this->mpMesh)->GetDistributedVectorFactory(),
                                      mOutputDirectory,
                                      mFilenamePrefix);

	// Set writer to output all nodes
	mpHdf5Writer->DefineFixedDimension((this->mpMesh)->GetNumNodes());

	// Only used to get an estimate of the number of timesteps below
	unsigned estimated_num_printing_timesteps = 1 + (mTend - mTstart)/(mIdealTimeStep*mPrintingTimestepMultiple);

	/**
	 * \todo allow user to specify units of time and names and units of
	 * dependent variables using DefineVariable() and DefineUnlimitedDimension()
	 * (#1841)
	 */
	assert(mVariableColumnIds.empty());
	for (unsigned i=0; i<PROBLEM_DIM; i++)
	{
		std::stringstream variable_name;
		variable_name << "Variable_" << i;
		mVariableColumnIds.push_back(mpHdf5Writer->DefineVariable(variable_name.str(),"undefined"));
	}
	mpHdf5Writer->DefineUnlimitedDimension("Time", "undefined", estimated_num_printing_timesteps);

    // End the define mode of the writer
	mpHdf5Writer->EndDefineMode();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractDynamicLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh),
      mTimesSet(false),
      mInitialCondition(NULL),
      mMatrixIsAssembled(false),
      mMatrixIsConstant(false),
      mIdealTimeStep(-1.0),
      mLastWorkingTimeStep(-1),
      mpTimeAdaptivityController(NULL),
      mOutputToMeshalyzer(false),
      mOutputToVtk(false),
      mOutputToParallelVtk(false),
      mOutputDirectory(""),
      mFilenamePrefix(""),
      mPrintingTimestepMultiple(1),
      mpHdf5Writer(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimes(double tStart, double tEnd)
{
    mTstart = tStart;
    mTend = tEnd;

    if (mTstart >= mTend)
    {
        EXCEPTION("Starting time has to less than ending time");
    }

    mTimesSet = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimeStep(double dt)
{
    if (dt <= 0)
    {
        EXCEPTION("Time step has to be greater than zero");
    }

    mIdealTimeStep = dt;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetInitialCondition(Vec initialCondition)
{
    assert(initialCondition != NULL);
    mInitialCondition = initialCondition;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::WriteOneStep(double time, Vec solution)
{
	mpHdf5Writer->PutUnlimitedVariable(time);
	if (PROBLEM_DIM == 1)
	{
		mpHdf5Writer->PutVector(mVariableColumnIds[0], solution);
	}
	else
	{
		mpHdf5Writer->PutStripedVector(mVariableColumnIds, solution);
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve()
{
	// Begin by checking that everything has been set up correctly
	if (!mTimesSet)
	{
		EXCEPTION("SetTimes() has not been called");
	}
	if ((mIdealTimeStep <= 0.0) && (mpTimeAdaptivityController==NULL))
	{
		EXCEPTION("SetTimeStep() has not been called");
	}
	if (mInitialCondition == NULL)
	{
		EXCEPTION("SetInitialCondition() has not been called");
	}

	// If required, initialise HDF5 writer and output initial condition to HDF5 file
	bool print_output = (mOutputToMeshalyzer || mOutputToVtk || mOutputToParallelVtk);
	if (print_output)
	{
		InitialiseHdf5Writer();
		WriteOneStep(mTstart, mInitialCondition);
		mpHdf5Writer->AdvanceAlongUnlimitedDimension();
	}

    this->InitialiseForSolve(mInitialCondition);

    if (mIdealTimeStep < 0) // hasn't been set, so a controller must have been given
    {
        mIdealTimeStep = mpTimeAdaptivityController->GetNextTimeStep(mTstart, mInitialCondition);
    }

    /*
     * Note: we use the mIdealTimeStep here (the original timestep that was passed in, or
     * the last timestep suggested by the controller), rather than the last timestep used
     * (mLastWorkingTimeStep), because the timestep will be very slightly altered by the
     * stepper in the final timestep of the last printing-timestep-loop, and these floating
     * point errors can add up and eventually cause exceptions being thrown.
     */
    TimeStepper stepper(mTstart, mTend, mIdealTimeStep, mMatrixIsConstant);

    Vec solution = mInitialCondition;
    Vec next_solution;

    while (!stepper.IsTimeAtEnd())
    {
        bool timestep_changed = false;

        PdeSimulationTime::SetTime(stepper.GetTime());

        // Determine timestep to use
        double new_dt;
        if (mpTimeAdaptivityController)
        {
            // Get the timestep the controller wants to use and store it as the ideal timestep
            mIdealTimeStep = mpTimeAdaptivityController->GetNextTimeStep(stepper.GetTime(), solution);

            // Tell the stepper to use this timestep from now on...
            stepper.ResetTimeStep(mIdealTimeStep);

            // ..but now get the timestep from the stepper, as the stepper might need
            // to trim the timestep if it would take us over the end time
            new_dt = stepper.GetNextTimeStep();

            timestep_changed = (fabs(mLastWorkingTimeStep-new_dt) > 1e-8);
        }
        else
        {
            new_dt = stepper.GetNextTimeStep();
        }

        // Save the timestep as the last one use, and also put it in PdeSimulationTime
        // so everyone can see it
        mLastWorkingTimeStep = new_dt;
        PdeSimulationTime::SetPdeTimeStep(new_dt);

        // Solve

        this->PrepareForSetupLinearSystem(solution);

        bool compute_matrix = (!mMatrixIsConstant || !mMatrixIsAssembled || timestep_changed);
//        if (compute_matrix) std::cout << " ** ASSEMBLING MATRIX!!! ** " << std::endl;
        this->SetupLinearSystem(solution, compute_matrix);

        this->FinaliseLinearSystem(solution);

        if (compute_matrix)
        {
            this->mpLinearSystem->ResetKspSolver();
        }

        next_solution = this->mpLinearSystem->Solve(solution);

        if (mMatrixIsConstant)
        {
            mMatrixIsAssembled = true;
        }

        this->FollowingSolveLinearSystem(next_solution);

        stepper.AdvanceOneTimeStep();

        // Avoid memory leaks
        if (solution != mInitialCondition)
        {
            HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
            VecDestroy(solution);
            HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
        }
        solution = next_solution;

        // If required, output next solution to HDF5 file
        if (print_output)
        {
        	WriteOneStep(stepper.GetTime(), solution);
            mpHdf5Writer->AdvanceAlongUnlimitedDimension();
        }
    }

    // Avoid memory leaks
    if (mpHdf5Writer != NULL)
    {
        delete mpHdf5Writer;
        mpHdf5Writer = NULL;
    }

    // Convert HDF5 output to other formats as required

///\todo allow output to Meshalyzer (#1841)
//	if (mOutputToMeshalyzer)
//    {
//        Hdf5ToMeshalyzerConverter<ELEMENT_DIM,SPACE_DIM> converter(mOutputDirectory, mFilenamePrefix, mpMesh);
//    }
    if (mOutputToVtk)
    {
        Hdf5ToVtkConverter<ELEMENT_DIM,SPACE_DIM> converter(mOutputDirectory, mFilenamePrefix, this->mpMesh, false, false);
    }
    if (mOutputToParallelVtk)
    {
        Hdf5ToVtkConverter<ELEMENT_DIM,SPACE_DIM> converter(mOutputDirectory, mFilenamePrefix, this->mpMesh, true, false);
    }

    return solution;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetMatrixIsNotAssembled()
{
	mMatrixIsAssembled = false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimeAdaptivityController(AbstractTimeAdaptivityController* pTimeAdaptivityController)
{
    assert(pTimeAdaptivityController != NULL);
    assert(mpTimeAdaptivityController == NULL);
    mpTimeAdaptivityController = pTimeAdaptivityController;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputToMeshalyzer(bool output)
{
	mOutputToMeshalyzer = output;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputToVtk(bool output)
{
	mOutputToVtk = output;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputToParallelVtk(bool output)
{
	mOutputToParallelVtk = output;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputDirectoryAndPrefix(std::string outputDirectory, std::string prefix)
{
	mOutputDirectory = outputDirectory;
	mFilenamePrefix = prefix;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetPrintingTimestepMultiple(unsigned multiple)
{
	mPrintingTimestepMultiple = multiple;
}

#endif /*ABSTRACTDYNAMICLINEARPDESOLVER_HPP_*/
