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

#include "AbstractCardiacProblem.hpp"

#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"
#include "HeartEventHandler.hpp"
#include "TimeStepper.hpp"
#include "PetscTools.hpp"
#include "DistributedVector.hpp"
#include "ProgressReporter.hpp"
#include "LinearSystem.hpp"
#include "PostProcessingWriter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "Hdf5ToVtkConverter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractCardiacProblem(
            AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory)
    : mMeshFilename(""),               // i.e. undefined
      mNodesPerProcessorFilename(""),  // i.e. undefined
      mUseMatrixBasedRhsAssembly(true),
      mAllocatedMemoryForMesh(false),
      mWriteInfo(false),
      mPrintOutput(true),
      mpCardiacPde(NULL),
      mpAssembler(NULL),
      mpCellFactory(pCellFactory),
      mpMesh(NULL),
      mSolution(NULL),
      mCurrentTime(0.0),
      mpWriter(NULL)
{
    assert(mNodesToOutput.empty());

    HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractCardiacProblem()
    // it doesn't really matter what we initialise these to, as they'll be overwritten by
    // the serialization methods
    : mMeshFilename(""),
      mNodesPerProcessorFilename(""),
      mUseMatrixBasedRhsAssembly(true),
      mAllocatedMemoryForMesh(false), // Handled by AbstractCardiacPde
      mWriteInfo(false),
      mPrintOutput(true),
      mVoltageColumnId(UINT_MAX),
      mTimeColumnId(UINT_MAX),
      mNodeColumnId(UINT_MAX),
      mpCardiacPde(NULL),
      mpAssembler(NULL),
      mpCellFactory(NULL),
      mpMesh(NULL),
      mSolution(NULL),
      mCurrentTime(0.0),
      mpWriter(NULL)
{
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~AbstractCardiacProblem()
{
    delete mpCardiacPde;
    if (mSolution)
    {
        VecDestroy(mSolution);
    }

    if (mAllocatedMemoryForMesh)
    {
        delete mpMesh;
    }
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::Initialise()
{
    HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
    if (mpMesh==NULL)
    {
        // If no mesh has been passed, we get it from the configuration file
        try
        {
            if(HeartConfig::Instance()->GetLoadMesh())
            {
                TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(HeartConfig::Instance()->GetMeshName());
                if (ELEMENT_DIM == 1)
                {
                    ///\todo We CAN currently instantiate the parallel mesh in 1D, but there's an archiving test which assumes that 1D meshes are sequential
                    mpMesh = new TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>();
                }
                else
                {
                    mpMesh = new ParallelTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>();
                }
                mpMesh->ConstructFromMeshReader(mesh_reader);
            }
            else
            {
                if(HeartConfig::Instance()->GetCreateMesh())
                {
                    assert(HeartConfig::Instance()->GetSpaceDimension()==SPACE_DIM);
                    double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();

                    switch(HeartConfig::Instance()->GetSpaceDimension())
                    {
                        case 1:
                        {
                            c_vector<double, 1> fibre_length;
                            HeartConfig::Instance()->GetFibreLength(fibre_length);

                            mpMesh = new TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>();

                            unsigned slab_nodes_x = (unsigned)round(fibre_length[0]/inter_node_space);

                            static_cast<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(mpMesh)->ConstructLinearMesh(slab_nodes_x);

                            break;
                        }
                        case 2:
                        {
                            c_vector<double, 2> sheet_dimensions; //cm
                            HeartConfig::Instance()->GetSheetDimensions(sheet_dimensions);

                            mpMesh = new TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>();

                            unsigned slab_nodes_x = (unsigned)round(sheet_dimensions[0]/inter_node_space);
                            unsigned slab_nodes_y = (unsigned)round(sheet_dimensions[1]/inter_node_space);

                            static_cast<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(mpMesh)->ConstructRectangularMesh(slab_nodes_x, slab_nodes_y);

                            break;
                        }
                        case 3:
                        {
                            c_vector<double, 3> slab_dimensions; //cm
                            HeartConfig::Instance()->GetSlabDimensions(slab_dimensions);

                            ///\todo This could be parallel execpt for a acceptance test which expects to be able to post-process.
                            mpMesh = new TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>();

                            unsigned slab_nodes_x = (unsigned)round(slab_dimensions[0]/inter_node_space);
                            unsigned slab_nodes_y = (unsigned)round(slab_dimensions[1]/inter_node_space);
                            unsigned slab_nodes_z = (unsigned)round(slab_dimensions[2]/inter_node_space);
                            ///\todo Do we still need a cast here?
                            ///\todo This could be parallel execpt for a acceptance test which expects to be able to post-process.
                            static_cast<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(mpMesh)->ConstructCuboid(slab_nodes_x,
                                                   slab_nodes_y,
                                                   slab_nodes_z);
                            break;
                        }
                        default:
                            NEVER_REACHED;
                    }
                    // scale
                    double mesh_scale_factor = inter_node_space;
                    mpMesh->Scale(mesh_scale_factor, mesh_scale_factor, mesh_scale_factor);

                }
                else
                {
                    NEVER_REACHED;
                }
            }

            mAllocatedMemoryForMesh = true;
        }
        catch (Exception& e)
        {
            EXCEPTION(std::string("No mesh given: define it in XML parameters file or call SetMesh()\n") + e.GetShortMessage());
        }
    }
    mpCellFactory->SetMesh( mpMesh );
    
    // if the user requested transmural stuff, we fill in the mCellHeterogeneityAreas here.
    if (HeartConfig::Instance()->AreCellularTransmuralHeterogeneitiesRequested())
    {
        mpCellFactory->FillInCellularTransmuralAreas();
    }
    
    if (mNodesPerProcessorFilename != "")
    {
        mpMesh->ReadNodesPerProcessorFile(mNodesPerProcessorFilename);
    }
    HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);

    ///\todo Should this method be rolled into the Solve() method or the PreSolveChecks()?
    delete mpCardiacPde; // In case we're called twice
    mpCardiacPde = CreateCardiacPde();
    ///\todo The above line isn't accounted for by the event handler, and can be expensive (e.g. conductivity tensors)

    // Delete any previous solution, so we get a fresh initial condition
    if (mSolution)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        VecDestroy(mSolution);
        mSolution = NULL;
        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
    }

    // Always start at time zero
    mCurrentTime = 0.0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetNodesPerProcessorFilename(const std::string& filename)
{
    mNodesPerProcessorFilename = filename;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetBoundaryConditionsContainer(boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> > pBcc)
{
    this->mpBoundaryConditionsContainer = pBcc;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PreSolveChecks()
{
    if ( mpCardiacPde == NULL ) // if pde is NULL, Initialise() probably hasn't been called
    {
        EXCEPTION("Pde is null, Initialise() probably hasn't been called");
    }
    if ( HeartConfig::Instance()->GetSimulationDuration() <= mCurrentTime)
    {
        EXCEPTION("End time should be in the future");
    }
    if (mPrintOutput==true)
    {
        if( (HeartConfig::Instance()->GetOutputDirectory()=="") || (HeartConfig::Instance()->GetOutputFilenamePrefix()==""))
        {
            EXCEPTION("Either explicitly specify not to print output (call PrintOutput(false)) or specify the output directory and filename prefix");
        }
    }

    double end_time = HeartConfig::Instance()->GetSimulationDuration();
    double pde_time = HeartConfig::Instance()->GetPdeTimeStep();

    // MatrixIsConstant stuff requires CONSTANT dt - do some checks to make sure the TimeStepper won't find
    // non-constant dt.
    // Note: printing_time does not have to divide end_time, but dt must divide printing_time and end_time.
    // HeartConfig checks pde_dt divides printing dt
    if( fabs( end_time - pde_time*round(end_time/pde_time)) > 1e-10 )
    {
        EXCEPTION("Pde timestep does not seem to divide end time - check parameters");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::CreateInitialCondition()
{
    DistributedVectorFactory* p_factory = mpMesh->GetDistributedVectorFactory();
    Vec initial_condition = p_factory->CreateVec(PROBLEM_DIM);
    DistributedVector ic = p_factory->CreateDistributedVector(initial_condition);
    std::vector<DistributedVector::Stripe> stripe;
    stripe.reserve(PROBLEM_DIM);

    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        stripe.push_back(DistributedVector::Stripe(ic, i));
    }

    for (DistributedVector::Iterator index = ic.Begin();
         index != ic.End();
         ++index)
    {
        stripe[0][index] = mpCardiacPde->GetCardiacCell(index.Global)->GetVoltage();
        if (PROBLEM_DIM==2)
        {
            stripe[1][index] = 0;
        }
    }

    ic.Restore();

    return initial_condition;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    // If this fails the mesh has already been set. We assert rather throw an exception
    // to avoid a memory leak when checking it throws correctly
    assert(mpMesh==NULL);
    assert(pMesh!=NULL);
    mAllocatedMemoryForMesh = false;
    mpMesh = pMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintOutput(bool printOutput)
{
    mPrintOutput = printOutput;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetWriteInfo(bool writeInfo)
{
    mWriteInfo = writeInfo;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetSolution()
{
    return mSolution;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
DistributedVector AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetSolutionDistributedVector()
{
    return mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(mSolution);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM> & AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::rGetMesh()
{
    assert (mpMesh);
    return *mpMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacPde<ELEMENT_DIM,SPACE_DIM>* AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetPde()
{
    return mpCardiacPde;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::Solve()
{
    PreSolveChecks();

    if (!mpBoundaryConditionsContainer) // the user didn't supply a bcc
    {
        // Set up the default bcc
        mpDefaultBoundaryConditionsContainer.reset(new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>);
        for (unsigned problem_index=0; problem_index<PROBLEM_DIM; problem_index++)
        {
            mpDefaultBoundaryConditionsContainer->DefineZeroNeumannOnMeshBoundary(mpMesh, problem_index);
        }
        mpBoundaryConditionsContainer = mpDefaultBoundaryConditionsContainer;
    }

    assert(mpAssembler==NULL);
    mpAssembler = CreateAssembler(); // passes mpBoundaryConditionsContainer to assember

    // If we have already run a simulation, use the old solution as initial condition
    Vec initial_condition;
    if (mSolution)
    {
        initial_condition = mSolution;
    }
    else
    {
        initial_condition = CreateInitialCondition();
    }

    TimeStepper stepper(mCurrentTime,
                        HeartConfig::Instance()->GetSimulationDuration(),
                        HeartConfig::Instance()->GetPrintingTimeStep());

    std::string progress_reporter_dir;

    if (mPrintOutput)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
        InitialiseWriter();

        // If we are resuming a simulation (i.e. mSolution already exists) there's no need
        // of writing the initial timestep,
        // since it was already written as the last timestep of the previous run
        if (!mSolution)
        {
            WriteOneStep(stepper.GetTime(), initial_condition);
            mpWriter->AdvanceAlongUnlimitedDimension();
        }
        HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);

        progress_reporter_dir = HeartConfig::Instance()->GetOutputDirectory();
    }
    else
    {
        progress_reporter_dir = ""; // progress printed to CHASTE_TEST_OUTPUT
    }

    // create a progress reporter so users can track how much has gone and
    // estimate how much time is left. (Note this has to be done after the
    // InitialiseWriter above (if mPrintOutput==true)
    ProgressReporter progress_reporter(progress_reporter_dir, mCurrentTime,
                                       HeartConfig::Instance()->GetSimulationDuration());
    progress_reporter.Update(0);


    while ( !stepper.IsTimeAtEnd() )
    {
        // solve from now up to the next printing time
        mpAssembler->SetTimes(stepper.GetTime(), stepper.GetNextTime(), HeartConfig::Instance()->GetPdeTimeStep());
        mpAssembler->SetInitialCondition( initial_condition );

        try
        {
            mSolution = mpAssembler->Solve();
        }
        catch (Exception &e)
        {
            // Free memory.
            delete mpAssembler;
            mpAssembler=NULL;
            //VecDestroy(initial_condition);
#ifndef NDEBUG
            PetscTools::ReplicateException(true);
#endif
            // Re-throw
            HeartEventHandler::Reset();//EndEvent(HeartEventHandler::EVERYTHING);

            CloseFilesAndPostProcess();
            throw e;
        }
#ifndef NDEBUG
        PetscTools::ReplicateException(false);
#endif

        // Free old initial condition
        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        VecDestroy(initial_condition);
        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

        // Initial condition for next loop is current solution
        initial_condition = mSolution;

        // update the current time
        stepper.AdvanceOneTimeStep();
        mCurrentTime = stepper.GetTime();

        if (mPrintOutput)
        {
            // print out details at current time if asked for
            if (mWriteInfo)
            {
                HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
                WriteInfo(stepper.GetTime());
                HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
            }

            // Writing data out to the file <FilenamePrefix>.dat
            HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
            WriteOneStep(stepper.GetTime(), mSolution);
            // Just flags that we've finished a time-step; won't actually 'extend' unless new data is written.
            mpWriter->AdvanceAlongUnlimitedDimension();

            HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
        }

		progress_reporter.Update(stepper.GetTime());

        OnEndOfTimestep(stepper.GetTime());
    }

    // Free assembler
    delete mpAssembler;
    mpAssembler=NULL;

    // close the file that stores voltage values
    progress_reporter.PrintFinalising();
    CloseFilesAndPostProcess();
    HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::CloseFilesAndPostProcess()
{
    // Close files
    if (!mPrintOutput)
    {
        //Nothing to do
        return;
    }
    mpWriter->Close();
    delete mpWriter;

    HeartEventHandler::BeginEvent(HeartEventHandler::USER2); //Temporarily using USER2 to instrument post-processing
    // Only if results files were written and we are outputting all nodes
    if (mNodesToOutput.empty())
    {
        if (HeartConfig::Instance()->GetVisualizeWithMeshalyzer())
        {
            //Convert simulation data to Meshalyzer format
            Hdf5ToMeshalyzerConverter<ELEMENT_DIM,SPACE_DIM> converter(HeartConfig::Instance()->GetOutputDirectory(), HeartConfig::Instance()->GetOutputFilenamePrefix(), mpMesh);
        }

        if (HeartConfig::Instance()->GetVisualizeWithCmgui())
        {
            //Convert simulation data to Cmgui format
            Hdf5ToCmguiConverter<ELEMENT_DIM,SPACE_DIM> converter(HeartConfig::Instance()->GetOutputDirectory(), HeartConfig::Instance()->GetOutputFilenamePrefix(), mpMesh);
        }

        if (HeartConfig::Instance()->GetVisualizeWithVtk())
        {

            //Convert simulation data to Cmgui format
            Hdf5ToVtkConverter<ELEMENT_DIM,SPACE_DIM> converter(HeartConfig::Instance()->GetOutputDirectory(), HeartConfig::Instance()->GetOutputFilenamePrefix(), mpMesh);
        }
    }

    if(HeartConfig::Instance()->IsPostProcessingRequested())
    {
        //Test that we have a tetrahedral mesh
        TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* p_tetmesh= dynamic_cast<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(mpMesh);

//        if (p_tetmesh == NULL)
//        {
//            //Assume that it's a parallel tetrahedral mesh - Note that every process throws together
//            ///\todo need to sort out issues with parallel meshes
//            EXCEPTION("Cannot post-process on a parallel mesh yet");
//        }
        assert(p_tetmesh != NULL);
        PostProcessingWriter<ELEMENT_DIM, SPACE_DIM> post_writer(*p_tetmesh, HeartConfig::Instance()->GetOutputDirectory(),
                        HeartConfig::Instance()->GetOutputFilenamePrefix(), true);
        post_writer.WritePostProcessingFiles();
    }

    HeartEventHandler::EndEvent(HeartEventHandler::USER2); //Temporarily using USER2 to instrument post-processing
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineWriterColumns(bool extending)
{
    if (!extending)
    {
        if ( mNodesToOutput.empty() )
        {
            //Set writer to output all nodes
            mpWriter->DefineFixedDimension( mpMesh->GetNumNodes() );
        }
        else
        {
            //Output only the nodes indicted
            mpWriter->DefineFixedDimension( mNodesToOutput, mpMesh->GetNumNodes() );
        }
        //mNodeColumnId = mpWriter->DefineVariable("Node", "dimensionless");
        mVoltageColumnId = mpWriter->DefineVariable("V","mV");

        mpWriter->DefineUnlimitedDimension("Time","msecs");
    }
    else
    {
        mVoltageColumnId = mpWriter->GetVariableByName("V");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineExtraVariablesWriterColumns(bool extending)
{
    // Check if any extra output variables have been requested
    if (HeartConfig::Instance()->GetOutputVariablesProvided())
    {
        // Get their names in a vector
        std::vector<std::string> output_variables;
        HeartConfig::Instance()->GetOutputVariables(output_variables);

        // Loop over them
        for (unsigned var_index=0; var_index<output_variables.size(); var_index++)
        {
            // Get variable name
            std::string var_name = output_variables[var_index];

            // Register it (or look it up) in the data writer
            unsigned column_id;
            if (extending)
            {
                column_id = this->mpWriter->GetVariableByName(var_name);
            }
            else
            {
                column_id = this->mpWriter->DefineVariable(var_name, "");
            }

            // Store column id
            mExtraVariablesId.push_back(column_id);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::WriteExtraVariablesOneStep()
{
    // Loop over the requested state variables
    for (unsigned var_index=0; var_index<mExtraVariablesId.size(); var_index++)
    {
        // Create vector for storing values over the local nodes
        Vec variable_data =  this->mpMesh->GetDistributedVectorFactory()->CreateVec();
        DistributedVector distributed_var_data = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(variable_data);

        // Loop over the local nodes and gather the data
        for (DistributedVector::Iterator index = distributed_var_data.Begin();
             index!= distributed_var_data.End();
             ++index)
        {
            // Store value for node "index"
            distributed_var_data[index] = this->mpCardiacPde->GetCardiacCell(index.Global)->GetStateVariableValueByNumber(mExtraVariablesId[var_index]);
        }
        distributed_var_data.Restore();

        // Write it to disc
        this->mpWriter->PutVector(mExtraVariablesId[var_index], variable_data);

        VecDestroy(variable_data);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::InitialiseWriter()
{
    bool extend_file = (mSolution != NULL);

    mpWriter = new Hdf5DataWriter(*mpMesh->GetDistributedVectorFactory(),
                                  HeartConfig::Instance()->GetOutputDirectory(),
                                  HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                  !extend_file, // don't clear directory if extension requested
                                  extend_file);

    // Define columns, or get the variable IDs from the writer
    DefineWriterColumns(extend_file);
    if (!extend_file)
    {
        mpWriter->EndDefineMode();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetOutputNodes(std::vector<unsigned> &nodesToOutput)
{
    mNodesToOutput = nodesToOutput;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Hdf5DataReader AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDataReader()
{
    if( (HeartConfig::Instance()->GetOutputDirectory()=="") || (HeartConfig::Instance()->GetOutputFilenamePrefix()==""))
    {
        EXCEPTION("Data reader invalid as data writer cannot be initialised");
    }
    return Hdf5DataReader(HeartConfig::Instance()->GetOutputDirectory(), HeartConfig::Instance()->GetOutputFilenamePrefix());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::UseMatrixBasedRhsAssembly(bool usematrix)
{
    mUseMatrixBasedRhsAssembly = usematrix;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

// Monodomain
template class AbstractCardiacProblem<1,1,1>;
template class AbstractCardiacProblem<1,2,1>;
template class AbstractCardiacProblem<1,3,1>;
template class AbstractCardiacProblem<2,2,1>;
template class AbstractCardiacProblem<3,3,1>;

// Bidomain
template class AbstractCardiacProblem<1,1,2>;
template class AbstractCardiacProblem<2,2,2>;
template class AbstractCardiacProblem<3,3,2>;
