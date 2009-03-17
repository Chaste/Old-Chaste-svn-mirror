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

#include "AbstractCardiacProblem.hpp"

#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"
#include "HeartEventHandler.hpp"
#include "TimeStepper.hpp"
#include "PetscTools.hpp"
#include "DistributedVector.hpp"
#include "ProgressReporter.hpp"

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::AbstractCardiacProblem(
            AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
    : mMeshFilename(""),               // i.e. undefined
      mNodesPerProcessorFilename(""),  // i.e. undefined
      mOutputDirectory(""),            // i.e. undefined
      mOutputFilenamePrefix(""),       // i.e. undefined
      mUseMatrixBasedRhsAssembly(true),
      mpBoundaryConditionsContainer(NULL),
      mpDefaultBoundaryConditionsContainer(NULL),
      mpCellFactory(pCellFactory),
      mpMesh(NULL),
      mpWriter(NULL)
{
    mWriteInfo = false;
    mPrintOutput = true;
    mCallChaste2Meshalyzer = false;
    mpCardiacPde = NULL;
    mpAssembler = NULL;
    mSolution = NULL;
    mAllocatedMemoryForMesh = false;
    assert(mNodesToOutput.empty());

    HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::~AbstractCardiacProblem()
{
    if(mpDefaultBoundaryConditionsContainer!=NULL)
    {
        delete mpDefaultBoundaryConditionsContainer;
    }
    
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

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::Initialise()
{
    mOutputDirectory = HeartConfig::Instance()->GetOutputDirectory();
    mOutputFilenamePrefix = HeartConfig::Instance()->GetOutputFilenamePrefix();
    
    if (mpMesh==NULL)
    {
        // If no mesh has been passed, we get it from the configuration file
        try
        {
            /// \todo: Only considering \<LoadMesh/> definition. Consider \<Slab/> too
            if(HeartConfig::Instance()->GetLoadMesh())
            {
                TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(HeartConfig::Instance()->GetMeshName());
                if (SPACE_DIM == 1)
                {
                    ///\todo We can't currently instantiate the parallel mesh in 1D
                    mpMesh = new TetrahedralMesh<SPACE_DIM, SPACE_DIM>();
                }
                else
                {
                    mpMesh = new ParallelTetrahedralMesh<SPACE_DIM, SPACE_DIM>();
                }
        
                HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
                mpMesh->ConstructFromMeshReader(mesh_reader);
                HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
            }
            else
            {
                if(HeartConfig::Instance()->GetCreateMesh())
                {
                    switch(HeartConfig::Instance()->GetSpaceDimension())
                    {
                        case 1:
                        {
                            c_vector<double, 1> fibre_length;
                            HeartConfig::Instance()->GetFibreLength(fibre_length);
                            double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();

                            mpMesh = new TetrahedralMesh<SPACE_DIM, SPACE_DIM>();
                            
                            unsigned slab_nodes_x = (unsigned)round(fibre_length[0]/inter_node_space);
                        
                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->ConstructLinearMesh(slab_nodes_x);
                            // place at origin
//                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->Translate(-(double)slab_nodes_x/2.0,
//                                             -(double)slab_nodes_y/2.0,
//                                             -(double)slab_nodes_z/2.0);
                                             
                            // scale
                            double mesh_scale_factor = inter_node_space;
                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->Scale(mesh_scale_factor, 1, 1);
                            break;
                        }
                        case 2:
                        {
                            c_vector<double, 2> sheet_dimensions; //cm                    
                            HeartConfig::Instance()->GetSheetDimensions(sheet_dimensions);
                            double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();
                            
                            mpMesh = new TetrahedralMesh<SPACE_DIM, SPACE_DIM>();
                            
                            unsigned slab_nodes_x = (unsigned)round(sheet_dimensions[0]/inter_node_space);
                            unsigned slab_nodes_y = (unsigned)round(sheet_dimensions[1]/inter_node_space);
                        
                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->ConstructRectangularMesh(slab_nodes_x, slab_nodes_y);
                            
                            // place at origin
//                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->Translate(-(double)slab_nodes_x/2.0,
//                                             -(double)slab_nodes_y/2.0,
//                                             -(double)slab_nodes_z/2.0);
                                             
                            // scale
                            double mesh_scale_factor = inter_node_space;
                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->Scale(mesh_scale_factor, mesh_scale_factor, mesh_scale_factor);
                            break;
                        }
                        case 3:
                        {
                            c_vector<double, 3> slab_dimensions; //cm                    
                            HeartConfig::Instance()->GetSlabDimensions(slab_dimensions);
                            double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();
                            
                            mpMesh = new TetrahedralMesh<SPACE_DIM, SPACE_DIM>();
                            
                            unsigned slab_nodes_x = (unsigned)round(slab_dimensions[0]/inter_node_space);
                            unsigned slab_nodes_y = (unsigned)round(slab_dimensions[1]/inter_node_space);
                            unsigned slab_nodes_z = (unsigned)round(slab_dimensions[2]/inter_node_space);
                        
                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->ConstructCuboid(slab_nodes_x,
                                                   slab_nodes_y,
                                                   slab_nodes_z,
                                                   true);
                            // place at origin
//                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->Translate(-(double)slab_nodes_x/2.0,
//                                             -(double)slab_nodes_y/2.0,
//                                             -(double)slab_nodes_z/2.0);
                                             
                            // scale
                            double mesh_scale_factor = inter_node_space;
                            static_cast<TetrahedralMesh<SPACE_DIM, SPACE_DIM>*>(mpMesh)->Scale(mesh_scale_factor, mesh_scale_factor, mesh_scale_factor);
                            break;
                        }
                        default:
                            NEVER_REACHED;
                    }
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
            EXCEPTION(std::string("No mesh given: define it in XML parameters file or call SetMesh()\n") + e.GetMessage());
        }
    }
    mpCellFactory->SetMesh( mpMesh );

    if (mNodesPerProcessorFilename != "")
    {
        mpMesh->ReadNodesPerProcessorFile(mNodesPerProcessorFilename);
    }
    ///\todo Should this method be rolled into the Solve() method or the PreSolveChecks()?
    delete mpCardiacPde; // In case we're called twice
    mpCardiacPde = CreateCardiacPde();
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::SetNodesPerProcessorFilename(const std::string& filename)
{
    mNodesPerProcessorFilename = filename;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::SetBoundaryConditionsContainer(BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM> *bcc)
{
    this->mpBoundaryConditionsContainer = bcc;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::PreSolveChecks()
{
    if ( mpCardiacPde == NULL ) // if pde is NULL, Initialise() probably hasn't been called
    {
        EXCEPTION("Pde is null, Initialise() probably hasn't been called");
    }
    if ( HeartConfig::Instance()->GetSimulationDuration() <= 0.0)
    {
        EXCEPTION("End time should be greater than 0");
    }
    if (mPrintOutput==true)
    {
        if( (mOutputDirectory=="") || (mOutputFilenamePrefix==""))
        {
            EXCEPTION("Either explicitly specify not to print output (call PrintOutput(false)) or specify the output directory and filename prefix");
        }
    }
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::CreateInitialCondition()
{
    //if (DistributedVector::GetProblemSize()==0)
    //{
    //    DistributedVector::SetProblemSize(mpMesh->GetNumNodes());
    //}
    Vec initial_condition = DistributedVector::CreateVec(PROBLEM_DIM);
    DistributedVector ic(initial_condition);
    std::vector<DistributedVector::Stripe> stripe;
    stripe.reserve(PROBLEM_DIM);

    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        stripe.push_back(DistributedVector::Stripe(ic, i));
    }

    for (DistributedVector::Iterator index = DistributedVector::Begin();
         index != DistributedVector::End();
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

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::ConvertOutputToMeshalyzerFormat(bool call)
{
    mCallChaste2Meshalyzer=call;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::SetMesh(AbstractMesh<SPACE_DIM,SPACE_DIM>* pMesh)
{
    // If this fails the mesh has already been set. We assert rather throw an exception
    // to avoid a memory leak when checking it throws correctly
    assert(mpMesh==NULL);
    assert(pMesh!=NULL);
    mAllocatedMemoryForMesh = false;
    mpMesh = pMesh;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::PrintOutput(bool printOutput)
{
    mPrintOutput = printOutput;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::SetWriteInfo(bool writeInfo)
{
    mWriteInfo = writeInfo;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::GetSolution()
{
    return mSolution;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractMesh<SPACE_DIM,SPACE_DIM> & AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::rGetMesh()
{
    assert (mpMesh);        
    return *mpMesh;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacPde<SPACE_DIM>* AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::GetPde()
{
    return mpCardiacPde;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::Solve()
{
    PreSolveChecks();

    if(mpBoundaryConditionsContainer == NULL) // the user didnt supply a bcc
    {
        // set up the default bcc
        mpDefaultBoundaryConditionsContainer = new BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>;
        for (unsigned problem_index=0; problem_index<PROBLEM_DIM; problem_index++)
        {
            mpDefaultBoundaryConditionsContainer->DefineZeroNeumannOnMeshBoundary(mpMesh, problem_index);
        }
        mpBoundaryConditionsContainer = mpDefaultBoundaryConditionsContainer;
    }

    mpAssembler = CreateAssembler(); // passes mpBoundaryConditionsContainer to assember
    Vec initial_condition = CreateInitialCondition();

    TimeStepper stepper(0.0, HeartConfig::Instance()->GetSimulationDuration(), 
                        HeartConfig::Instance()->GetPrintingTimeStep());

    std::string progress_reporter_dir;

    if (mPrintOutput)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
        InitialiseWriter();
        WriteOneStep(stepper.GetTime(), initial_condition);
        HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
        
        progress_reporter_dir = mOutputDirectory;
    }
    else
    {
        progress_reporter_dir = ""; // progress printed to CHASTE_TEST_OUTPUT
    }

    // create a progress reporter so users can track how much has gone and 
    // estimate how much time is left. (Note this has to be done after the
    // InitialiseWriter above (if mPrintOutput==true)
    ProgressReporter progress_reporter(progress_reporter_dir, 0.0, HeartConfig::Instance()->GetSimulationDuration());
    progress_reporter.Update(0);

    // If we have already run a simulation, free the old solution vec
    if (mSolution)
    {
        VecDestroy(mSolution);
    }

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
            //VecDestroy(initial_condition);
            
            PetscTools::ReplicateException(true);
            // Re-throw
            HeartEventHandler::Reset();//EndEvent(HeartEventHandler::EVERYTHING);
            
            CloseFilesAndPostProcess();
            throw e;
        }
        PetscTools::ReplicateException(false);

        // Free old initial condition
        VecDestroy(initial_condition);

        // Initial condition for next loop is current solution
        initial_condition = mSolution;

        // update the current time
        stepper.AdvanceOneTimeStep();

        if (mPrintOutput)
        {
            // print out details at current time if asked for
            if (mWriteInfo)
            {
                WriteInfo(stepper.GetTime());
            }

            // Writing data out to the file <mOutputFilenamePrefix>.dat
            HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
            mpWriter->AdvanceAlongUnlimitedDimension(); //creates a new file
            WriteOneStep(stepper.GetTime(), mSolution);
            HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
        }
        
        progress_reporter.Update(stepper.GetTime());
        
        OnEndOfTimestep(stepper.GetTime());
    }

    // Free assembler
    delete mpAssembler;

    // close the file that stores voltage values
    progress_reporter.PrintFinalising();
    CloseFilesAndPostProcess();
    HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::CloseFilesAndPostProcess()
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
    if (mCallChaste2Meshalyzer && mNodesToOutput.empty()) 
    {
        //Convert simulation data to Meshalyzer format
        std::string output_directory =  mOutputDirectory + "/output";
        Hdf5ToMeshalyzerConverter converter(mOutputDirectory, output_directory, mOutputFilenamePrefix);
        
        //Write mesh in a suitable form for meshalyzer
        if (PetscTools::AmMaster())
        {
            //Write the mesh
            MeshalyzerMeshWriter<SPACE_DIM,SPACE_DIM> mesh_writer(output_directory, mOutputFilenamePrefix+"_mesh", false);

            try
            {
                // If this mesh object has been constructed from a mesh reader we can get reference to it            
                TrianglesMeshReader<SPACE_DIM,SPACE_DIM> mesh_reader(mpMesh->GetMeshFileBaseName());
                mesh_writer.WriteFilesUsingMeshReader(mesh_reader, mpMesh->rGetNodePermutation());
            }
            catch(Exception& e)
            {
                //If there isn't a MeshReader available we will use the data contained in the actual mesh object.                
                ///\todo: WriteFilesUsingMesh cannot handle ParallelTetrahedralMesh objects. Abort if so.                            
                mesh_writer.WriteFilesUsingMesh(*mpMesh);
            }
            
            //Write the parameters out
            HeartConfig::Instance()->Write(output_directory, mOutputFilenamePrefix+"_parameters.xml");
        }
    }
    HeartEventHandler::EndEvent(HeartEventHandler::USER2); //Temporarily using USER2 to instrument post-processing
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::DefineWriterColumns()
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

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::WriteOneStep(double time, Vec voltageVec)
{
    mpWriter->PutUnlimitedVariable(time);

    //DistributedVector::Stripe transmembrane(voltageVec, 0);
    mpWriter->PutVector(mVoltageColumnId, voltageVec);
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::InitialiseWriter()
{
    mpWriter = new Hdf5DataWriter(mOutputDirectory,mOutputFilenamePrefix);
    DefineWriterColumns();
    mpWriter->EndDefineMode();
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::SetOutputNodes(std::vector<unsigned> &nodesToOutput)
{
    mNodesToOutput = nodesToOutput;
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Hdf5DataReader AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::GetDataReader()
{
    if( (mOutputDirectory=="") || (mOutputFilenamePrefix==""))
    {
        EXCEPTION("Data reader invalid as data writer cannot be initialised");
    }
    return Hdf5DataReader(mOutputDirectory, mOutputFilenamePrefix);
}

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<SPACE_DIM,PROBLEM_DIM>::UseMatrixBasedRhsAssembly(bool usematrix)
{
    mUseMatrixBasedRhsAssembly = usematrix;
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

// Monodomain
template class AbstractCardiacProblem<1,1>;
template class AbstractCardiacProblem<2,1>;
template class AbstractCardiacProblem<3,1>;

// Bidomain
template class AbstractCardiacProblem<1,2>;
template class AbstractCardiacProblem<2,2>;
template class AbstractCardiacProblem<3,2>;
