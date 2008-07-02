/*

Copyright (C) University of Oxford, 2008

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


#ifndef ABSTRACTCARDIACPROBLEM_HPP_
#define ABSTRACTCARDIACPROBLEM_HPP_

#include "ConformingTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "Hdf5DataWriter.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "TimeStepper.hpp"
#include "DistributedVector.hpp"
#include "AbstractCardiacPde.hpp"
#include "AbstractDynamicAssemblerMixin.hpp"
#include "EventHandler.hpp"
#include "PetscTools.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"

#include "OrthotropicConductivityTensors.hpp"
#include "AxisymmetricConductivityTensors.hpp"

template<unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractCardiacProblem
{
protected:
    std::string mMeshFilename;
    bool mAllocatedMemoryForMesh;
    std::string mNodesPerProcessorFilename;

    bool mWriteInfo;
    bool mPrintOutput;
    bool mCallChaste2Meshalyzer;
    std::vector<unsigned> mNodesToOutput;

    AbstractCardiacPde<SPACE_DIM>* mpCardiacPde;

    /** data is not written if output directory or output file prefix are not set*/
    std::string  mOutputDirectory, mOutputFilenamePrefix;

protected:
    BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* mpBoundaryConditionsContainer;
    AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* mpAssembler;

    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;

    Vec mVoltage; // Current solution
//    double mLinearSolverTolerance;
//    bool mUseLinearSolverAbsoluteTolerance;


    /**
     * Subclasses must override this method to create a PDE object of the appropriate type.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual AbstractCardiacPde<SPACE_DIM>* CreateCardiacPde() =0;

    /**
     * Subclasses must override this method to create a suitable assembler object.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual AbstractDynamicAssemblerMixin<SPACE_DIM, SPACE_DIM, PROBLEM_DIM>* CreateAssembler() =0;

    unsigned mVoltageColumnId;
    unsigned mTimeColumnId;
    unsigned mNodeColumnId;

public:
    // This (and things in MonodomainProblem) being public are hacks for
    // AbstractCardiacElectroMechanicsWriter to work.
    // TODO AbstractCardiacElectroMechanicsWriter should be a friend, but not sure
    // how to get friends to work when both friends are templated and abstract.
    Hdf5DataWriter* mpWriter;

public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     * create cells.
     */
    AbstractCardiacProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : mMeshFilename(""),     // i.e. undefined
              mNodesPerProcessorFilename(""),     // i.e. undefined
              mOutputDirectory(""),  // i.e. undefined
              mOutputFilenamePrefix(""),   // i.e. undefined
              mpBoundaryConditionsContainer(NULL),
              mpCellFactory(pCellFactory),
              mpMesh(NULL),
              mpWriter(NULL)
    {
        mWriteInfo = false;
        mPrintOutput = true;
        mCallChaste2Meshalyzer = false;
        mpCardiacPde = NULL;
        mpAssembler = NULL;
        mVoltage = NULL;
//        mLinearSolverTolerance=1e-6;
//        mUseLinearSolverAbsoluteTolerance = false;
        mAllocatedMemoryForMesh = false;
        assert(mNodesToOutput.empty());

        EventHandler::BeginEvent(EVERYTHING);
    }

    virtual ~AbstractCardiacProblem()
    {
        delete mpCardiacPde;
        if (mVoltage)
        {
            VecDestroy(mVoltage);
        }

        if(mAllocatedMemoryForMesh)
        {
            delete mpMesh;
        }
    };

    /*
     *  Initialise the system. Must be called before Solve()
     */
    void Initialise()
    {
        if (mpMesh==NULL)
        {
            EXCEPTION("SetMesh() or SetMeshFilename() was not set");
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

    void SetNodesPerProcessorFilename(const std::string& filename)
    {
        mNodesPerProcessorFilename = filename;
    }

    void SetBoundaryConditionsContainer(BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM> *bcc)
    {
        this->mpBoundaryConditionsContainer = bcc;
    }

    void PreSolveChecks()
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

    // Perhaps this should be a method of AbstractCardiacPde??)

    Vec CreateInitialCondition()
    {
        //if (DistributedVector::GetProblemSize()==0)
        //{
        //    DistributedVector::SetProblemSize(mpMesh->GetNumNodes());
        //}
        Vec initial_condition=DistributedVector::CreateVec(PROBLEM_DIM);
        DistributedVector ic(initial_condition);
        std::vector< DistributedVector::Stripe > stripe;
        stripe.reserve(PROBLEM_DIM);

        for (unsigned i=0; i<PROBLEM_DIM; i++)
        {
            stripe.push_back(DistributedVector::Stripe(ic, i));
        }

        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
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

    /** Set whether to call the Chaste2Meshalyzer script.
     * This script gets everything ready to visualize the results with meshalyser
     * and is useful in testing. By default the script is called.
     * In performance testing for example it desirable to disable the script.
     */
    void ConvertOutputToMeshalyzerFormat(bool call = true)
    {
        mCallChaste2Meshalyzer=call;
    }

    void SetMeshFilename(const std::string &rMeshFilename)
    {
        // If this fails the mesh has already been set. We assert rather throw an exception
        // to avoid a memory leak when checking it throws correctly
        assert(mpMesh==NULL);

        if ( rMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was passed in empty");
        }

        mMeshFilename = rMeshFilename;

        TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(mMeshFilename);
        mpMesh = new ConformingTetrahedralMesh<SPACE_DIM, SPACE_DIM>();
        mAllocatedMemoryForMesh = true;

        EventHandler::BeginEvent(READ_MESH);
        mpMesh->ConstructFromMeshReader(mesh_reader);
        EventHandler::EndEvent(READ_MESH);
    }

    void SetMesh(ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
    {
        // If this fails the mesh has already been set. We assert rather throw an exception
        // to avoid a memory leak when checking it throws correctly
        assert(mpMesh==NULL);
        mAllocatedMemoryForMesh = false;
        assert(pMesh!=NULL);
        mpMesh = pMesh;
    }

    void SetOutputDirectory(const std::string& rOutputDirectory)
    {
        mOutputDirectory = rOutputDirectory;
    }

    void SetOutputFilenamePrefix(const std::string& rOutputFilenamePrefix)
    {
        mOutputFilenamePrefix = rOutputFilenamePrefix;
    }

    /**
     *  Set whether the simulation will generate results files.
     */
    void PrintOutput(bool rPrintOutput)
    {
        mPrintOutput = rPrintOutput;
    }

    /**
     *  Set whether extra info will be written to stdout during computation.
     */
    void SetWriteInfo(bool writeInfo = true)
    {
        mWriteInfo = writeInfo;
    }

    /**
     * Get the final solution vector. This vector is distributed over all processes.
     *
     * In case of Bidomain, this is of length 2*numNodes, and of the form
     *  (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N).
     *  where V_j is the voltage at node j and phi_j is the
     *  extracellular potential at node j.
     *
     *  Use with caution since we don't want to alter the state of the PETSc vector.
     */
    Vec GetVoltage()
    {
        return mVoltage;
    }

    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh()
    {
        return *mpMesh;
    }

    AbstractCardiacPde<SPACE_DIM>* GetPde()
    {
        return mpCardiacPde;
    }

    void Solve()
    {
        PreSolveChecks();

        // set default bcc if required
        BoundaryConditionsContainer<SPACE_DIM, SPACE_DIM, PROBLEM_DIM> default_bcc;
        if(mpBoundaryConditionsContainer == NULL) // the user didnt supply a bcc
        {
            for (unsigned problem_index=0; problem_index<PROBLEM_DIM; problem_index++)
            {
                default_bcc.DefineZeroNeumannOnMeshBoundary(mpMesh, problem_index);
            }
            mpBoundaryConditionsContainer = &default_bcc;
        }

        mpAssembler = CreateAssembler(); // passes mpBoundaryConditionsContainer to assember
        Vec initial_condition = CreateInitialCondition();

        TimeStepper stepper(0.0, HeartConfig::Instance()->GetSimulationDuration(), 
                            HeartConfig::Instance()->GetPrintingTimeStep());

        if (mPrintOutput)
        {
            EventHandler::BeginEvent(WRITE_OUTPUT);
            InitialiseWriter();
            WriteOneStep(stepper.GetTime(), initial_condition);
            EventHandler::EndEvent(WRITE_OUTPUT);
        }

        // If we have already run a simulation, free the old solution vec
        if (mVoltage)
        {
            VecDestroy(mVoltage);
        }

        while ( !stepper.IsTimeAtEnd() )
        {
            // solve from now up to the next printing time
            mpAssembler->SetTimes(stepper.GetTime(), stepper.GetNextTime(), HeartConfig::Instance()->GetPdeTimeStep());
            mpAssembler->SetInitialCondition( initial_condition );

            try
            {
                mVoltage = mpAssembler->Solve();
            }
            catch (Exception &e)
            {
                // Free memory.
                delete mpAssembler;
                //VecDestroy(initial_condition);
                // Close files
                if (mPrintOutput)
                {
                    mpWriter->Close();
                    delete mpWriter;
                }
                PetscTools::ReplicateException(true);
                // Re-throw
                //EventHandler::EndEvent(EVERYTHING);
                throw e;
            }
            PetscTools::ReplicateException(false);

            // Free old initial condition
            VecDestroy(initial_condition);

            // Initial condition for next loop is current solution
            initial_condition = mVoltage;

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
                EventHandler::BeginEvent(WRITE_OUTPUT);
                mpWriter->AdvanceAlongUnlimitedDimension(); //creates a new file
                WriteOneStep(stepper.GetTime(), mVoltage);
                EventHandler::EndEvent(WRITE_OUTPUT);
            }
        }

        // Free assembler
        delete mpAssembler;

        // close the file that stores voltage values
        if (mPrintOutput)
        { 
            mpWriter->Close();
            delete mpWriter;

            bool am_master = PetscTools::AmMaster();

            // Only if results files were written and we are outputting all nodes
            if (mCallChaste2Meshalyzer && mNodesToOutput.empty()) 
            {
                
                if (am_master)
                {
                    // call shell script which converts the data to meshalyzer format
                    std::string chaste_2_meshalyzer;
                    std::stringstream space_dim;
                    space_dim << SPACE_DIM;
                    chaste_2_meshalyzer = "anim/chaste2meshalyzer "     // the executable.
                                      + space_dim.str() + " "       // argument 1 is the dimension.
                                      + mMeshFilename + " "         // arg 2 is mesh prefix, path relative to Chaste directory
                                      + "last_simulation";          // arg 3 is the output prefix, relative to
                                                                    // anim folder.
                    system(chaste_2_meshalyzer.c_str());
                }
                //Convert simulation data to Meshalyzer format
                Hdf5ToMeshalyzerConverter converter(mOutputDirectory, mOutputFilenamePrefix);
            }
        }
        EventHandler::EndEvent(EVERYTHING);
    }

    virtual void WriteInfo(double time) =0;

    virtual void DefineWriterColumns()
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

    virtual void WriteOneStep(double time, Vec voltageVec)
    {
        mpWriter->PutUnlimitedVariable(time);

        //DistributedVector::Stripe transmembrane(voltageVec, 0);
        mpWriter->PutVector(mVoltageColumnId, voltageVec);
    }

    void InitialiseWriter()
    {
        mpWriter = new Hdf5DataWriter(mOutputDirectory,mOutputFilenamePrefix);
        DefineWriterColumns();
        mpWriter->EndDefineMode();
    }

    void SetOutputNodes(std::vector<unsigned> &nodesToOutput)
    {
        mNodesToOutput = nodesToOutput;
    }
};
#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
