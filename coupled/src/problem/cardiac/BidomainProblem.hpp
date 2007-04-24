#ifndef BIDOMAINPROBLEM_HPP_
#define BIDOMAINPROBLEM_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "BidomainDg0Assembler.hpp"
#include "TrianglesMeshReader.cpp"
#include "ParallelColumnDataWriter.hpp"
#include "BidomainPde.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"

/**
 * Class which specifies and solves a bidomain problem.
 *
 * The solution vector is of the form:
 * (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N),
 * where V_j is the voltage at node j and phi_j is the
 * extracellular potential at node j.
 */
template<unsigned SPACE_DIM>
class BidomainProblem
{
private:
    std::string mMeshFilename;
    
    /**
     *  Start time defaults to 0, pde timestep defaults to 0.01 (ms), the
     *  end time is not defaulted and must be set
     */
    double mStartTime;
    double mEndTime;
    double mPdeTimeStep;
    double mPrintingTimeStep;
    
    bool mWriteInfo;
    bool mPrintOutput;
    
    /** data is not written if output directory or output file prefix are not set*/
    std::string  mOutputDirectory, mOutputFilenamePrefix;
    
    BidomainPde<SPACE_DIM>* mpBidomainPde;
    
    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    
    Vec mVoltage; // Current solution
    unsigned mLo, mHi; // Ownership range for the double size potential vector mVoltage
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mMesh;
    
    std::vector<unsigned> mFixedExtracellularPotentialNodes; /** nodes at which the extracellular voltage is fixed to zero (replicated) */
    
    double mLinearSolverRelativeTolerance;
public:

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    BidomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : mMeshFilename(""),     // i.e. undefined
            mOutputDirectory(""),  // i.e. undefined
            mOutputFilenamePrefix(""),   // i.e. undefined
            mpBidomainPde(NULL),
            mpCellFactory(pCellFactory)
    {
        mStartTime        = 0.0;  // ms
        mPdeTimeStep      = 0.01; // ms
        mEndTime          = -1;   // negative so can check has been set
        mPrintingTimeStep = -1;   // negative so can check has been set
        
        //initialise these to 0
        mLo = 0;
        mHi = 0;
        
        mWriteInfo = false;
        
        mPrintOutput = true;   // We want some output by default
        
        mFixedExtracellularPotentialNodes.resize(0);
        
        mLinearSolverRelativeTolerance=1e-6;
    }
    
    /**
     * Destructor
     */
    ~BidomainProblem()
    {
        if (mpBidomainPde)
        {
            delete mpBidomainPde;
        }
    }
    
    
    /** Initialise the system. Must be called before Solve() */
    void Initialise()
    {
        if ( mMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was not set");
        }
        
        mpCellFactory->SetMesh( &mMesh );
        
        if (mpBidomainPde)
        {
            delete mpBidomainPde;
        }
        mpBidomainPde = new BidomainPde<SPACE_DIM>( mpCellFactory );
    }
    
    
    /**
     * Solve the problem
     */
    void Solve()
    {
        if ( mpBidomainPde == NULL ) // if pde is NULL, Initialise() probably hasn't been called
        {
            EXCEPTION("Bidomain pde is null, Initialise() probably hasn't been called");
        }
        
        if ( mStartTime >= mEndTime )
        {
            EXCEPTION("Start time should be less than end time");
        }
        
        // Assembler
        BidomainDg0Assembler<SPACE_DIM,SPACE_DIM> bidomain_assembler(
            &mMesh, mpBidomainPde,
            2, mLinearSolverRelativeTolerance);
            
        if (mFixedExtracellularPotentialNodes.size()>0)
        {
            bidomain_assembler.SetFixedExtracellularPotentialNodes(mFixedExtracellularPotentialNodes);
        }
        
        
        DistributedVector::SetProblemSize(mMesh.GetNumNodes());
                
        
        Vec initial_condition=DistributedVector::CreateVec(2);
        
        DistributedVector striped_vec(initial_condition);
        DistributedVector::Stripe intra(striped_vec,0);
        DistributedVector::Stripe extra(striped_vec,1);
                
        mLo=DistributedVector::Begin().Global;
        mHi=DistributedVector::End().Global;
        
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            intra[index]= mpBidomainPde->GetCardiacCell(index.Global)->GetVoltage();
            extra[index] =0;
        }        
        
        striped_vec.Restore();

        ParallelColumnDataWriter *p_test_writer = NULL;
        unsigned time_var_id = 0;
        unsigned voltage_var_id = 0;
        bool write_files = false;
        double current_time = mStartTime;
        if (mPrintOutput)
        {
            if (mOutputFilenamePrefix.length() > 0)
            {
                write_files = true;
                
                p_test_writer = new ParallelColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);
                
                p_test_writer->DefineFixedDimension("Node", "dimensionless", 2*mMesh.GetNumNodes() );
                time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
                
                voltage_var_id = p_test_writer->DefineVariable("Vm_And_Phi_e","mV");
                p_test_writer->EndDefineMode();
            }
            
            if (write_files)
            {
                p_test_writer->PutVariable(time_var_id, current_time);
                p_test_writer->PutVector(voltage_var_id, initial_condition);
            }
            
            if ( mPrintingTimeStep < 0) //ie if it hasn't been set
            {
                mPrintingTimeStep = mPdeTimeStep; ///\todo: pick good default
            }
            
            
            // at the moment the tend-tstart must be a multiple of the printing timestep
            // because the 'next_printing_time = mEndTime' line below has been commented
            // out (see comment below). remove this assert when this is sorted out
            assert( fabs(        ( (mEndTime-mStartTime)/mPrintingTimeStep )
                                 - round( (mEndTime-mStartTime)/mPrintingTimeStep ) ) < 1e-10 );
                                 
                                 
            // check the printing time step is a multiple of the pde timestep.
            assert( fabs(        (mPrintingTimeStep/mPdeTimeStep)
                                 -round(mPrintingTimeStep/mPdeTimeStep) ) < 1e-10 );
        }
        
        double next_printing_time;
        
        while ( current_time < mEndTime )
        {
            if (mPrintOutput)
            {
                // compute the next printing time
                next_printing_time = current_time + mPrintingTimeStep;
                if (next_printing_time > mEndTime)
                {
                    // this line is needed but it's use leads to assertions tripping due to
                    // floating point errors.
                    ///\todo: sort this out! then change TestPrintsOnlyAtRequestedTimes in
                    // TestMonodomainDg0Assembler and TestBidomainProblem to have an endtime that
                    // is not a multiple of the printing time...
                    //                next_printing_time = mEndTime;
                }
            }
            else
            {
                next_printing_time = current_time + mPdeTimeStep;
            }
            
            // solve from now up to the next printing time
            bidomain_assembler.SetTimes(current_time, next_printing_time, mPdeTimeStep);
            bidomain_assembler.SetInitialCondition( initial_condition );
            
            try
            {
                mVoltage = bidomain_assembler.Solve();
            }
            //Ill-conditioned solutions are covered in Monodomain problem
            //(and possibly in Nightly/Weekly) so we don't insist on it
            //in the coverage test.
            
#define COVERAGE_IGNORE
            catch (Exception &e)
            {
                if (mPrintOutput)
                {
                    if (write_files)
                    {
                        p_test_writer->Close();
                        delete p_test_writer;
                    }
                }
                
                throw e;
            }
#undef COVERAGE_IGNORE
            
            // Free old initial condition
            VecDestroy(initial_condition);
            
            // Initial condition for next loop is current solution
            initial_condition = mVoltage;
            
            // update the current time
            current_time = next_printing_time;
            
            if (mPrintOutput)
            {
                // print out details at current time if asked for
                
                if (mWriteInfo)
                {
                    WriteInfo(current_time);
                }
                
                // Writing data out to the file <mOutputFilenamePrefix>.dat
                if (write_files)
                {
                    p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
                    p_test_writer->PutVariable(time_var_id, current_time);
                    p_test_writer->PutVector(voltage_var_id, mVoltage);
                }
            }
        }
        
        // close the file that stores voltage values
        
        if (mPrintOutput)
        {
            if (write_files)
            {
                p_test_writer->Close();
                delete p_test_writer;
            }
            
            PetscInt my_rank;
            MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
            if ((my_rank==0) && (write_files)) // ie only if master process and results files were written
            {
                // call shell script which converts the data to meshalyzer format
                std::string chaste_2_meshalyzer;
                std::stringstream space_dim;
                space_dim << SPACE_DIM;
                chaste_2_meshalyzer = "anim/chaste2meshalyzer "         // the executable.
                                      + space_dim.str() + " "       // argument 1 is the dimension.
                                      + mMeshFilename + " "         // arg 2 is mesh prefix, path relative to
                                      // the main chaste directory.
                                      + mOutputDirectory + "/"
                                      + mOutputFilenamePrefix + " " // arg 3 is the results folder and prefix,
                                      // relative to the testoutput folder.
                                      + "last_simulation";          // arg 4 is the output prefix, relative to
                // anim folder.
                
                system(chaste_2_meshalyzer.c_str());
            }
        }
    }
    
    
    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to 
     *  zero. This does not necessarily have to be called. If it is not, phi_e 
     *  is only defined up to a constant.
     * 
     *  @param the nodes to be fixed.
     * 
     *  NOTE: currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> nodes)
    {
        if (nodes.size() == 0)
        {
            EXCEPTION("Number of fixed nodes should be greater than zero");
        }
        
        mFixedExtracellularPotentialNodes.resize(nodes.size());
        for (unsigned i=0; i<nodes.size(); i++)
        {
            // the assembler checks that the nodes[i] is less than
            // the number of nodes in the mesh so this is not done here
            mFixedExtracellularPotentialNodes[i] = nodes[i];
        }
    }
    
    
    
    
    
    void SetStartTime(const double &rStartTime)
    {
        mStartTime = rStartTime;
    }
    
    void SetLinearSolverRelativeTolerance(const double &rRelTol)
    {
        mLinearSolverRelativeTolerance = rRelTol;
    }
    
    void SetEndTime(const double &rEndTime)
    {
        mEndTime = rEndTime;
    }
    
    void SetPdeTimeStep(double pdeTimeStep)
    {
        if (pdeTimeStep <= 0)
        {
            EXCEPTION("Pde time step should be positive");
        }
        mPdeTimeStep = pdeTimeStep;
    }
    
    /**
     *  Set the times to print output. The printing time step must be 
     *  a multiple of the pde timestep 
     */
    void SetPrintingTimeStep(double printingTimeStep)
    {
        if (printingTimeStep <= 0.0)
        {
            EXCEPTION("Printing time step should be positive");
        }
        mPrintingTimeStep = printingTimeStep;
    }
    
    /**
     *  Set the simulation to print every n timesteps. Only set this
     *  AFTER setting the pde timestep
     */
    void PrintEveryNthTimeStep(unsigned n)
    {
        mPrintingTimeStep = n*mPdeTimeStep;
    }
    
    
    /**
     *  Set the simulation to print every n timesteps. Only set this
     *  AFTER setting the pde timestep
     */
    void PrintOutput(const bool& rPrintOutput)
    {
        mPrintOutput = rPrintOutput;
    }
    
    
    double GetPdeTimeStep()
    {
        return mPdeTimeStep;
    }
    
    void SetMeshFilename(const std::string &rMeshFilename)
    {
        if ( mMeshFilename!="" )
        {
            EXCEPTION("Mesh filename was already set");
        }
        if ( rMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was passed in empty");
        }
        
        
        mMeshFilename = rMeshFilename;
        
        TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(mMeshFilename);
        mMesh.ConstructFromMeshReader(mesh_reader);
        
    }
    
    void SetOutputDirectory(const std::string &rOutputDirectory)
    {
        mOutputDirectory = rOutputDirectory;
    }
    
    void SetOutputFilenamePrefix(const std::string &rOutputFilenamePrefix)
    {
        mOutputFilenamePrefix = rOutputFilenamePrefix;
    }
    
    
    /**
     *  Get the final solution vector. This is of length 2*numNodes, and of the form
     *  (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N). 
     *  where V_j is the voltage at node j and phi_j is the
     *  extracellular potential at node j.
     * 
     *  This vector is distributed over all processes.
     */

    
    Vec GetVoltage()
    {
        return mVoltage;
    }
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh()
    {
        return mMesh;
    }
    
    /**
     *  Get the pde. Can only be called after Initialise()
     */
    BidomainPde<SPACE_DIM>* GetBidomainPde()
    {
        assert(mpBidomainPde!=NULL);
        return mpBidomainPde;
    }
    
    /**
     *  Set info to be printed during computation. 
     */
    void SetWriteInfo(bool writeInfo = true)
    {
        mWriteInfo = writeInfo;
    }
    
    
    /**
     *  Print out time and max/min voltage/phi_e values at current time.
     */
    void WriteInfo(double time)
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
        
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(mVoltage);
        
        double v_max = -1e5, v_min = 1e5, phi_max = -1e5, phi_min = 1e5;
        for (unsigned i=0; i<mMesh.GetNumNodes(); i++)
        {
            if ( voltage_replicated[2*i] > v_max)
            {
                v_max = voltage_replicated[2*i];
            }
            if ( voltage_replicated[2*i] < v_min)
            {
                v_min = voltage_replicated[2*i];
            }
            if ( voltage_replicated[2*i+1] > phi_max)
            {
                phi_max = voltage_replicated[2*i+1];
            }
            if ( voltage_replicated[2*i+1] < phi_min)
            {
                phi_min = voltage_replicated[2*i+1];
            }
            
        }
        
        std::cout << " max/min V, phi_e = "
        << v_max << " "
        << v_min << " "
        << phi_max << " "
        << phi_min << "\n" << std::flush;
    }
};


#endif /*BIDOMAINPROBLEM_HPP_*/
