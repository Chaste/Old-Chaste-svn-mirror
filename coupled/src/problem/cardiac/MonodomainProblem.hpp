#ifndef MONODOMAINPROBLEM_HPP_
#define MONODOMAINPROBLEM_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.cpp"
#include "ParallelColumnDataWriter.hpp"
#include "MonodomainPde.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "TimeStepper.hpp"
#include "AbstractCardiacProblem.hpp"


/**
 * Class which specifies and solves a monodomain problem.
 */
template<unsigned SPACE_DIM>
class MonodomainProblem : public AbstractCardiacProblem<SPACE_DIM>
{
private:
    MonodomainPde<SPACE_DIM> *mpMonodomainPde;    

public:

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    MonodomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : AbstractCardiacProblem<SPACE_DIM>(pCellFactory),
            mpMonodomainPde(NULL)
    {
    }
    
    /**
     * Destructor
     */
    ~MonodomainProblem()
    {
        delete mpMonodomainPde;
    }
    
    /** Initialise the system. Must be called before Solve() */
    void Initialise()
    {
        AbstractCardiacProblem<SPACE_DIM>::Initialise(mpMonodomainPde);
        mpMonodomainPde = new MonodomainPde<SPACE_DIM>(this->mpCellFactory);
    }
    
    /**
     * Solve the problem
     */
    void Solve()
    {

        AbstractCardiacProblem<SPACE_DIM>::PreSolveChecks(mpMonodomainPde);
        // Assembler
        MonodomainDg0Assembler<SPACE_DIM,SPACE_DIM> monodomain_assembler(&this->mMesh, mpMonodomainPde);
        
        // initial condition;     
        Vec initial_condition = AbstractCardiacProblem<SPACE_DIM>::CreateInitialCondition(mpMonodomainPde, 1);

        
        //  Write data to a file <this->mOutputFilenamePrefix>_xx.dat, 'xx' refers to
        //  'xx'th time step using ColumnDataWriter
        ParallelColumnDataWriter *p_test_writer = NULL;
        
        unsigned time_var_id = 0;
        unsigned voltage_var_id = 0;
        bool write_files = false;
        if (this->mOutputFilenamePrefix.length() > 0)
        {
            write_files = true;
            
            p_test_writer = new ParallelColumnDataWriter(this->mOutputDirectory,this->mOutputFilenamePrefix);
            
            p_test_writer->DefineFixedDimension("Node", "dimensionless", this->mMesh.GetNumNodes() );
            time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
            
            voltage_var_id = p_test_writer->DefineVariable("V","mV");
            p_test_writer->EndDefineMode();
        }
            
        TimeStepper stepper(this->mStartTime, this->mEndTime, this->mPrintingTimeStep);

        if (write_files)
        {
            p_test_writer->PutVariable(time_var_id, stepper.GetTime());
            p_test_writer->PutVector(voltage_var_id, initial_condition);
                              
        }
        // check the printing time step is a multiple of the pde timestep.
        assert (this->mPrintingTimeStep >= this->mPdeTimeStep);
        assert(  fabs( (this->mPrintingTimeStep/this->mPdeTimeStep)
                           -round(this->mPrintingTimeStep/this->mPdeTimeStep) ) < 1e-10 );
                       
                 
        while ( !stepper.IsTimeAtEnd() )
        {
            // solve from now up to the next printing time
            monodomain_assembler.SetTimes(stepper.GetTime(), stepper.GetNextTime(), this->mPdeTimeStep);
            monodomain_assembler.SetInitialCondition( initial_condition );
            
            try
            {
                this->mVoltage = monodomain_assembler.Solve();
            }
            catch (Exception &e)
            {
                if (write_files)
                {
                    p_test_writer->Close();
                    delete p_test_writer;
                }
                throw e;
            }
            
            // Free old initial condition
            VecDestroy(initial_condition);
            
            // Initial condition for next loop is current solution
            initial_condition = this->mVoltage;
            
            // advance to next time
            stepper.AdvanceOneTimeStep();
            
            // print out details at current time if asked for
            if (this->mWriteInfo)
            {
                WriteInfo(stepper.GetTime());
            }
            
            // Writing data out to the file <this->mOutputFilenamePrefix>.dat
            if (write_files)
            {
                p_test_writer->AdvanceAlongUnlimitedDimension(); // creates a new file
                p_test_writer->PutVariable(time_var_id, stepper.GetTime());
                p_test_writer->PutVector(voltage_var_id, this->mVoltage);
            }
        }
        
        // close the file that stores voltage values
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
                                  + this->mMeshFilename + " "         // arg 2 is mesh prefix, path relative to
                                  // the main chaste directory.
                                  + this->mOutputDirectory + "/"
                                  + this->mOutputFilenamePrefix + " " // arg 3 is the results folder and prefix,
                                  // relative to the testoutput folder.
                                  + "last_simulation";          // arg 4 is the output prefix, relative to
            // anim folder.
            
            system(chaste_2_meshalyzer.c_str());
        }
        
        
    }
    
    
    
    void SetStartTime(const double &rStartTime)
    {
        this->mStartTime = rStartTime;
    }
    
    void SetEndTime(const double &rEndTime)
    {
        this->mEndTime = rEndTime;
    }
    
    /**
     * Set the PDE time step.
     * \todo SetPdeAndPrintingTimeStep
     * Note that the printing time step should also set with this call.
     * Move assertion
     */
    
    void SetPdeTimeStep(double pdeTimeStep)
    {
        if (pdeTimeStep <= 0)
        {
            EXCEPTION("Pde time step should be positive");
        }
        this->mPdeTimeStep = pdeTimeStep;
    }
    
    /** Set the times to print output. The printing time step must be
     *  a multiple of the pde timestep. If SetPdeTimeStep is used it should be called
     * before SetPrintingTimeStep.
     */
    void SetPrintingTimeStep(double printingTimeStep)
    {
        if (printingTimeStep <= 0.0)
        {
            EXCEPTION("Printing time step should be positive");
        }
        this->mPrintingTimeStep = printingTimeStep;
    }
    
    /** Set the simulation to print every n timesteps. Only set this
     *  AFTER setting the pde timestep
     */
    void PrintEveryNthTimeStep(unsigned n)
    {
        this->mPrintingTimeStep = n*this->mPdeTimeStep;
    }
    
    
    double GetPdeTimeStep()
    {
        return this->mPdeTimeStep;
    }
    
    void SetMeshFilename(const std::string &rMeshFilename)
    {
        if ( this->mMeshFilename!="" )
        {
            EXCEPTION("Mesh filename was already set");
        }
        if ( rMeshFilename=="" )
        {
            EXCEPTION("Mesh filename was passed in empty");
        }
        this->mMeshFilename = rMeshFilename;
        
        TrianglesMeshReader<SPACE_DIM, SPACE_DIM> mesh_reader(this->mMeshFilename);
        this->mMesh.ConstructFromMeshReader(mesh_reader);
    }
    
    void SetOutputDirectory(const std::string &rOutputDirectory)
    {
        this->mOutputDirectory = rOutputDirectory;
    }
    
    void SetOutputFilenamePrefix(const std::string &rOutputFilenamePrefix)
    {
        this->mOutputFilenamePrefix = rOutputFilenamePrefix;
    }
    
    /**
     * Get the final solution vector. This vector is distributed over all processes.
     */
    
    Vec GetVoltage()
    {
        //Use with caution since we don't want to alter the state of the PETSc vector
        return this->mVoltage;
    }
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh()
    {
        return this->mMesh;
    }
    
    MonodomainPde<SPACE_DIM> * GetMonodomainPde()
    {
        return mpMonodomainPde;
    }
    
    /**
     *  Set info to be printed during computation. 
     */
    void SetWriteInfo(bool writeInfo = true)
    {
        this->mWriteInfo = writeInfo;
    }
    
    
    /**
     *  Print out time and max/min voltage values at current time.
     
     */
    void WriteInfo(double time)
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
        
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(this->mVoltage);
        
        double v_max = -1e5, v_min = 1e5;
        for (unsigned i=0; i<this->mMesh.GetNumNodes(); i++)
        {
            if ( voltage_replicated[i] > v_max)
            {
                v_max = voltage_replicated[i];
            }
            if ( voltage_replicated[i] < v_min)
            {
                v_min = voltage_replicated[i];
            }
        }
        
        std::cout << " max/min V = "
        <<   v_max << " "
        <<   v_min << "\n" << std::flush;
    }
};

#endif /*MONODOMAINPROBLEM_HPP_*/
