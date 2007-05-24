#ifndef BIDOMAINPROBLEM_HPP_
#define BIDOMAINPROBLEM_HPP_


#include "ConformingTetrahedralMesh.cpp"
#include "BidomainDg0Assembler.hpp"
#include "TrianglesMeshReader.cpp"
#include "ParallelColumnDataWriter.hpp"
#include "BidomainPde.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "TimeStepper.hpp"
#include "AbstractCardiacProblem.hpp"


/**
 * Class which specifies and solves a bidomain problem.
 *
 * The solution vector is of the form:
 * (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N),
 * where V_j is the voltage at node j and phi_j is the
 * extracellular potential at node j.
 */
template<unsigned SPACE_DIM>
class BidomainProblem : public AbstractCardiacProblem<SPACE_DIM>
{
private:    
    BidomainPde<SPACE_DIM>* mpBidomainPde;

            
    std::vector<unsigned> mFixedExtracellularPotentialNodes; /** nodes at which the extracellular voltage is fixed to zero (replicated) */
    
    double mLinearSolverRelativeTolerance;
public:

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    BidomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
            : AbstractCardiacProblem<SPACE_DIM>(pCellFactory),
            mpBidomainPde(NULL)
    {        
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
        AbstractCardiacProblem<SPACE_DIM>::Initialise(mpBidomainPde);
        mpBidomainPde = new BidomainPde<SPACE_DIM>( this->mpCellFactory );
    }
    
    
    /**
     * Solve the problem
     */
    void Solve()
    {
        AbstractCardiacProblem<SPACE_DIM>::PreSolveChecks(mpBidomainPde);
        
        // Assembler
        BidomainDg0Assembler<SPACE_DIM,SPACE_DIM> bidomain_assembler(
            &this->mMesh, mpBidomainPde,
            2, mLinearSolverRelativeTolerance);
            
        if (mFixedExtracellularPotentialNodes.size()>0)
        {
            bidomain_assembler.SetFixedExtracellularPotentialNodes(mFixedExtracellularPotentialNodes);
        }

        Vec initial_condition = AbstractCardiacProblem<SPACE_DIM>::CreateInitialCondition(mpBidomainPde, 2);

        ParallelColumnDataWriter *p_test_writer = NULL;
        unsigned time_var_id = 0;
        unsigned voltage_var_id = 0;
        bool write_files = false;

        TimeStepper stepper(this->mStartTime, this->mEndTime, this->mPrintingTimeStep);

        if (this->mPrintOutput)
        {
            if (this->mOutputFilenamePrefix.length() > 0)
            {
                write_files = true;
                
                p_test_writer = new ParallelColumnDataWriter(this->mOutputDirectory,this->mOutputFilenamePrefix);
                
                p_test_writer->DefineFixedDimension("Node", "dimensionless", 2*this->mMesh.GetNumNodes() );
                time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
                
                voltage_var_id = p_test_writer->DefineVariable("Vm_And_Phi_e","mV");
                p_test_writer->EndDefineMode();
            }
            
            if (write_files)
            {
                p_test_writer->PutVariable(time_var_id, stepper.GetTime());
                p_test_writer->PutVector(voltage_var_id, initial_condition);
            }
        }
        
        while ( !stepper.IsTimeAtEnd() )
        {
            // solve from now up to the next printing time
            bidomain_assembler.SetTimes(stepper.GetTime(), stepper.GetNextTime(), this->mPdeTimeStep);
            bidomain_assembler.SetInitialCondition( initial_condition );

            try
            {
                this->mVoltage = bidomain_assembler.Solve();
            }
            //Ill-conditioned solutions are covered in Monodomain problem
            //(and possibly in Nightly/Weekly) so we don't insist on it
            //in the coverage test.
            #define COVERAGE_IGNORE
            catch (Exception &e)
            {
                if (this->mPrintOutput)
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
            initial_condition = this->mVoltage;
            
            // update the current time
            stepper.AdvanceOneTimeStep();
            
            if (this->mPrintOutput)
            {
                // print out details at current time if asked for
                if (this->mWriteInfo)
                {
                    WriteInfo(stepper.GetTime());
                }
                
                // Writing data out to the file <this->mOutputFilenamePrefix>.dat
                if (write_files)
                {
                    p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
                    p_test_writer->PutVariable(time_var_id, stepper.GetTime());
                    p_test_writer->PutVector(voltage_var_id, this->mVoltage);
                }
            }
        }

        // close the file that stores voltage values
        if (this->mPrintOutput)
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
        this->mWriteInfo = writeInfo;
    }
    
    
    /**
     *  Print out time and max/min voltage/phi_e values at current time.
     */
    void WriteInfo(double time)
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
        
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(this->mVoltage);
        
        double v_max = -1e5, v_min = 1e5, phi_max = -1e5, phi_min = 1e5;
        for (unsigned i=0; i<this->mMesh.GetNumNodes(); i++)
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
