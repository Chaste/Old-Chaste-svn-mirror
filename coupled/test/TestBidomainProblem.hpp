#ifndef TESTBIDOMAINPROBLEM_HPP_
#define TESTBIDOMAINPROBLEM_HPP_


// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <vector>
//#include <iostream>
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"


class PointStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PointStimulusCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(-600, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus, mTimeStep);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
        }
    }
        
    ~PointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};  


// new cell factory inheriting from PointStimulusCellFactory but adding
// zero extracellular stimuli.
class BidomainPointStimulusCellFactory : public PointStimulusCellFactory
{
public:
    // constructor calls base constructor
    BidomainPointStimulusCellFactory() : PointStimulusCellFactory() {}
        
    //\todo: come up with a better way of adding extra_cell stimuli than using finalise()
    void FinaliseCellCreation(std::vector< AbstractCardiacCell* >* pCellsDistributed, int lo, int hi)
    {
        for(int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            (*pCellsDistributed)[local_index]->SetExtracellularStimulusFunction(mpZeroStimulus);
        }
    }
};    
  
    

class TestBidomainDg0Assembler : public CxxTest::TestSuite 
{
public:
    void TestBidomainDg01D()
    {
        BidomainPointStimulusCellFactory bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(1);   // 1 ms
        bidomain_problem.SetOutputDirectory("testoutput/bidomainDg01d");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_1d");

        bidomain_problem.Initialise();
        try
        {
            bidomain_problem.Solve();
        }
        catch(Exception e)
        {
            TS_TRACE(e.GetMessage());
        }
            
        double* p_voltage_array;
        int lo, hi;
        bidomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 

        //lo and hi not set correctly in problem
//        for (int global_index=lo; global_index<hi; global_index++)
        for(int global_index=0; global_index<11; global_index++)
        {
            int local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;   // mV 
            double Ek    = -77.0;   // mV

            TS_ASSERT_LESS_THAN_EQUALS( p_voltage_array[local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage_array[local_index] + (Ek-30), 0);

            std::vector<double> odeVars = bidomain_problem.GetBidomainPde()->GetCardiacCell(global_index)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1
                if ((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);
                }
            }
        }
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        if (num_procs == 1)
        {
            /// \todo This is work in progress
            //TS_ASSERT_DELTA(p_voltage_array[0], 28.2462, 1e-3);
        }

        bidomain_problem.RestoreVoltageArray(&p_voltage_array);       
    } 
    
    
    // unfinished test - compare bidomain and monodomain
    void TestCompareBidomainProblemWithMonodomain()
    {
        PointStimulusCellFactory cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        monodomain_problem.SetEndTime(1);   // 1 ms
        monodomain_problem.SetOutputDirectory("testoutput/Monodomain1d");
        monodomain_problem.SetOutputFilenamePrefix("monodomain1d");
        
        monodomain_problem.Initialise();
        monodomain_problem.Solve();


        BidomainPointStimulusCellFactory bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        bidomain_problem.SetEndTime(1);   // 1 ms
        bidomain_problem.SetOutputDirectory("testoutput/Bidomain1d");
        bidomain_problem.SetOutputFilenamePrefix("bidomain1d");

        bidomain_problem.Initialise();
        bidomain_problem.Solve();
        
        
        double* p_mono_voltage_array;
        double* p_bi_voltage_array;
        int lo, hi;

        bidomain_problem.GetVoltageArray(&p_bi_voltage_array, lo, hi); 
        monodomain_problem.GetVoltageArray(&p_mono_voltage_array, lo, hi); 
        
        for(int global_index=lo; global_index<hi; global_index++)
        {
            //int local_index = global_index - lo;
            
            //std::cout << p_mono_voltage_array[local_index] << " " << p_bi_voltage_array[local_index] << "\n";
        } 
    }
};

#endif /*TESTBIDOMAINPROBLEM_HPP_*/
