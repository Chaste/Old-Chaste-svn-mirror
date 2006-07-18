#ifndef TESTBIDOMAIN3D_HPP_
#define TESTBIDOMAIN3D_HPP_




#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"


class BidomainFaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    BidomainFaceStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.01)
    {
        mpStimulus = new InitialStimulus(-600.0*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            //std::cout << node+1 << "\n";
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
        }
    }
        
    ~BidomainFaceStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestBidomain3D :  public CxxTest::TestSuite 
{
public:
    
    void TestBidomain3d()
    {
        BidomainFaceStimulusCellFactory bidomain_cell_factory;
        
        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        bidomain_problem.SetEndTime(4);   // ms
        bidomain_problem.SetOutputDirectory("Bidomain3d");
        bidomain_problem.SetOutputFilenamePrefix("bidomain3d");

        bidomain_problem.Initialise();

        bidomain_problem.Solve();
        
        Vec voltage=bidomain_problem.GetVoltage();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);
        
        /*
         * Test the top right node against the right one in the 1D case, 
         * comparing voltage, and then test all the nodes on the right hand 
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage;

        need_initialisation = true;

        // Test the RHF of the mesh
        for (int i = 0; i < bidomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (bidomain_problem.rGetMesh().GetNodeAt(i)->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to 
                //       to test the value of the nodes
                
                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[2*i];
                    need_initialisation = false;
                }
                else
                {
                    // the voltage at the end face varies a little because
                    // of drift due to the orientation of the tets in the mesh,
                    // hence the tolerance of 0.2
                    TS_ASSERT_DELTA(voltage_replicated[2*i], probe_voltage, 0.2);
                }
                
                // if a 1D simulation is run for 4ms on the 0_1mm_10elements mesh
                // the result at the end node is 20.0755
                TS_ASSERT_DELTA(voltage_replicated[2*i], 20.0755, 1.3);
            }
        }        

    }
            
    
    
    
    ////////////////////////////////////////////////////////////
    // Compare Mono and Bidomain Simulations
    ////////////////////////////////////////////////////////////
    void TestCompareBidomainProblemWithMonodomain3D()
    {
        ///////////////////////////////////////////////////////////////////
        // monodomain
        ///////////////////////////////////////////////////////////////////
        BidomainFaceStimulusCellFactory cell_factory;
        MonodomainProblem<3> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        monodomain_problem.SetEndTime(1);   // 1 ms
        monodomain_problem.SetOutputDirectory("Monodomain3d");
        monodomain_problem.SetOutputFilenamePrefix("monodomain3d");
        
        monodomain_problem.Initialise();
       
        // now solve       
        monodomain_problem.Solve();


        ///////////////////////////////////////////////////////////////////
        // bidomain
        ///////////////////////////////////////////////////////////////////
        BidomainProblem<3> bidomain_problem( &cell_factory );

        bidomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        bidomain_problem.SetEndTime(1);   // 1 ms
        bidomain_problem.SetOutputDirectory("Bidomain3d");
        bidomain_problem.SetOutputFilenamePrefix("bidomain3d");

        bidomain_problem.Initialise();
 
        // the bidomain equations reduce to the monodomain equations 
        // if sigma_e is infinite (equivalent to saying the extra_cellular
        // space is grounded. sigma_e is set to be very large here:
        c_matrix<double,3,3> sigma_i = bidomain_problem.GetBidomainPde()->GetIntracellularConductivityTensor();
        bidomain_problem.GetBidomainPde()->SetExtracellularConductivityTensor(10000*sigma_i);
                
        // now solve
        bidomain_problem.Solve();
        
                
        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        double* p_mono_voltage_array;
        double* p_bi_voltage_array;
        unsigned bi_lo, bi_hi, mono_lo, mono_hi;

        bidomain_problem.GetVoltageArray(&p_bi_voltage_array, bi_lo, bi_hi); 
        monodomain_problem.GetVoltageArray(&p_mono_voltage_array, mono_lo, mono_hi); 
        
        for (unsigned global_index=mono_lo; global_index<mono_hi; global_index++)
        {
            unsigned local_index = global_index - mono_lo;
            
            double monodomain_voltage      =   p_mono_voltage_array[local_index];
            double   bidomain_voltage      =   p_bi_voltage_array  [2*local_index];
            double extracellular_potential =   p_bi_voltage_array  [2*local_index+1];
            // std::cout << p_mono_voltage_array[local_index] << " " << p_bi_voltage_array[2*local_index] << "\n";

            // the mono and bidomains should agree closely 
            TS_ASSERT_DELTA(monodomain_voltage, bidomain_voltage, 0.5);
            
            // the extracellular potential should be uniform 
            TS_ASSERT_DELTA(extracellular_potential, 0, 0.1);
        } 
    }
    
};


#endif /*TESTBIDOMAIN3D_HPP_*/
