#ifndef TEST3DBIDOMAINWITHPERMFOREFFICIENCY_HPP_
#define TEST3DBIDOMAINWITHPERMFOREFFICIENCY_HPP_




#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "RegularStimulus.hpp"
#include "RandomNumberGenerator.hpp"
#include "BidomainFaceStimulusCellFactory.hpp"

class Test3dBidomainProblemWithPermForEfficiency :  public CxxTest::TestSuite
{
public:

    void TestBidomain3d() throw (Exception)
    {
        BidomainFaceStimulusCellFactory bidomain_cell_factory;
        
        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );
        
        bidomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");
        bidomain_problem.SetEndTime(150);   // ms
        bidomain_problem.PrintOutput(false);
        bidomain_problem.SetLinearSolverRelativeTolerance(1e-6);
        
        RandomNumberGenerator::Instance();
        bidomain_problem.rGetMesh().PermuteNodes();
        RandomNumberGenerator::Destroy();
        
        //PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "jacobi");
        PetscOptionsSetValue("-options_table", "");
        PetscOptionsSetValue("-log_summary", "");
        
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
        double probe_voltage=-9999.;
        
        need_initialisation = true;
        
        // Test the RHF of the mesh
        for (unsigned i = 0; i < bidomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (bidomain_problem.rGetMesh().GetNode(i)->GetPoint()[0] == 0.05)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.05cm is the other end of the mesh and where we want to
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
                    // hence the tolerance of 0.02
                    TS_ASSERT_DELTA(voltage_replicated[2*i], probe_voltage, 0.02);
                }
                
                // Check against hard coded value
                // For 50 ms test TS_ASSERT_DELTA(voltage_replicated[2*i],  7.3, 0.2);
                // For 150 ms test
                TS_ASSERT_DELTA(voltage_replicated[2*i],  -1.735, 0.1);
            }
        }
    }
};


#endif /*TEST3DBIDOMAINWITHPERMFOREFFICIENCY_HPP_*/
