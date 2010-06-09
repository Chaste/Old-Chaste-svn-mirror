#ifndef TESTBENCHMARK_HPP_
#define TESTBENCHMARK_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"



class BenchmarkCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BenchmarkCellFactory()
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-50000.0, 2))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        double y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];
        double z = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[2];

        if ( (x<0.15+1e-6) && (y<0.15+1e-6) && (z<0.15+1e-6) ) 
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }
    }
};
    



class TestBenchmark : public CxxTest::TestSuite
{

public:
    void TestMonodomain( void ) throw(Exception)
    {
        TetrahedralMesh<3,3> mesh;
        double h=0.01; //cm
        unsigned num_elem_x = (unsigned)(2/h);
        unsigned num_elem_y = (unsigned)(0.7/h);
        unsigned num_elem_z = (unsigned)(0.3/h);
        
        std::cout << num_elem_x << " " << num_elem_y << " " << num_elem_z << "\n";  
        
        mesh.ConstructCuboid(num_elem_x, num_elem_y, num_elem_z);
        mesh.Scale(h, h, h); 
        c_vector<double,6> extremes = mesh.GetExtremes();
        std::cout << extremes(0) << " "
                  << extremes(1) << " "
                  << extremes(2) << " "
                  << extremes(3) << " "
                  << extremes(4) << " "
                  << extremes(5) << "\n";
                  
        HeartConfig::Instance()->SetOutputDirectory("TestBenchmark");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");          
 
        HeartConfig::Instance()->SetSimulationDuration(50); //ms
        
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400); // 1400 1/cm
        HeartConfig::Instance()->SetCapacitance(1); // 1uF/cm^2
 
        double long_conductance = 0.17 * 0.62/(0.17+0.62) * 10; // harmonic mean of 0.17,0.62 S/m convert mS/cm
        double trans_conductance = 0.019 * 0.24/(0.019+0.24) * 10; // harmonic mean of 0.019,0.24 S/m convert mS/cm
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(long_conductance, trans_conductance, trans_conductance));
 
        BenchmarkCellFactory cell_factory;
        
        MonodomainProblem<3> problem( &cell_factory );

        problem.SetMesh(&mesh);
        
        problem.Initialise();
        problem.SetWriteInfo();
        
        problem.Solve();
    }
};


#endif /*TESTBENCHMARK_HPP_*/
