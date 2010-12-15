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

#ifndef TESTMONODOMAINWITHSVI_HPP_
#define TESTMONODOMAINWITHSVI_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "MonodomainProblem.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "MatrixBasedMonodomainSolver.hpp"
#include "TenTusscher2006Epi.hpp"
#include "Mahajan2008.hpp"

// stimulate a block of cells (an interval in 1d, a block in a corner in 2d)
template<unsigned DIM>
class BlockCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BlockCellFactory()
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
        assert(DIM<3);
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNodeOrHaloNode(nodeIndex)->rGetLocation()[0];
        double y;
        if(DIM==2)
        {
            y = this->GetMesh()->GetNodeOrHaloNode(nodeIndex)->rGetLocation()[1];
        }

        if (    (DIM==1 && fabs(x)<0.02+1e-6)
             || (DIM==2 && fabs(x)<0.02+1e-6 && fabs(y)<0.02+1e-6) )  
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};


// non-identical cell models
class HeterogeneousCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    HeterogeneousCellFactory()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5)) 
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNodeOrHaloNode(nodeIndex)->rGetLocation()[0];
        if ( x<0.15 )  
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else if (x < 0.65 )
        {
            return new CellMahajan2008FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
        else
        {
            return new CellTenTusscher2006EpiFromCellML(this->mpSolver, this->mpZeroStimulus);
        }            
    }
};


class TestMonodomainWithSvi : public CxxTest::TestSuite
{
public:
    void TestConductionVelocityConvergesFasterWithSvi1d() throw(Exception)
    {
        EXIT_IF_PARALLEL;
        
        double h[3] = {0.001,0.01,0.02};
        std::vector<double> conduction_vel_nci(3);
        std::vector<double> conduction_vel_svi(3);

        ReplicatableVector final_voltage_nci;
        ReplicatableVector final_voltage_svi;

        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        for(unsigned i=0; i<3; i++)
        {
            // NCI - nodal current interpolation - the default
            {
                TetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);
    
                std::stringstream output_dir;
                output_dir << "MonodomainNci_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");
    
                // need to have this for i=1,2 cases!!
                HeartConfig::Instance()->SetUseStateVariableInterpolation(false);
                    
                BlockCellFactory<1> cell_factory;
                MonodomainProblem<1> monodomain_problem( &cell_factory );
                monodomain_problem.SetMesh(&mesh);
                monodomain_problem.Initialise();
            
                monodomain_problem.Solve();
                    
                final_voltage_nci.ReplicatePetscVector(monodomain_problem.GetSolution());

//// see #1633
//// end time needs to be increased for these (say, to 7ms)
//                Hdf5DataReader simulation_data(OutputFileHandler::GetChasteTestOutputDirectory() + output_dir.str(),
//                                               "results", false);
//                PropagationPropertiesCalculator ppc(&simulation_data);    
//                unsigned node_at_0_04 = (unsigned)round(0.04/h[i]);
//                unsigned node_at_0_40 = (unsigned)round(0.40/h[i]);
//                assert(fabs(mesh.GetNode(node_at_0_04)->rGetLocation()[0]-0.04)<1e-6);
//                assert(fabs(mesh.GetNode(node_at_0_40)->rGetLocation()[0]-0.40)<1e-6);
//                conduction_vel_nci[i] = ppc.CalculateConductionVelocity(node_at_0_04,node_at_0_40,0.36);        
//                std::cout << "conduction_vel_nci = " << conduction_vel_nci[i] << "\n";
            }
       
            // SVI - state variable interpolation
            {
                DistributedTetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);
    
                std::stringstream output_dir;
                output_dir << "MonodomainSvi_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");
    
                HeartConfig::Instance()->SetUseStateVariableInterpolation();
    
                BlockCellFactory<1> cell_factory;
                MonodomainProblem<1> monodomain_problem( &cell_factory );
                monodomain_problem.SetMesh(&mesh);
                monodomain_problem.Initialise();

                if(i==0)
                {
                    monodomain_problem.UseMatrixBasedRhsAssembly(false);
                    TS_ASSERT_THROWS_CONTAINS(monodomain_problem.Solve(),"State variable interpolation only available");
                    monodomain_problem.UseMatrixBasedRhsAssembly(true);
                }
                            
                monodomain_problem.Solve();
                    
                final_voltage_svi.ReplicatePetscVector(monodomain_problem.GetSolution());

//                Hdf5DataReader simulation_data(OutputFileHandler::GetChasteTestOutputDirectory() + output_dir.str(),
//                                               "results", false);
//                PropagationPropertiesCalculator ppc(&simulation_data);
//                unsigned node_at_0_04 = (unsigned)round(0.04/h[i]);
//                unsigned node_at_0_40 = (unsigned)round(0.40/h[i]);
//                assert(fabs(mesh.GetNode(node_at_0_04)->rGetLocation()[0]-0.04)<1e-6);
//                assert(fabs(mesh.GetNode(node_at_0_40)->rGetLocation()[0]-0.40)<1e-6);    
//                conduction_vel_svi[i] = ppc.CalculateConductionVelocity(node_at_0_04,node_at_0_40,0.36); 
//                std::cout << "conduction_vel_svi = " << conduction_vel_svi[i] << "\n";
            }

            double voltage_at_0_03_finest_mesh;
            if(i==0) // finest mesh
            {        
                for(unsigned j=0; j<final_voltage_nci.GetSize(); j++)
                {
                    // visually checked they agree at this mesh resolution, and chosen tolerance from results
                    TS_ASSERT_DELTA(final_voltage_nci[j], final_voltage_svi[j], 0.3);
                    
                    if(final_voltage_nci[j]>-80)
                    {
                        // shouldn't be exactly equal, as long as away from resting potential
                        TS_ASSERT_DIFFERS(final_voltage_nci[j], final_voltage_svi[j]);
                    }
                }
                
                voltage_at_0_03_finest_mesh = final_voltage_nci[300];
                TS_ASSERT_DELTA(voltage_at_0_03_finest_mesh, 11.1182, 1e-3); //hardcoded value
            }
            else if(i==1)
            {
                double nci_voltage_at_0_03_middle_mesh = final_voltage_nci[30];
                double svi_voltage_at_0_03_middle_mesh = final_voltage_svi[30];
                // NCI conduction velocity > SVI conduction velocity
                // and both should be greater than CV on finesh mesh 
                TS_ASSERT_DELTA(nci_voltage_at_0_03_middle_mesh, 19.8924, 1e-3);
                TS_ASSERT_DELTA(svi_voltage_at_0_03_middle_mesh, 14.9579, 1e-3);
            }
            else
            {
                double nci_voltage_at_0_03_coarse_mesh = final_voltage_nci[15];
                double svi_voltage_at_0_03_coarse_mesh = final_voltage_svi[15];
                // NCI conduction velocity even greater than SVI conduction 
                // velocity
                TS_ASSERT_DELTA(nci_voltage_at_0_03_coarse_mesh, 24.4938, 1e-3);
                TS_ASSERT_DELTA(svi_voltage_at_0_03_coarse_mesh, 17.3131, 1e-3);
            }
        }
    }
    
    void TestConductionVelocityInCrossFibreDirection2d() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        ReplicatableVector final_voltage_nci;
        ReplicatableVector final_voltage_svi;

        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.01);
        
        // much lower conductivity in cross-fibre direction - NCI will struggle
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.17));
        
        // NCI - nodal current interpolation - the default
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("MonodomainNci2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation(false);
                
            BlockCellFactory<2> cell_factory;
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();
            monodomain_problem.Solve();
        
            final_voltage_nci.ReplicatePetscVector(monodomain_problem.GetSolution());
        }
   
        // SVI - state variable interpolation
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("MonodomainSvi2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation();

            BlockCellFactory<2> cell_factory;
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();

            monodomain_problem.Solve();
                
            final_voltage_svi.ReplicatePetscVector(monodomain_problem.GetSolution());
        }
        
        // Visualised results with h=0.02 and h=0.01 - results looks sensible according to 
        // paper:
        // 1. SVI h=0.01 and h=0.02 match more closely than NCI results - ie SVI converges faster
        // 2. CV in fibre direction faster for NCI (both values of h)
        // 3. CV in cross fibre direction: (i) h=0.01, faster for NCI; h=0.02 slower for NCI. 
        // (Matches results in paper)
        
        // node 20 (for h=0.02) is on the x-axis (fibre direction), SVI CV is slower
        TS_ASSERT_DELTA(final_voltage_nci[20], -9.2269, 1e-3);
        TS_ASSERT_DELTA(final_voltage_svi[20], -60.8510, 1e-3);
        // node 130 (for h=0.02) is on the y-axis (cross-fibre direction), NCI CV is slower
        TS_ASSERT_DELTA(final_voltage_nci[130], 14.7918, 1e-3);
        TS_ASSERT_DELTA(final_voltage_svi[130], 30.5281, 1e-3);
    }
    
    void TestCoverage3d() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        HeartConfig::Instance()->SetSimulationDuration(0.1); //ms
        HeartConfig::Instance()->SetUseStateVariableInterpolation(true);

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructRegularSlabMesh(0.02, 0.02, 0.02, 0.02);

        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML,3> cell_factory;
        MonodomainProblem<3> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();
        monodomain_problem.Solve();
    }

    void TestWithHeterogeneousCellModels() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
        HeartConfig::Instance()->SetUseStateVariableInterpolation(true);

        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.01, 1.0);

        HeterogeneousCellFactory cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();
        
        //// really very difficult to get hold of the correction assembler from here.
        //// The following, which requires 2 classes to be friends of this test compiles
        //// but fails asserts in the first line as bcc is not set up. 
        //AbstractDynamicLinearPdeSolver<1,1,1>* p_solver = monodomain_problem.CreateSolver();
        //MatrixBasedMonodomainSolver<1,1>* p_mono_solver = dynamic_cast<MatrixBasedMonodomainSolver<1,1>*>(p_solver);
        //MonodomainCorrectionTermAssembler<1,1>* p_assembler = p_mono_solver->mpMonodomainCorrectionTermAssembler;
        //TS_ASSERT_EQUALS(p_assembler->mElementsHasIdenticalCellModels.size(), 10u);        
        
        // therefore, we just test that calling Solve() runs (without the checking that
        // cell models are identical, this fails).

//#1462
// currently this fails (as expected) with an isnan assert, as the checking of identical cell models is
// not yet implemented (see constructor of AbstractCorrectionTermAssembler).
// Implement and uncomment
//        monodomain_problem.Solve();
    }
};

#endif /*TESTMONODOMAINWITHSVI_HPP_*/
