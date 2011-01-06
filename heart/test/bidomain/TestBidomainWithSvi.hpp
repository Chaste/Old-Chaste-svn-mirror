/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTBIDOMAINWITHSVI_HPP_
#define TESTBIDOMAINWITHSVI_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "BidomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PropagationPropertiesCalculator.hpp"

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
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        double y;
        if(DIM==2)
        {
            y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];
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

class TestBidomainWithSvi : public CxxTest::TestSuite
{
public:
    void TestConductionVelocityConvergesFasterWithSvi1d() throw(Exception)
    {
        EXIT_IF_PARALLEL;
        
        double h[3] = {0.001,0.01,0.02};
        std::vector<double> conduction_vel_nci(3);
        std::vector<double> conduction_vel_svi(3);

        ReplicatableVector final_solution_nci;
        ReplicatableVector final_solution_svi;

        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        for(unsigned i=0; i<3; i++)
        {
            // NCI - nodal current interpolation - the default
            {
                TetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);
    
                std::stringstream output_dir;
                output_dir << "BidomainNci_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");
    
                // need to have this for i=1,2 cases!!
                HeartConfig::Instance()->SetUseStateVariableInterpolation(false);
                    
                BlockCellFactory<1> cell_factory;
                BidomainProblem<1> bidomain_problem( &cell_factory );
                bidomain_problem.SetMesh(&mesh);
                bidomain_problem.Initialise();
            
                bidomain_problem.Solve();
                    
                final_solution_nci.ReplicatePetscVector(bidomain_problem.GetSolution());
            }
       
            // SVI - state variable interpolation
            {
                TetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);
    
                std::stringstream output_dir;
                output_dir << "BidomainSvi_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");
    
                HeartConfig::Instance()->SetUseStateVariableInterpolation();
    
                BlockCellFactory<1> cell_factory;
                BidomainProblem<1> bidomain_problem( &cell_factory );
                bidomain_problem.SetMesh(&mesh);
                bidomain_problem.Initialise();

                if(i==0)
                {
                    bidomain_problem.UseMatrixBasedRhsAssembly(false);
                    TS_ASSERT_THROWS_CONTAINS(bidomain_problem.Solve(),"State variable interpolation only available");
                    bidomain_problem.UseMatrixBasedRhsAssembly(true);
                }
                            
                bidomain_problem.Solve();
                    
                final_solution_svi.ReplicatePetscVector(bidomain_problem.GetSolution());
            }

            double voltage_at_0_03_finest_mesh;
            if(i==0) // finest mesh
            {        
                for(unsigned j=0; j<final_solution_nci.GetSize(); j++)
                {
                    // visually checked they agree at this mesh resolution, and chosen tolerance from results
                    TS_ASSERT_DELTA(final_solution_nci[j], final_solution_svi[j], 0.35);
                    
                    if(j%2==0 /* just look at voltage */ && final_solution_nci[j]>-80) 
                    {
                        // shouldn't be exactly equal, as long as away from resting potential
                        TS_ASSERT_DIFFERS(final_solution_nci[j], final_solution_svi[j]);
                    }
                }
                
                voltage_at_0_03_finest_mesh = final_solution_nci[600];
                TS_ASSERT_DELTA(voltage_at_0_03_finest_mesh, -65.2218, 1e-3); //hardcoded value
            }
            else if(i==1)
            {
                double nci_voltage_at_0_03_middle_mesh = final_solution_nci[60];
                double svi_voltage_at_0_03_middle_mesh = final_solution_svi[60];
                // NCI conduction velocity > SVI conduction velocity
                // and both should be greater than CV on finesh mesh 
                TS_ASSERT_DELTA(nci_voltage_at_0_03_middle_mesh, -44.3111, 1e-3);
                TS_ASSERT_DELTA(svi_voltage_at_0_03_middle_mesh, -60.7765, 1e-3);
            }
            else
            {
                double nci_voltage_at_0_03_coarse_mesh = final_solution_nci[30];
                double svi_voltage_at_0_03_coarse_mesh = final_solution_svi[30];
                // NCI conduction velocity even greater than SVI conduction 
                // velocity
                TS_ASSERT_DELTA(nci_voltage_at_0_03_coarse_mesh,  -6.5622, 1e-3);
                TS_ASSERT_DELTA(svi_voltage_at_0_03_coarse_mesh, -51.8848, 1e-3);
            }
        }
    }
    
    void TestConductionVelocityInCrossFibreDirection2d() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        ReplicatableVector final_solution_nci;
        ReplicatableVector final_solution_svi;

        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.01);
        
        // much lower conductivity in cross-fibre direction - NCI will struggle
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.17));
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(7.0, 0.7));
        
        // NCI - nodal current interpolation - the default
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("BidomainNci2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation(false);
                
            BlockCellFactory<2> cell_factory;
            BidomainProblem<2> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();
        
            final_solution_nci.ReplicatePetscVector(bidomain_problem.GetSolution());
        }
   
        // SVI - state variable interpolation
        {
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRegularSlabMesh(0.02 /*h*/, 0.5, 0.3);

            HeartConfig::Instance()->SetOutputDirectory("BidomainSvi2d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            HeartConfig::Instance()->SetUseStateVariableInterpolation();

            BlockCellFactory<2> cell_factory;
            BidomainProblem<2> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();

            bidomain_problem.Solve();
                
            final_solution_svi.ReplicatePetscVector(bidomain_problem.GetSolution());
        }
        
        // See comments in equivalent part of test/monodomain/TestMonodomainWithSvi.hpp
        // For bidomain with h=0.02, SVI CV is slower in both fibre and cross-fibre 
        // directions

        // node 20 (for h=0.02) is on the x-axis (fibre direction)
        TS_ASSERT_DELTA(final_solution_nci[20*2], 22.1924, 1e-3);
        TS_ASSERT_DELTA(final_solution_svi[20*2], 14.7526, 1e-3);
        // node 234 (for h=0.02) is on the y-axis (cross-fibre direction)
        TS_ASSERT_DELTA(final_solution_nci[234*2], 19.5676, 1e-3);
        TS_ASSERT_DELTA(final_solution_svi[234*2], 13.6111, 1e-3);
    }
};

#endif /*TESTBIDOMAINWITHSVI_HPP_*/
