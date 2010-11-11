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
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PropagationPropertiesCalculator.hpp"

// stimulate an interval of cells
class BlockCellFactory1d : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BlockCellFactory1d()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        if ( fabs(x)<0.02+1e-6 )  
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
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
                    
                BlockCellFactory1d cell_factory;
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
                TetrahedralMesh<1,1> mesh;
                mesh.ConstructRegularSlabMesh(h[i], 1.0);
    
                std::stringstream output_dir;
                output_dir << "MonodomainSvi_" << h[i];
                HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
                HeartConfig::Instance()->SetOutputFilenamePrefix("results");
    
                HeartConfig::Instance()->SetUseStateVariableInterpolation();
    
                BlockCellFactory1d cell_factory;
                MonodomainProblem<1> monodomain_problem( &cell_factory );
                monodomain_problem.SetMesh(&mesh);
                monodomain_problem.Initialise();
            
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
};

#endif /*TESTMONODOMAINWITHSVI_HPP_*/
