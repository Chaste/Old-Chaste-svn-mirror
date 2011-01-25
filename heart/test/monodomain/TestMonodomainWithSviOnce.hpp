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
#include "MatrixBasedMonodomainSolver.hpp"


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


//This test passes on 1,2 processes
//Fail on 3: Note that x=0.3 is the probe position, which is near the proc 0/1 boundary

///\todo #1462 Test a bit more thoroughly -- then delete this suite and move onto the ICI/SVI test suites.
class TestMonodomainWithSvi : public CxxTest::TestSuite
{
public:
    void TestWithSvi1d() throw(Exception)
    {
        double h = 0.02;
        unsigned probe_node_index =  15;

        ReplicatableVector final_voltage_svi;

        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-8);
        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(h, 1.0);

        //Double check (for later) that the indexing is as expected
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal( probe_node_index ))
        {
            TS_ASSERT_DELTA(mesh.GetNode( probe_node_index )->rGetLocation()[0], 0.3, 1e-8);
        }
        std::stringstream output_dir;
        output_dir << "MonodomainSvi_" << h<<"_"<<PetscTools::GetNumProcs();
        HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        HeartConfig::Instance()->SetUseStateVariableInterpolation();

        BlockCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();
                    
        monodomain_problem.Solve();
            
        final_voltage_svi.ReplicatePetscVector(monodomain_problem.GetSolution());


        double svi_voltage_at_0_3 = final_voltage_svi[ probe_node_index ];
        double hard_coded_at_0_3 = 17.313147483511354352; //hardcoded value from this mesh with svi sequential
        std::cout<<std::setprecision(20)<<svi_voltage_at_0_3<<"\t"<<svi_voltage_at_0_3 - hard_coded_at_0_3<<"\n";
        TS_ASSERT_DELTA(svi_voltage_at_0_3 - hard_coded_at_0_3, 0.0, 1e-13);
    }
    void TestWithIci1d() throw(Exception)
    {
        //This gives and indication of the drift that more processes adds
        double h = 0.02;
        unsigned probe_node_index =  15;

        ReplicatableVector final_voltage_svi;

        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-8);
        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);

        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(h, 1.0);

        //Double check (for later) that the indexing is as expected
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal( probe_node_index ))
        {
            TS_ASSERT_DELTA(mesh.GetNode( probe_node_index )->rGetLocation()[0], 0.3, 1e-8);
        }
        std::stringstream output_dir;
        output_dir << "MonodomainSvi_" << h<<"_"<<PetscTools::GetNumProcs();
        HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        HeartConfig::Instance()->SetUseStateVariableInterpolation(false);

        BlockCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();
                    
        monodomain_problem.Solve();
            
        final_voltage_svi.ReplicatePetscVector(monodomain_problem.GetSolution());


        double svi_voltage_at_0_3 = final_voltage_svi[ probe_node_index ];
        double hard_coded_at_0_3 = 24.493851957613916426; //hardcoded value from this mesh with svi sequential
        std::cout<<std::setprecision(20)<<svi_voltage_at_0_3<<"\t"<<svi_voltage_at_0_3 - hard_coded_at_0_3<<"\n";
        TS_ASSERT_DELTA(svi_voltage_at_0_3 - hard_coded_at_0_3, 0.0, 1e-13);
    }
    
};

#endif /*TESTMONODOMAINWITHSVI_HPP_*/
