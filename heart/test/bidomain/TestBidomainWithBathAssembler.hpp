/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef TESTBIDOMAINWITHBATHASSEMBLER_HPP_
#define TESTBIDOMAINWITHBATHASSEMBLER_HPP_


#include <cxxtest/TestSuite.h>
#include <vector>


#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "BidomainWithBathAssembler.hpp"
#include "TetrahedralMesh.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "HeartEventHandler.hpp"
#include "HeartRegionCodes.hpp"
#include "Timer.hpp"

typedef BidomainWithBathAssembler<1,1> ASSEMBLER_1D;


template<unsigned DIM>
class BathCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    c_vector<double,DIM> mStimulatedPoint;

public:
    BathCellFactory(double stimulusMagnitude, c_vector<double,DIM> stimulatedPoint)
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(stimulusMagnitude, 0.5)),
          mStimulatedPoint(stimulatedPoint)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        // paranoia - check this is really a tissue node
        assert(this->GetMesh()->GetNode(node)->GetRegion() == HeartRegionCode::TISSUE);

        // stimulate centre node normally..
        bool is_centre;

        if (DIM==1)
        {
            is_centre = (fabs(this->GetMesh()->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6);
        }
        else if (DIM==2)
        {
            is_centre = (    (fabs(this->GetMesh()->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6)
                          && (fabs(this->GetMesh()->GetNode(node)->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6) );
        }
        else
        {
            is_centre = (    (fabs(this->GetMesh()->GetNode(node)->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6)
                          && (fabs(this->GetMesh()->GetNode(node)->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6)
                          && (fabs(this->GetMesh()->GetNode(node)->GetPoint()[2]-mStimulatedPoint(2)) < 1e-6) );
        }

        if (is_centre)
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, this->mpZeroStimulus);
        }
    }

};




class TestBidomainWithBathAssembler : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Reset();
    }

    void TestLabellingNodes() throw (Exception)
    {
        // all this is just to create a mesh, pde and bcc to pass to the new
        // type of assembler
        HeartConfig::Instance()->SetSimulationDuration(0.01);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_10_elements_with_two_attributes");
        HeartConfig::Instance()->SetOutputDirectory("bidomainDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory, true );
        bidomain_problem.Initialise();

        AbstractTetrahedralMesh<1,1>* p_mesh = &(bidomain_problem.rGetMesh());
        BidomainPde<1>* p_pde = bidomain_problem.GetBidomainPde();
        BoundaryConditionsContainer<1,1,2> bcc;

        // Create the bidomain with bath assembler.
        BidomainWithBathAssembler<1,1> assembler(p_mesh, p_pde, &bcc);

        // the middle 4 elements are 'heart' elements (ie region=0),
        // so the middle 5 nodes should be heart nodes
        TS_ASSERT_EQUALS(p_mesh->GetNode(0)->GetRegion(), HeartRegionCode::BATH);
        TS_ASSERT_EQUALS(p_mesh->GetNode(1)->GetRegion(), HeartRegionCode::BATH);
        TS_ASSERT_EQUALS(p_mesh->GetNode(2)->GetRegion(), HeartRegionCode::BATH);
        TS_ASSERT_EQUALS(p_mesh->GetNode(3)->GetRegion(), HeartRegionCode::TISSUE);
        TS_ASSERT_EQUALS(p_mesh->GetNode(4)->GetRegion(), HeartRegionCode::TISSUE);
        TS_ASSERT_EQUALS(p_mesh->GetNode(5)->GetRegion(), HeartRegionCode::TISSUE);
        TS_ASSERT_EQUALS(p_mesh->GetNode(6)->GetRegion(), HeartRegionCode::TISSUE);
        TS_ASSERT_EQUALS(p_mesh->GetNode(7)->GetRegion(), HeartRegionCode::TISSUE);
        TS_ASSERT_EQUALS(p_mesh->GetNode(8)->GetRegion(), HeartRegionCode::BATH);
        TS_ASSERT_EQUALS(p_mesh->GetNode(9)->GetRegion(), HeartRegionCode::BATH);
        TS_ASSERT_EQUALS(p_mesh->GetNode(10)->GetRegion(), HeartRegionCode::BATH);

        // we need to call solve as other an EventHandling exception is thrown
        bidomain_problem.Solve();
    }


    void TestFailsIfNoBathElements() throw (Exception)
    {
        // all this is just to create a mesh, pde and bcc to pass to the new
        // type of assembler
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("bidomainDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory, true );
        // Fails because no bath
        TS_ASSERT_THROWS_THIS(bidomain_problem.Initialise(), "No bath element found");

        // Prevent an EventHandling exception in later tests
        HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
    }


    void TestBathIntracellularStimulation() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(10.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_1d");

        c_vector<double,1> centre;
        centre(0) = 0.5;
        BathCellFactory<1> cell_factory(-1e6, centre); // stimulates x=0.5 node

        BidomainProblem<1> bidomain_problem( &cell_factory, true );

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        // set the x<0.25 and x>0.75 regions as the bath region
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            if( (x<0.25) || (x>0.75) )
            {
                mesh.GetElement(i)->SetRegion(HeartRegionCode::BATH);
            }
        }

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        // test V = 0 for all bath nodes
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if(mesh.GetNode(i)->GetRegion()==HeartRegionCode::BATH) // bath
            {
                TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
            }
        }

        // test symmetry of V and phi_e
        for(unsigned i=0; i<=(mesh.GetNumNodes()-1)/2; i++)
        {
            unsigned opposite = mesh.GetNumNodes()-i-1;
            TS_ASSERT_DELTA(sol_repl[2*i], sol_repl[2*opposite], 1e-4);      // V
            TS_ASSERT_DELTA(sol_repl[2*i+1], sol_repl[2*opposite+1], 1e-4);  // phi_e
        }

        // a couple of hardcoded values
        TS_ASSERT_DELTA(sol_repl[2*50], 3.7684, 1e-3);
        TS_ASSERT_DELTA(sol_repl[2*70], 5.1777, 1e-3);
    }


    // In this test we have no cardiac tissue, so that the equations are just sigma * phi_e''=0
    // throughout the domain (with a Neumann boundary condition on x=1 and a dirichlet boundary
    // condition (ie grounding) on x=0), so the exact solution can be calculated and compared
    // against.
    void Test1dProblemOnlyBathGroundedOneSide() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(0.5);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBathOnlyBath");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> bidomain_cell_factory;

        c_vector<double,1> centre;
        centre(0) = 0.5;
        BathCellFactory<1> cell_factory(-1e6, centre);

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            mesh.GetElement(i)->SetRegion(HeartRegionCode::BATH);
        }

        // create boundary conditions container
        double boundary_val = 1.0;
        boost::shared_ptr<BoundaryConditionsContainer<1,1,2> > p_bcc(new BoundaryConditionsContainer<1,1,2>);
        ConstBoundaryCondition<1>* p_bc_stim = new ConstBoundaryCondition<1>(boundary_val);
        ConstBoundaryCondition<1>* p_zero_stim = new ConstBoundaryCondition<1>(0.0);

        // loop over boundary elements and set (sigma\gradphi).n = 1.0 on RHS edge
        for(TetrahedralMesh<1,1>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if (((*iter)->GetNodeLocation(0))[0]==1.0)
            {
                /// \todo: I think you need to provide a boundary condition for unknown#1 if you are gonig to provide one for unknown#2?
                p_bcc->AddNeumannBoundaryCondition(*iter, p_zero_stim, 0);
                p_bcc->AddNeumannBoundaryCondition(*iter, p_bc_stim,   1);
            }
        }

        BidomainProblem<1> bidomain_problem( &cell_factory, true );

        bidomain_problem.SetBoundaryConditionsContainer(p_bcc);
        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        // fix phi=0 on LHS edge
        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed_nodes);

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        // test phi = x*boundary_val/sigma (solution of phi''=0, phi(0)=0, sigma*phi'(1)=boundary_val
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double bath_cond = HeartConfig::Instance()->GetBathConductivity();
            double x = mesh.GetNode(i)->rGetLocation()[0];
            TS_ASSERT_DELTA(sol_repl[2*i],   0.0,   1e-12);               // V
            TS_ASSERT_DELTA(sol_repl[2*i+1], x*boundary_val/bath_cond, 1e-4);   // phi_e
        }
    }

    void Test2dBathIntracellularStimulation() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d");

        c_vector<double,2> centre;
        centre(0) = 0.05;
        centre(1) = 0.05;
        BathCellFactory<2> cell_factory(-5e6, centre); // stimulates x=0.05 node

        BidomainProblem<2> bidomain_problem( &cell_factory, true );

        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        ParallelTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Set everything outside a central circle (radius 0.4) to be bath
        //for(unsigned i=0; i<mesh.GetNumElements(); i++)
        for (ParallelTetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            double x = it->CalculateCentroid()[0];
            double y = it->CalculateCentroid()[1];
            if ( sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05)) > 0.04 )
            {
                it->SetRegion(HeartRegionCode::BATH);
            }
        }

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        // test V = 0 for all bath nodes
        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter=mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd(); ++iter)
        {
            if ((*iter).GetRegion()==HeartRegionCode::BATH) // bath
            {
                unsigned index=(*iter).GetIndex();

                TS_ASSERT_DELTA(sol_repl[2*index], 0.0, 1e-12);
            }
        }

        const std::vector<unsigned>& permutation = mesh.rGetNodePermutation();

        unsigned node_50;
        unsigned node_70;

        // In sequential the permutation vector is empty
        if (permutation.size() > 0)
        {
            node_50 = permutation[50];
            node_70 = permutation[70];
        }
        else
        {
            node_50 = 50;
            node_70 = 70;
        }

        // a couple of hardcoded value
        TS_ASSERT_DELTA(sol_repl[2*node_50], 28.3912, 1e-3);
        TS_ASSERT_DELTA(sol_repl[2*node_70], 28.3912, 1e-3);
    }

    void Test2dBathInputFluxEqualsOutputFlux() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(3.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath2dFluxCompare");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_2d_fluxes");

        HeartConfig::Instance()->SetOdeTimeStep(0.001);  //ms

        // need to create a cell factory but don't want any intra stim, so magnitude
        // of stim is zero.
        c_vector<double,2> centre;
        centre(0) = 0.05; // cm
        centre(1) = 0.05; // cm
        BathCellFactory<2> cell_factory( 0.0, centre);

        BidomainProblem<2> bidomain_problem( &cell_factory, true );

        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Set everything outside a central circle (radius 0.04cm) to be bath
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            double y = mesh.GetElement(i)->CalculateCentroid()[1];
            if( sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05)) > 0.02 )
            {
                mesh.GetElement(i)->SetRegion(HeartRegionCode::BATH);
            }
        }

        //boundary flux for Phi_e. -10e3 is under thershold, -14e3 crashes the cell model
        double boundary_flux = -11.0e3;
        double duration = 1.9; // of the stimulus, in ms

        Electrodes<2> electrodes(mesh,false,0,0.0,0.1,boundary_flux, duration);

        // Cover an exception
        {
            // Avoid event handler exception
            HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
            BidomainProblem<2> no_bath_problem( &cell_factory, false );
            TS_ASSERT_THROWS_THIS(no_bath_problem.SetElectrodes(electrodes),
                    "Cannot set electrodes when problem has been defined to not have a bath");
        }

        bidomain_problem.SetElectrodes(electrodes);

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        bool ap_triggered = false;
        /*
         * We are checking the last time step. This test will only make sure that an upstroke is triggered.
         * We ran longer simulation for 350 ms and a nice AP was observed.
         */

        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            // test V = 0 for all bath nodes and that an AP is triggered in the tissue
            if (mesh.GetNode(i)->GetRegion() == HeartRegionCode::BATH) // bath
            {
                TS_ASSERT_DELTA(sol_repl[2*i], 0.0, 1e-12);
            }
            else if (sol_repl[2*i] > 0.0)//at the last time step
            {
                ap_triggered = true;
            }
        }

        TS_ASSERT_EQUALS(electrodes.mAreActive, false); // should be switched off by now..
        TS_ASSERT(ap_triggered);
    }

    void TestMatrixBasedAssembledBath(void)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms

        // need to create a cell factory but don't want any intra stim, so magnitude
        // of stim is zero.
        c_vector<double,2> centre;
        centre(0) = 0.05;
        centre(1) = 0.05;
        BathCellFactory<2> cell_factory( 0.0, centre);

        //boundary flux for Phi_e
        double boundary_flux = -4e2;
        double duration = 0.2; //ms

        TrianglesMeshReader<2,2> reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        // Set everything outside a central circle (radius 0.4) to be bath
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            double y = mesh.GetElement(i)->CalculateCentroid()[1];
            if( sqrt((x-0.05)*(x-0.05) + (y-0.05)*(y-0.05)) > 0.04 )
            {
                mesh.GetElement(i)->SetRegion(HeartRegionCode::BATH);
            }
        }

        ///////////////////////////////////////////////////////////////////
        // matrix based
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("BidomainBathMatrixBased");
        HeartConfig::Instance()->SetOutputFilenamePrefix("matrix_based");

        BidomainProblem<2> matrix_based_bido( &cell_factory, true );

        {
            Timer::Reset();

            Electrodes<2> electrodes(mesh,false,0,0.0,0.1,boundary_flux, duration);

            matrix_based_bido.SetElectrodes(electrodes);
            matrix_based_bido.SetMesh(&mesh);
            matrix_based_bido.Initialise();
            matrix_based_bido.Solve();

            Timer::Print("2D Matrix based");
        }

        ///////////////////////////////////////////////////////////////////
        // non matrix based
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("BidomainBathNonMatrixBased");
        HeartConfig::Instance()->SetOutputFilenamePrefix("non_matrix_based");

        BidomainProblem<2> non_matrix_based_bido( &cell_factory, true);

        {
            Timer::Reset();

            Electrodes<2> electrodes(mesh,false,0,0.0,0.1,boundary_flux, duration);

            non_matrix_based_bido.SetElectrodes(electrodes);
            non_matrix_based_bido.SetMesh(&mesh);
            non_matrix_based_bido.UseMatrixBasedRhsAssembly(false);
            non_matrix_based_bido.Initialise();
            non_matrix_based_bido.Solve();

            Timer::Print("2D non matrix based");
        }

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector matrix_based_solution = matrix_based_bido.GetSolutionDistributedVector();
        DistributedVector non_matrix_based_solution = non_matrix_based_bido.GetSolutionDistributedVector();

        DistributedVector::Stripe matrix_based_voltage(matrix_based_solution, 0);
        DistributedVector::Stripe non_matrix_based_voltage(non_matrix_based_solution, 0);

        DistributedVector::Stripe matrix_based_ex_pot(matrix_based_solution, 1);
        DistributedVector::Stripe non_matrix_based_ex_pot(non_matrix_based_solution, 1);

        for (DistributedVector::Iterator index = matrix_based_solution.Begin();
             index != matrix_based_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(matrix_based_voltage[index], non_matrix_based_voltage[index], 1e-7);
            //TS_ASSERT_DELTA(matrix_based_ex_pot[index], non_matrix_based_ex_pot[index], 1e-7);
            //std::cout << matrix_based_voltage[index] << std::endl;
        }
    }

};

#endif /*TESTBIDOMAINWITHBATHASSEMBLER_HPP_*/
