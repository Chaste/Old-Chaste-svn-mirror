/*

Copyright (C) University of Oxford, 2008

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
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "BidomainWithBathAssembler.hpp"

typedef BidomainWithBathAssembler<1,1> ASSEMBLER_1D;


template<unsigned DIM>
class BathCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    SimpleStimulus* mpStimulus;

public:
    BathCellFactory() : AbstractCardiacCellFactory<DIM>()
    {
        // set the new stimulus
        mpStimulus = new SimpleStimulus(-600, 0.5);
    }

    BathCellFactory(double stimulusMagnitude) : AbstractCardiacCellFactory<DIM>()
    {
        // set the new stimulus
        mpStimulus = new SimpleStimulus(stimulusMagnitude, 0.5);
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        // stimulate centre node normally.. 
        if (fabs(this->mpMesh->GetNode(node)->GetPoint()[0]-0.5)<1e-6)
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, mpStimulus, this->mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver, this->mpZeroStimulus, this->mpZeroStimulus);
        }
    }

    ~BathCellFactory(void)
    {
        delete mpStimulus;
    }
};







class TestBidomainWithBathAssembler : public CxxTest::TestSuite
{
public:
    void TestLabellingNodes() throw (Exception)
    {
        // all this is just to create a mesh, pde and bcc to pass to the new 
        // type of assembler
        HeartConfig::Instance()->SetSimulationDuration(0.01);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_10_elements_with_two_attributes");
        HeartConfig::Instance()->SetOutputDirectory("bidomainDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");
                
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );
        bidomain_problem.Initialise();
        
        AbstractMesh<1,1>* p_mesh = bidomain_problem.mpMesh;
        BidomainPde<1>* p_pde = bidomain_problem.mpBidomainPde;
        BoundaryConditionsContainer<1,1,2>  bcc;
        
        // Create the bidomain with bath assembler.
        BidomainWithBathAssembler<1,1> assembler(p_mesh, p_pde, &bcc);
        
        // the middle 4 elements are 'heart' elements (ie region=0),
        // so the middle 5 nodes should be heart nodes
        TS_ASSERT_EQUALS(p_mesh->GetNode(0)->GetRegion(), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(1)->GetRegion(), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(2)->GetRegion(), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(3)->GetRegion(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(4)->GetRegion(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(5)->GetRegion(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(6)->GetRegion(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(7)->GetRegion(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(8)->GetRegion(), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(9)->GetRegion(), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(10)->GetRegion(), 1u);

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
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );
        bidomain_problem.Initialise();
        
        AbstractMesh<1,1>* p_mesh = bidomain_problem.mpMesh;
        BidomainPde<1>* p_pde = bidomain_problem.mpBidomainPde;
        BoundaryConditionsContainer<1,1,2>  bcc;
        
        // Create the bidomain with bath assembler.
        // Fails because this mesh has no bath elements               
        TS_ASSERT_THROWS_ANYTHING(  ASSEMBLER_1D assembler(p_mesh, p_pde, &bcc) );
        
        // we need to call solve as other an EventHandling exception is thrown
        bidomain_problem.Solve();
    }
 

    void doesntpass_Test1dProblemOnlyBath() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(10.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBathOnlyBath");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath");
        
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> bidomain_cell_factory;
                
        BathCellFactory<1> cell_factory(-1e6);

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            mesh.GetElement(i)->SetRegion(1);
        }

        // create boundary conditions container
        BoundaryConditionsContainer<1,1,2> bcc;
        ConstBoundaryCondition<1>* p_bc_stim = new ConstBoundaryCondition<1>(1);
        ConstBoundaryCondition<1>* p_zero_stim = new ConstBoundaryCondition<1>(0.0);

        // loop over boundary elements
        for(TetrahedralMesh<1,1>::BoundaryElementIterator iter 
              = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            // if the element is on the left of the mesh, add a stimulus to the bcc
            if (((*iter)->GetNodeLocation(0))[0]==0.0)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_zero_stim); //note: I think you need to provide a boundary condition for unknown#1 if you are gonig to provide one for unknown#2? (todo)
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_stim, 1);
            }
        }

        BidomainProblem<1> bidomain_problem( &cell_factory, true );
        
        bidomain_problem.SetBoundaryConditionsContainer(&bcc);
        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Solve(); // diverges..?
    }
 
    
    void doesntwork_Test1DBathProblem() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(10.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_1d");
        
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> bidomain_cell_factory;
                
        BathCellFactory<1> cell_factory(-1e6);
        BidomainProblem<1> bidomain_problem( &cell_factory, true );

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        for(unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            if( (x<0.25) || (x>0.75) )
            {
                mesh.GetElement(i)->SetRegion(1);
            }
        }

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Solve();
        
        ///\todo test..
    }
};
    
#endif /*TESTBIDOMAINWITHBATHASSEMBLER_HPP_*/
