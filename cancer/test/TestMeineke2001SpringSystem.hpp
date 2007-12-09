#ifndef TESTMEINEKE2001SPRINGSYSTEM_HPP_
#define TESTMEINEKE2001SPRINGSYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <cmath>
#include <vector>
#include "Meineke2001SpringSystem.hpp"
#include "TrianglesMeshReader.cpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "WntGradient.hpp"
#include "VoronoiTessellation.cpp"
#include "CellsGenerator.hpp"

class TestMeineke2001SpringSystem : public CxxTest::TestSuite
{
public:
    void TestForceCalculations() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        CancerParameters *p_params = CancerParameters::Instance();
 
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;
        
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellMutationState mutation_state = HEALTHY;
            if (i==60)
            {
                mutation_state = APC_TWO_HIT;
            }
            
            TissueCell cell(STEM, mutation_state, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(-10);
            cells.push_back(cell);
        }

        Tissue<2> tissue(*p_mesh, cells);
        tissue.SetGhostNodes(ghost_node_indices);

        Meineke2001SpringSystem<2> meineke_spring_system(tissue);
        
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test Calculate Velocities on each node
         ************************************************************************
         ************************************************************************ 
         */
                
        std::vector<c_vector<double, 2> >& velocities_on_each_node = meineke_spring_system.rCalculateVelocitiesOfEachNode();
 
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            std::set<unsigned>::iterator iter = ghost_node_indices.find(i);
            bool is_a_ghost_node = (iter!=ghost_node_indices.end());

            if (!is_a_ghost_node)
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], 0.0, 1e-4);
                TS_ASSERT_DELTA(velocities_on_each_node[i][1], 0.0, 1e-4);
            }
        }
        
        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point = p_mesh->GetNode(59)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];
  
        p_mesh->SetNode(59, new_point, false);
        velocities_on_each_node = meineke_spring_system.rCalculateVelocitiesOfEachNode();
        
        TS_ASSERT_DELTA(velocities_on_each_node[60][0], 0.5*p_params->GetSpringStiffness()/p_params->GetDampingConstantMutant(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[59][0], (-3+4.0/sqrt(7))*p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[58][0], 0.5*p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[58][1], 0.0, 1e-4);
        
        
        /*
         ************************************************************************
         ************************************************************************ 
         *  Test Calculate force on a spring
         ************************************************************************
         ************************************************************************ 
         */
        
        c_vector<double,2> force_on_spring ; // between nodes 59 and 60
        
        // Find one of the elements that nodes 59 and 60 live on
        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01 ;
        
        unsigned elem_index = p_mesh->GetContainingElementIndex(new_point2,false);
        Element<2,2>* p_element = p_mesh->GetElement(elem_index);
        
        force_on_spring = meineke_spring_system.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),p_element->GetNodeGlobalIndex(0));
        TS_ASSERT_DELTA(force_on_spring[0], 0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);        
        
        /*
         ******************************************
         ******************************************
         *  test force with cutoff point
         ******************************************
         ******************************************
         */
        double dist = norm_2( p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(), p_element->GetNode(1)->rGetLocation()) );   

        meineke_spring_system.UseCutoffPoint(dist-0.1);
        
        force_on_spring = meineke_spring_system.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),p_element->GetNodeGlobalIndex(0));
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestEdgeLengthBasedSpring() throw (Exception)
    {        
        CancerParameters::Instance()->Reset();

        // Set up the simulation time
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
                
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;
       
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false,crypt_width/cells_across);
        ConformingTetrahedralMesh<2,2>* p_mesh =generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();      
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Tissue<2> tissue(*p_mesh, cells);               
        tissue.SetGhostNodes(ghost_node_indices);

        Meineke2001SpringSystem<2> meineke_spring_system(tissue);
               
        // check that the force between nodes is correctly calculated when the spring constant is constant (!)
        meineke_spring_system.SetEdgeBasedSpringConstant(false);
                      
        for(Tissue<2>::SpringIterator spring_iterator=tissue.SpringsBegin();
            spring_iterator!=tissue.SpringsEnd();
            ++spring_iterator)
        {        
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = meineke_spring_system.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],6.25,1e-3);
        }
        
        // check that the force between nodes is correctly calculated when the spring constant 
        // is proportional to the length of the edge between adjacent cells  
        meineke_spring_system.SetEdgeBasedSpringConstant(true); 
        tissue.CreateVoronoiTessellation();  // normally done in a simulation loop
        
        
        for(Tissue<2>::SpringIterator spring_iterator=tissue.SpringsBegin();
            spring_iterator!=tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = meineke_spring_system.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],4.34027778,1e-3);
        }
        
        // choose two interior neighbour nodes
        c_vector<double, 2> force = meineke_spring_system.CalculateForceBetweenNodes(20u,21u);
        TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],4.34027778,1e-3);
        
        // now move node 21 a bit and check that the force calculation changes correctly
        c_vector<double,2> shift;
        shift[0] = 0.01;
        shift[1] = 0.0;                 
        ChastePoint<2> new_point(p_mesh->GetNode(21u)->rGetLocation() + shift);
        p_mesh->SetNode(21u, new_point, false);
        
        // check that the new force between nodes is correctly calculated
        tissue.CreateVoronoiTessellation();  
        c_vector<double, 2> new_force = meineke_spring_system.CalculateForceBetweenNodes(20u,21u);
        
        // force calculation: shift is along x-axis so we should have
        // new_edge_length = (5/6 + shift[0])*tan(0.5*arctan(5*sqrt(3)/(5 + 12*shift[0]))),
        // force^2 = mu^2 * (new_edge_length*sqrt(3))^2 * (1 - 5/6 - shift[0])^2
        TS_ASSERT_DELTA(new_force[0]*new_force[0] + new_force[1]*new_force[1], 3.83479824,1e-3);
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();    
    }
    

    void TestEdgeBasedSpringsOnPeriodicMesh() throw (Exception)
    {     
        // Test on a periodic mesh
        CancerParameters::Instance()->Reset();

        // Set up the simulation time
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

                
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;
       
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Tissue<2> tissue(*p_mesh, cells);               
        tissue.SetGhostNodes(ghost_node_indices);

        Meineke2001SpringSystem<2> meineke_spring_system(tissue);
       
         // check that the force between nodes is correctly calculated when the spring constant is constant (!)
        meineke_spring_system.SetEdgeBasedSpringConstant(false);
                      
        for(Tissue<2>::SpringIterator spring_iterator=tissue.SpringsBegin();
            spring_iterator!=tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = meineke_spring_system.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
                        
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],6.25,1e-3);
        }
        
        // check that the force between nodes is correctly calculated when the spring constant 
        // is proportional to the length of the edge between adjacenet cells  
        meineke_spring_system.SetEdgeBasedSpringConstant(true); 
        tissue.CreateVoronoiTessellation();  
        for(Tissue<2>::SpringIterator spring_iterator=tissue.SpringsBegin();
            spring_iterator!=tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = meineke_spring_system.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1],4.34027778,1e-3);
        }              
            
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();                
    }
    

    void TestAreaBasedVisocity() throw (Exception)
    {        
        CancerParameters::Instance()->Reset();

        // Set up the simulation time
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
             
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 2;
       
        // Test a non-periodic mesh        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        
        // scale the mesh in one direction
        p_mesh->Scale(1.0,1.2);
        
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();      
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Tissue<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        crypt.CreateVoronoiTessellation();  // normally done in a simulation loop

        Meineke2001SpringSystem<2> meineke_spring_system(crypt);
                             
        std::vector<c_vector<double,2> >& velocities = meineke_spring_system.rCalculateVelocitiesOfEachNode();
        std::vector<double> norm_vel;
        
        for(unsigned i=0; i<velocities.size(); i++)
        {
            //check if this is a real cell
            if(ghost_node_indices.find(i)==ghost_node_indices.end())
            {
                norm_vel.push_back(norm_2(velocities[i]));
            }  
        }

        // now check that the velocities scale correctly when the viscosity is area-dependent
        meineke_spring_system.SetAreaBasedViscosity(true);
        
        velocities = meineke_spring_system.rCalculateVelocitiesOfEachNode();
        
        std::vector<double> norm_vel_area;
        
        for(unsigned i=0; i<velocities.size(); i++)
        {
            //check if this is a real cell
            if(ghost_node_indices.find(i)==ghost_node_indices.end())
            {
                norm_vel_area.push_back(norm_2(velocities[i]));
            }  
        }
                
        TS_ASSERT(norm_vel.size() > 0);
        
        // note that d0 and d1 are hardcoded in TissueSimulation::mpMechanicsSystem->rCalculateVelocitiesOfEachNode()  
        for(unsigned i=0; i<norm_vel.size(); i++)
        {
            TS_ASSERT_DELTA(norm_vel_area[i], norm_vel[i]/(0.1 +  1.2*0.9), 1e-3);            
        }
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();    
    }   
        
    void TestAreaBasedVisocityOnAPeriodicMesh() throw (Exception)
    {        
        CancerParameters::Instance()->Reset();

        // Set up the simulation time
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
                
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 2;
        
        double scale_factor = 1.2;
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, scale_factor);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();      
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Tissue<2> tissue(*p_mesh, cells);               
        tissue.SetGhostNodes(ghost_node_indices);

        Meineke2001SpringSystem<2> meineke_spring_system(tissue);
    
        // seems quite difficult to test this on a periodic mesh, so just check the areas 
        // of all the cells are correct 
        tissue.CreateVoronoiTessellation();
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            //check if this is a real cell
            if(ghost_node_indices.find(i)==ghost_node_indices.end())
            {
                double area = tissue.rGetVoronoiTessellation().GetFaceArea(i);
                TS_ASSERT_DELTA(area, sqrt(3)*scale_factor*scale_factor/2, 1e-6);
            }
        }        
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();    
    } 


    void TestSpringConstantsForMutantCells()
    {
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // create a small tissue
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0));
        nodes.push_back(new Node<2>(1, false, 0, 2));
        nodes.push_back(new Node<2>(2, false, 2, 2));
        nodes.push_back(new Node<2>(3, false, 2, 0));
        
        ConformingTetrahedralMesh<2,2> mesh(nodes);
        
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, mesh);
        
        Tissue<2> tissue(mesh, cells);
        
        Meineke2001SpringSystem<2> meineke_spring_system(tissue);
        
        // set cells mutation states
        tissue.rGetCellAtNodeIndex(0).SetMutationState(HEALTHY);
        tissue.rGetCellAtNodeIndex(1).SetMutationState(LABELLED);
        tissue.rGetCellAtNodeIndex(2).SetMutationState(APC_TWO_HIT);
        tissue.rGetCellAtNodeIndex(3).SetMutationState(BETA_CATENIN_ONE_HIT);
        
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(0,1)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(1,2)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(2,3)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(3,0)), 15.0, 1e-10);
        
        meineke_spring_system.SetMutantSprings(true);
        
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(0,1)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(1,2)), 22.5, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(2,3)), 30.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(3,0)), 22.5, 1e-10);
        
        meineke_spring_system.SetMutantSprings(true, 4.0, 3.0);
        
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(0,1)), 15.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(1,2)), 45.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(2,3)), 60.0, 1e-10);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(3,0)), 45.0, 1e-10);
        
        SimulationTime::Destroy();
    }
    

    void TestSpringConstantsForIngeBCatCells()
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        HoneycombMeshGenerator generator(6, 12, 0, true, 1.1);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        // Set up cells 
        std::vector<TissueCell> cells;                      
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, INGE_WNT_SWAT_HYPOTHESIS_TWO, false);
        
        Tissue<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);  
        
        WntGradient::Instance()->SetType(LINEAR);  
        WntGradient::Instance()->SetTissue(crypt);
        
        Meineke2001SpringSystem<2> meineke_spring_system(crypt);
        
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(20,21)), 1.50, 1e-10);
        
        meineke_spring_system.SetBCatSprings(true);
        crypt.CreateVoronoiTessellation();  // normally done in a simulation loop
        
        // Note this is just a crap test to check that you get some dependency on BCat of both cells
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(20,21)), 1.5*8.59312/18.14, 1e-5);
        p_params->SetBetaCatSpringScaler(20/6.0);
        TS_ASSERT_DELTA( norm_2(meineke_spring_system.CalculateForceBetweenNodes(20,21)), 1.5*8.59312/20.0, 1e-5);
        
        SimulationTime::Destroy();
        WntGradient::Destroy(); 
    }


    void TestCalcForcesIn3d() throw (Exception)
    {
        CancerParameters::Instance()->Reset();

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        std::vector<TissueCell> cells;
        TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            cell.SetNodeIndex(i);
            cell.SetBirthTime(-50.0);
            cells.push_back(cell);
        }
        
        Tissue<3> tissue(mesh,cells);
        Meineke2001SpringSystem<3> meineke_spring_system(tissue);

        // Test forces on springs
        unsigned nodeA = 0, nodeB = 1;
        Element<3,3>* p_element = mesh.GetElement(0);
        c_vector<double, 3> force = meineke_spring_system.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(nodeA),p_element->GetNodeGlobalIndex(nodeB));
        for(unsigned i=0; i < 3;i++)
        {
            TS_ASSERT_DELTA(force[i],0.0,1e-6);
        }
        
        //  Test forces on nodes        
        for (unsigned i=0 ; i<1 ; i++)
        {
            std::vector<c_vector<double,3> >& velocities = meineke_spring_system.rCalculateVelocitiesOfEachNode();
            
            for (unsigned j=0; j<4; j++)
            {
                for(unsigned k=0;k<3;k++)
                {
                    TS_ASSERT_DELTA(velocities[j](k),0.0,1e-6);
                }
            }            
        }

        //
        // Scale entire mesh and check that forces are correctly calculated
        //
        CancerParameters *p_params = CancerParameters::Instance();
        double scale_factor = 1.5;
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> old_point = mesh.GetNode(i)->rGetLocation();
            ChastePoint<3> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            new_point.rGetLocation()[1] = scale_factor*old_point[1];
            new_point.rGetLocation()[2] = scale_factor*old_point[2];
            mesh.SetNode(i, new_point, false);            
        }
               
        std::vector<c_vector<double,3> >& new_velocities = meineke_spring_system.rCalculateVelocitiesOfEachNode();
        
        for (unsigned j=0; j<4; j++)
        {
            for(unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(fabs(new_velocities[j](k)),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(scale_factor-1)*sqrt(2),1e-6);
            }
        }
        
        //
        //  Move one node and check that forces are correctly calculated
        //
        ConformingTetrahedralMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        Tissue<3> tissue2(mesh2,cells);
        Meineke2001SpringSystem<3> meineke_spring_system2(tissue2);

        c_vector<double,3> old_point = mesh2.GetNode(0)->rGetLocation();
        ChastePoint<3> new_point;
        new_point.rGetLocation()[0] = 0.0;
        new_point.rGetLocation()[1] = 0.0;
        new_point.rGetLocation()[2] = 0.0;
        mesh2.SetNode(0, new_point, false);   
          
        unsigned nodeA2 = 0, nodeB2 = 1;
        Element<3,3>* p_element2 = mesh2.GetElement(0);
        c_vector<double,3> force2 = meineke_spring_system2.CalculateForceBetweenNodes(p_element2->GetNodeGlobalIndex(nodeA2),p_element2->GetNodeGlobalIndex(nodeB2));
        
        for(unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(fabs(force2[i]),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
        
        new_velocities = meineke_spring_system2.rCalculateVelocitiesOfEachNode();
        
        for(unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(new_velocities[0](i),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
  
        SimulationTime::Destroy();
    } 

    
    void TestArchiving() throw (Exception)
    {   
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "meineke_spring_system.arch";

        unsigned num_nodes;
        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        
            ConformingTetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            num_nodes = mesh.GetNumNodes();

            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            std::vector<TissueCell> cells;
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                cell.SetNodeIndex(i);
                cell.SetBirthTime(-50.0);
                cells.push_back(cell);
            }
        
            Tissue<2> tissue(mesh,cells);
            Meineke2001SpringSystem<2> meineke_spring_system(tissue);
         
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            Meineke2001SpringSystem<2> * const p_meineke_spring_system = &meineke_spring_system;  
            
            p_meineke_spring_system->UseCutoffPoint(1.1);
            p_meineke_spring_system->SetAreaBasedViscosity(true);
            p_meineke_spring_system->SetMutantSprings(true,0.2,0.3);
            p_meineke_spring_system->SetBCatSprings(true);

            output_arch << p_meineke_spring_system;
        }
       
        {
            Tissue<2>::meshPathname = "mesh/test/data/square_2_elements";
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            Meineke2001SpringSystem<2>* p_meineke_spring_system;
            
            // restore from the archive
            input_arch >> p_meineke_spring_system;
            
            // test the member data
            TS_ASSERT_EQUALS(p_meineke_spring_system->mUseCutoffPoint,true);
            TS_ASSERT_DELTA(p_meineke_spring_system->mCutoffPoint,1.1,1e-12);            
            TS_ASSERT_EQUALS(p_meineke_spring_system->mUseEdgeBasedSpringConstant, false);
            TS_ASSERT_EQUALS(p_meineke_spring_system->mUseAreaBasedViscosity, true);
            TS_ASSERT_EQUALS(p_meineke_spring_system->mUseMutantSprings, true);
            TS_ASSERT_DELTA(p_meineke_spring_system->mMutantMutantMultiplier, 0.2, 1e-12);
            TS_ASSERT_DELTA(p_meineke_spring_system->mNormalMutantMultiplier, 0.3, 1e-12);
            TS_ASSERT_EQUALS(p_meineke_spring_system->mUseBCatSprings, true);

        }
    } 
};

#endif /*TESTMEINEKE2001SPRINGSYSTEM_HPP_*/

