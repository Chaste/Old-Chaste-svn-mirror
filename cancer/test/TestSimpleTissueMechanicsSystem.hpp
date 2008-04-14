#ifndef TESTSIMPLETISSUEMECHANICSSYSTEM_HPP_
#define TESTSIMPLETISSUEMECHANICSSYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <vector>

#include "SimpleTissueMechanicsSystem.hpp"
#include "SimpleTissue.hpp"
#include "CellsGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "CellsGenerator.hpp"

class TestSimpleTissueMechanicsSystem : public AbstractCancerTestSuite
{    
public:
    void TestForceCalculationsBasic() throw (Exception)
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        double stiffness = CancerParameters::Instance()->GetSpringStiffness();
        
        // Create a mesh and get its nodes
        HoneycombMeshGenerator generator(3, 3, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<Node<2> > nodes;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            nodes.push_back(*(p_mesh->GetNode(i)));
        }        
        
        // Create cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells,*p_mesh);

        // Create simple tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // Create simple tissue mechanics system (with default cutoff=1.5)
        SimpleTissueMechanicsSystem<2> mechanics_system(simple_tissue);
        
        // Test rCalculateVelocitiesOfEachNode() method
        std::vector<c_vector<double,2> >& velocities_on_each_node = mechanics_system.rCalculateVelocitiesOfEachNode();
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(velocities_on_each_node[i][0], 0.0, 1e-4);
            TS_ASSERT_DELTA(velocities_on_each_node[i][1], 0.0, 1e-4);
        }
               
        // Test rCalculateVelocitiesOfEachNode() further with an increased cutoff
        mechanics_system.UseCutoffPoint(1.9);
        
        // Now each node will experience a non-zero force from any node 
        // a distance sqrt(3) away
        velocities_on_each_node = mechanics_system.rCalculateVelocitiesOfEachNode();
        
        double force = stiffness*(sqrt(3)-1)*exp(-5.0*(sqrt(3)-1)); //magnitude of force between any two nodes

        double inv_damping = 1.0/CancerParameters::Instance()->GetDampingConstantNormal();

        // Node 0 is at (0,0), has force contributions from two nodes (nodes 4 and 6)
        // that are sqrt(3) away.
        // All nodes are connected to two other nodes, sqrt(3) away
        // Total force on a node is force*sqrt(3) in magnitude
        TS_ASSERT_DELTA(velocities_on_each_node[0][0], inv_damping*force*sqrt(3)/2.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[0][1], inv_damping*force*3.0/2.0, 1e-4);
        
        TS_ASSERT_DELTA(velocities_on_each_node[1][0], inv_damping*force*sqrt(3)/2.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[1][1], inv_damping*force*3.0/2.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[2][0], -inv_damping*force*sqrt(3)/2.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[2][1],  inv_damping*force*3.0/2.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[3][0], inv_damping*force*sqrt(3), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[3][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[4][0], -inv_damping*force*sqrt(3), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[4][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[5][0], -inv_damping*force*sqrt(3), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[5][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[6][0], inv_damping*force*sqrt(3)/2.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[6][1], -inv_damping*force*3.0/2.0, 1e-4);
        
        TS_ASSERT_DELTA(velocities_on_each_node[7][0], inv_damping*force*sqrt(3)/2.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[7][1], -inv_damping*force*3.0/2.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[8][0], -inv_damping*force*sqrt(3)/2.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[8][1], -inv_damping*force*3.0/2.0, 1e-4);
    }


    void TestForceCalculationsWithBirthAndDeath() throw (Exception)
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,10);
        
        // Create mesh
        unsigned cells_across = 2;
        unsigned cells_up = 2;
        HoneycombMeshGenerator generator(cells_across, cells_up, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        
        // Get nodes
        std::vector<Node<2> > nodes;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            nodes.push_back(*(p_mesh->GetNode(i)));
        }        
        
        // Create cells, making cells 0 and 1 young
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellMutationState mutation_state = HEALTHY;
            TissueCell cell(STEM, mutation_state, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            
            if (i==0 || i==1)
            {
                cell.SetBirthTime(-0.1);
            }
            else
            {
                cell.SetBirthTime(-10);
            }
            cells.push_back(cell);
        }

       
        // Create simple tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // Create simple tissue mechanics system (with default cutoff=1.5)
        SimpleTissueMechanicsSystem<2> mechanics_system(simple_tissue);
        
        // Test rCalculateVelocitiesOfEachNode() method
        
        /// \todo Improve this test! (see #642)
        std::vector<c_vector<double,2> >& velocities_on_each_node = mechanics_system.rCalculateVelocitiesOfEachNode();
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            // Cells 0 and 1 are young and so the spring between them 
            // temporarily has a smaller natural rest length
            if (i==0)
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], 0.7114, 1e-4);
            }
            else if (i==1)
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], -0.7114, 1e-4);
            }
            else
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], 0.0, 1e-4);
            }
            TS_ASSERT_DELTA(velocities_on_each_node[i][1], 0.0, 1e-4);
        }

        // Get cells 0 and 3 and start apoptosis on them          
        SimpleTissue<2>::Iterator iter = simple_tissue.Begin();
        iter->StartApoptosis(); //cell 0
        ++iter;
        ++iter;
        ++iter;
        iter->StartApoptosis(); //cell 3
        SimulationTime::Instance()->IncrementTimeOneStep();
        
        velocities_on_each_node = mechanics_system.rCalculateVelocitiesOfEachNode();

        // The velocities of all cells are affected by the fact that cells 0 and 3 are apoptosing.
        TS_ASSERT_DELTA(velocities_on_each_node[0][0], 1.1311, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[0][1], 0.9557, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[1][0], -0.0275, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[1][1], 0.9557, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[2][0], 0.5518, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[2][1], -0.9557, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[3][0], -1.6554, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[3][1], -0.9557, 1e-4);    
    }
    
    void TestArchiving() throw (Exception)
    {   
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "simple_mech_system.arch";

        {
            // Set up SimulationTime
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
            
            // Create mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");        
            ConformingTetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            
            // Get nodes
            std::vector<Node<2> > nodes;
            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                nodes.push_back(*(mesh.GetNode(i)));
            }

            std::vector<TissueCell> cells;
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                cell.SetNodeIndex(i);
                cell.SetBirthTime(-50.0);
                cells.push_back(cell);                
            }
        
            // Create simple tissue
            SimpleTissue<2> simple_tissue(nodes, cells);
            
            // Create simple tissue mechanics system (with default cutoff=1.5)
            SimpleTissueMechanicsSystem<2> mechanics_system(simple_tissue);
                     
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            SimpleTissueMechanicsSystem<2> * const p_mech_system = &mechanics_system;  
            
            p_mech_system->UseCutoffPoint(1.1);
            
            output_arch << p_mech_system;
        }
       
        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            SimpleTissueMechanicsSystem<2>* p_mech_system;
            
            // Restore from the archive
            input_arch >> p_mech_system;
            
            // Test the member data
            TS_ASSERT_DELTA(p_mech_system->mCutoffPoint,1.1,1e-12);
            
            delete p_mech_system->mpTissue;
            delete p_mech_system;
        }
    } 
    
};

#endif /*TESTSIMPLETISSUEMECHANICSSYSTEM_HPP_*/
