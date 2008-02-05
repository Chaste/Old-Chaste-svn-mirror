#ifndef TESTSIMPLETISSUEMECHANICSSYSTEM_HPP_
#define TESTSIMPLETISSUEMECHANICSSYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <vector>

#include "SimpleTissueMechanicsSystem.hpp"
#include "SimpleTissue.cpp"
#include "CellsGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestSimpleTissueMechanicsSystem : public AbstractCancerTestSuite
{    
public:

    void TestForceCalculations() throw (Exception)
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        // Create mesh
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        HoneycombMeshGenerator generator(cells_across, cells_up, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        
        // Get nodes
        std::vector<Node<2> > nodes;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            nodes.push_back(*(p_mesh->GetNode(i)));
        }        
        
        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellMutationState mutation_state = HEALTHY;
            TissueCell cell(STEM, mutation_state, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            
            if (i==10 || i==11)
            {
                cell.SetBirthTime(-0.1);
            }
            else
            {
                cell.SetBirthTime(-10);
            }
            cells.push_back(cell);
        }

        cells[23].StartApoptosis();
        
        // Create simple tissue
        SimpleTissue<2> simple_tissue(nodes, cells);
        
        // Create simple tissue mechanics system (with default cutoff=1.5)
        SimpleTissueMechanicsSystem<2> mechanics_system(simple_tissue);        
        
        // Test rCalculateVelocitiesOfEachNode() method
        std::vector<c_vector<double,2> >& velocities_on_each_node = mechanics_system.rCalculateVelocitiesOfEachNode();
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            // Cells 10 and 11 are young and so the spring between them 
            // temporarily has a smaller natural rest length
            if (i==10)
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], 6.7499, 1e-4);
            }
            else if (i==11)
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], -6.7499, 1e-4);
            }
            else
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], 0.0, 1e-4);
            }
            TS_ASSERT_DELTA(velocities_on_each_node[i][1], 0.0, 1e-4);
        }
               
        // Test rCalculateVelocitiesOfEachNode() further with an increased cutoff
        mechanics_system.SetCutoffPoint(1.9);
        
        // Now each node will experience a non-zero force from any node 
        // a distance sqrt(3) away
        velocities_on_each_node = mechanics_system.rCalculateVelocitiesOfEachNode();
        
        // Node 0 is at (0,0), the bottom left corner, so experiences force 
        // contributions from node 14 at (0,sqrt(3)) and node 8 at (1.5,sqrt(3)/2)
        TS_ASSERT_DELTA(velocities_on_each_node[0][0], 16.4711, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[0][1], 28.5288, 1e-4);
        
        // Node 7 is at (0.5,sqrt(3)/2), on the left, so experiences force 
        // contributions from node 2 at (2,0), node 16 at (2,sqrt(3))
        // and node 21 at (0.5,1.5*sqrt(3))
        TS_ASSERT_DELTA(velocities_on_each_node[7][0], 32.9422, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[7][1], 19.0192  , 1e-4);
        
        // Node 8 is at (1.5,sqrt(3)/2), on the left, so experiences force
        // contributions from node 0 at (0,0), node 14 at (0,sqrt(3)) and 
        // node 22 at (1.5,1.5*sqrt(3)), node 17 at (3,sqrt(3)) and node 3 at (3,0)
        TS_ASSERT_DELTA(velocities_on_each_node[8][0], 0.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[8][1], 19.0192, 1e-4);
    }

};

#endif /*TESTSIMPLETISSUEMECHANICSSYSTEM_HPP_*/
