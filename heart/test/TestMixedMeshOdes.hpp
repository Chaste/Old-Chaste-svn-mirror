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


#ifndef TESTMIXEDMESHODES_HPP_
#define TESTMIXEDMESHODES_HPP_

#include <cxxtest/TestSuite.h>

#include "MixedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "DistributedVector.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "FastSlowLuoRudyIModel1991.hpp"
#include "ZeroStimulus.hpp"
#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"


class TestMixedMeshOdes : public CxxTest::TestSuite
{
public:

    /*
     * Set up a mixed mesh
     * Set up a vector of cells associated with the fine nodes of the mesh
     * Run the ODE's on the coarse-mesh cells
     * Interpolate the slow currents from the coarse-mesh cells to the fine-mesh cells
     * Run the ODE's forward on the fine-mesh cells
     * 
     * The mesh is 3 by 3 square, with the coarse mesh just the 4 corner nodes. 
     * 
     */
    void Test2DSerial(void)
    {
        ConformingTetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRectangularMesh(2, 2, false);
        double half=1.0L/2.0L;
        fine_mesh.Scale(half, half, 0.0);
        
        // create coarse mesh as RTM
        MixedTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(1, 1, false);
        
        // give fine mesh to coarse mesh and calculate node map
        coarse_mesh.SetFineMesh(&fine_mesh);
        
        unsigned num_fine_nodes=fine_mesh.GetNumNodes();
        
        // create a vector of bools telling us which of the fine nodes
        // have a coarse counter part
        std::vector<bool> fine_has_coarse_counterpart;
        fine_has_coarse_counterpart.resize(num_fine_nodes);
        for (unsigned fine_node_index=0; fine_node_index<num_fine_nodes; fine_node_index++)
        {
            fine_has_coarse_counterpart[fine_node_index] = false;
        }

        for (unsigned coarse_node_index=0;
             coarse_node_index<coarse_mesh.GetNumNodes();
             coarse_node_index++)
        {
            unsigned fine_node_index = coarse_mesh.rGetCoarseFineNodeMap().GetNewIndex(coarse_node_index);
            fine_has_coarse_counterpart[fine_node_index] = true;
        }

        // print the mesh
        /*
        std::cout << "Coarse\n";
        for(unsigned i=0; i<coarse_mesh.GetNumNodes(); i++)
        {
            std::cout << i << " -- " << coarse_mesh.GetNode(i)->rGetLocation() << "\n";
        }
        std::cout << "Fine\n";
        for(unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            std::cout << i << " -- " << fine_mesh.GetNode(i)->rGetLocation() << "\n";
        }
        */
        
        // create vector of cells - slow at the coarse nodes and fast at the fine-and-not-coarse nodes    
        std::vector<FastSlowLuoRudyIModel1991* > cells;
        cells.resize(num_fine_nodes);
        
        ZeroStimulus zero_stimulus;
        InitialStimulus full_stimulus(-600, 0.5);
        InitialStimulus half_stimulus(-300, 0.5);
        EulerIvpOdeSolver solver;
        double time_step = 0.01;
        double start_time = 0.0;
        double end_time = 0.1;
        
        for (unsigned fine_node_index=0; fine_node_index<num_fine_nodes; fine_node_index++)
        {
            double x_location = fine_mesh.GetNode(fine_node_index)->rGetLocation()[0];
            
            AbstractStimulusFunction* p_stimulus;
                      
            if(x_location == 0.0)
            {
                p_stimulus = &zero_stimulus;
            }
            else if(x_location == 0.5)
            {
                p_stimulus = &half_stimulus;             
            }
            else //xlocation == 1.0
            {
                p_stimulus = &full_stimulus;
            }
            
            // fast_model if fine has no coarse counterpart, slow if it does
            // (ie fast if fine only node, slow if coarse node)
            bool is_fast_model = !fine_has_coarse_counterpart[fine_node_index];
            
            cells[fine_node_index] = new FastSlowLuoRudyIModel1991(is_fast_model,
                                                                        &solver,
                                                                        time_step,
                                                                        p_stimulus);
        }
        
        // run odes on coarse mesh cells
        for (unsigned fine_node_index=0; fine_node_index<num_fine_nodes; fine_node_index++)
        {
            if (fine_has_coarse_counterpart[fine_node_index])
            {
                assert(cells[fine_node_index]->IsFast()==false);
                cells[fine_node_index]->Compute(start_time, end_time);
            }
            
        }
        
        // interpolate slow currents to fine mesh cells
        for (unsigned index=0; index<num_fine_nodes; index++)
        {
            if (!fine_has_coarse_counterpart[index])
            {
                Element<2,2>* p_coarse_element = coarse_mesh.GetACoarseElementForFineNodeIndex(index);
                
                const ChastePoint<2>& r_position_of_fine_node = fine_mesh.GetNode(index)->rGetLocation();

                c_vector<double,2+1> weights = p_coarse_element->CalculateInterpolationWeights(r_position_of_fine_node);
         
                unsigned num_slow_values = cells[p_coarse_element->GetNodeGlobalIndex(0)]->GetNumSlowValues();
                std::vector<double> interpolated_slow_values(num_slow_values, 0.0); 
                for (unsigned i=0; i<2+1/*num_nodes*/; i++)
                {
                    unsigned coarse_cell_index = p_coarse_element->GetNodeGlobalIndex(i);
                    unsigned corresponding_fine_mesh_index = coarse_mesh.rGetCoarseFineNodeMap().GetNewIndex(coarse_cell_index);
                    FastSlowLuoRudyIModel1991* p_coarse_node_cell = cells [ corresponding_fine_mesh_index ];
                    
                    TS_ASSERT_EQUALS(fine_has_coarse_counterpart [corresponding_fine_mesh_index], true); 
                    TS_ASSERT_EQUALS(p_coarse_node_cell->IsFast(), false);
                    
                    std::vector<double> nodal_slow_values;
                    p_coarse_node_cell->GetSlowValues(nodal_slow_values);
                    assert(nodal_slow_values.size() == num_slow_values);
                    for(unsigned j=0; j<nodal_slow_values.size(); j++)
                    {                    
                        interpolated_slow_values[j] += nodal_slow_values[j]*weights(i);
                    }
                }

                cells[index]->SetSlowValues(interpolated_slow_values);
            }
        }        
        
        // solve ODEs on fine-only nodes
        for (unsigned fine_node_index=0; fine_node_index<num_fine_nodes; fine_node_index++)
        {
            if (!fine_has_coarse_counterpart[fine_node_index])
            {
                assert(cells[fine_node_index]->IsFast());
                cells[fine_node_index]->Compute(start_time, end_time);
            }            
        }
        
        // destroy vector of cells
        for (unsigned fine_node_index=0; fine_node_index<num_fine_nodes; fine_node_index++)
        {
            delete cells[fine_node_index];
        }
    }
};


#endif /*TESTMIXEDMESHODES_HPP_*/
