#ifndef TESTTISSUEFORCEEXPERIMENTS_HPP_
#define TESTTISSUEFORCEEXPERIMENTS_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"

//////////////////////////////////////////////////////////////////////////////////////
// This is Pras experimental work.. - to be moved to the appropriate place when the 
// repository restucture is done..
//
//  - presumably if the TissueSimulation class is altered and keeping this working 
//    would a problem (as a class is defined below that inherits from 
//    TissueSimulation), this test and the class should just be commented out, or
//    something
//
//////////////////////////////////////////////////////////////////////////////////////
class TissueSimulationForForceExperiments : public TissueSimulation<2>
{
private :
    bool mFixXNotY;
    double mMaxHeight;
    
public :
    TissueSimulationForForceExperiments(Tissue<2>& rTissue,
                                        bool fixXNotY)
        : TissueSimulation<2>(rTissue),
          mFixXNotY(fixXNotY)
    {   
        // this class is hardcoded for a particular honeycomb mesh! check num nodes is 
        // as expected  
        assert(mrTissue.rGetMesh().GetNumNodes()==360);
        
        // calc max value in the fixed direction.
        unsigned x_or_y = mFixXNotY ? 0 : 1;

        double max = 0;
        for(Tissue<2>::Iterator iter = rGetTissue().Begin();
            iter != rGetTissue().End();
            ++iter)
        {
            double val = iter.rGetLocation()[x_or_y];
            if(val > max)
            {
                max = val;
            }
        }

        mMaxHeight = max;
    }
    
    virtual ~TissueSimulationForForceExperiments()
    {
    }

//////////
// doesn't do anything useful anymore as tissue sim doesn't loop over ghost edges...
//
// could be removed
//
// need to overload appropriate method in Tissue...
//////////    
//    /*
//     *  Overloaded CalcForceBetweenNodes
//     *  
//     *  only a minor change to the method - the rest length between a ghost node and a real node
//     *  is set to be 2 - this pushes ghost nodes away and stops them interfering this real nodes
//     *  and possibly breaking bonds
//     * 
//     *  Note this means the initial honeycomb mesh isn't the equilibrium solution if stretch=1,
//     *  so not not much point running if stretch=1.
//     */
//    c_vector<double,2> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex)
//    {
//        assert(nodeAGlobalIndex!=nodeBGlobalIndex);
//        c_vector<double,2> unit_difference;
//        c_vector<double,2> node_a_location = mrTissue.rGetMesh().GetNode(nodeAGlobalIndex)->rGetLocation();
//        c_vector<double,2> node_b_location = mrTissue.rGetMesh().GetNode(nodeBGlobalIndex)->rGetLocation();
//        
//        // there is reason not to substract one position from the other (cyclidrical meshes). clever gary
//        unit_difference = mrTissue.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);   
//        double distance_between_nodes = norm_2(unit_difference);
//        unit_difference /= distance_between_nodes;
//        
//        
//        double rest_length = 1.0;
//
//        if(!mrTissue.rGetGhostNodes()[nodeAGlobalIndex] && mrTissue.rGetGhostNodes()[nodeBGlobalIndex])
//        {
//            assert(0);
//            rest_length = 2;
//        }
//        if(mrTissue.rGetGhostNodes()[nodeAGlobalIndex] && !mrTissue.rGetGhostNodes()[nodeBGlobalIndex])
//        {
//            assert(0);
//            rest_length = 2;
//        }
//
//        return mpParams->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
//    }


    /*
     *  Overloaded method.
     *   
     *  Doesn't allow nodes to move into
     *    - y<0 or y>max if pushing/pulling top/bottom
     *    - x<0 or x>max if pushing/pulling sides
     * 
     *  If pushing comment out appropriate 2 lines below if want to allow slip
     */
    void UpdateNodePositions(const std::vector<c_vector<double,2> >& rDrDt)
    {
        mrTissue.UpdateGhostPositions(mDt);

        for (Tissue<2>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            TissueCell& cell = *cell_iter;
            unsigned index = cell.GetNodeIndex();
            
            ChastePoint<2> new_point(mrTissue.rGetMesh().GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
            
            if(mFixXNotY)
            { 
                // all nodes must stay in x>0
                if(new_point.rGetLocation()[0] < 0)
                {
                    new_point.rGetLocation()[0] = 0;
                }

                // all nodes must stay in x<max height
                if(new_point.rGetLocation()[0] > mMaxHeight)
                {
                    new_point.rGetLocation()[0] = mMaxHeight;
                }

                // right hand side
                if(index>=57 && index<=327 && ((index-27)%30==0) )
                {
                    new_point.rGetLocation()[0] = mMaxHeight;
                    new_point.rGetLocation()[1] = mrTissue.rGetMesh().GetNode(index)->rGetLocation()[1];
                }
              
                // left hand side
                if(index>=32 && index<=302 && ((index-2)%30==0) )
                {
                    new_point.rGetLocation()[0] = 0;
                    new_point.rGetLocation()[1] = mrTissue.rGetMesh().GetNode(index)->rGetLocation()[1];
                }
            }
            else
            {
                // all nodes must stay in y>0
                if(new_point.rGetLocation()[1] < 0)
                {
                    new_point.rGetLocation()[1] = 0;
                }                
                // all nodes must stay in y<max height
                if(new_point.rGetLocation()[1] > mMaxHeight)
                {
                    new_point.rGetLocation()[1] = mMaxHeight;
                }

                // top surface
                if(index>=317 && index<=327)
                {
                    // comment this first line out to allow slip on top surface
                    new_point.rGetLocation()[0] = mrTissue.rGetMesh().GetNode(index)->rGetLocation()[0];
                    new_point.rGetLocation()[1] = mMaxHeight;
                }
                
                // bottom surface
                if(index>=32 && index<=42)
                {
                    // comment this first line out to allow slip on bottom surface
                    new_point.rGetLocation()[0] = mrTissue.rGetMesh().GetNode(index)->rGetLocation()[0];
                    new_point.rGetLocation()[1] = 0;
                }
            }
            
            mrTissue.MoveCell(cell_iter, new_point);
        }
    }


    /*
     *  Calculate the force on the pulled/pushed boundary
     *  
     *  Estimated by looping over all cells on that boundary, and summing the absolute value of
     *  the force in each spring connected those cells
     */
    c_vector<double,2> CalculateTotalForce()
    {
        c_vector<double,2> total_force = zero_vector<double>(2);

        std::set<std::set<unsigned> > node_pairs_checked;
        for (unsigned elem_index = 0; elem_index<mrTissue.rGetMesh().GetNumAllElements(); elem_index++)
        {
            Element<2,2>* p_element = mrTissue.rGetMesh().GetElement(elem_index);
            if (!p_element->IsDeleted())
            {
                for (unsigned k=0; k<2+1; k++)
                {
                    unsigned nodeA = k;
                    
                    for(unsigned l=k+1; l<k+2+1; l++)
                    {
                        unsigned nodeB = l%(2+1);
                    
                        assert(!p_element->GetNode(nodeA)->IsDeleted());
                        assert(!p_element->GetNode(nodeB)->IsDeleted());
                                        
                        unsigned nodeA_global_index = p_element->GetNode(nodeA)->GetIndex();
                        unsigned nodeB_global_index = p_element->GetNode(nodeB)->GetIndex();
                        
                        if(    mrTissue.rGetGhostNodes()[nodeA_global_index] 
                            || mrTissue.rGetGhostNodes()[nodeB_global_index] )
                        {
                            break;
                        }
                        
                        unsigned index = mFixXNotY ? 0 : 1;

                        double val_at_A = p_element->GetNode(nodeA)->rGetLocation()[index];
                        double val_at_B = p_element->GetNode(nodeB)->rGetLocation()[index];
    
                        assert(val_at_A <= mMaxHeight );
                        assert(val_at_B <= mMaxHeight );
    
                        if( (fabs(val_at_A-mMaxHeight) < 1e-6) || (fabs(val_at_B-mMaxHeight) < 1e-6) )
                        {
                            // check whether we have already worked out the force between these two...
                            bool is_force_already_calculated = false;
                            
                            std::set<unsigned> current_node_pair;
                            current_node_pair.insert(nodeA_global_index);
                            current_node_pair.insert(nodeB_global_index);
                            
                            // see if the node pair is in the set of node pairs done
                            std::set<std::set<unsigned> >::iterator set_iter = node_pairs_checked.find(current_node_pair);                    
                            if(set_iter!=node_pairs_checked.end())
                            {
                                // node pair found
                                is_force_already_calculated = true;
                            }
                            else
                            {
                                is_force_already_calculated = false;
                                // add the node pair to the list of node pairs
                                node_pairs_checked.insert(current_node_pair);
                            }
                            
                            if(!is_force_already_calculated)
                            {
                                c_vector<double,2> force = CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(nodeA),p_element->GetNodeGlobalIndex(nodeB));
                                if ((!mrTissue.rGetGhostNodes()[nodeA_global_index]) && (!mrTissue.rGetGhostNodes()[nodeB_global_index]))
                                {
                                    for(unsigned i=0; i<2; i++)
                                    {
                                        total_force(i) += fabs(force[i]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return total_force;
    }
};

class TissueSimulationForForceExperimentsShearing : public TissueSimulation<2>
{
private :
    double mShear;
    double mMaxHeight;
    bool mFirst;
    
public :
    TissueSimulationForForceExperimentsShearing(Tissue<2>& rTissue, double shear)
        : TissueSimulation<2>(rTissue),
          mShear(shear)
    {   
        mFirst = true;
        // this class is hardcoded for a particular honeycomb mesh! check num nodes is 
        // as expected  
        assert(mrTissue.rGetMesh().GetNumNodes()==360);
        
        double max = 0;
        for(Tissue<2>::Iterator iter = rGetTissue().Begin();
            iter != rGetTissue().End();
            ++iter)
        {
            double val = iter.rGetLocation()[1];
            if(val > max)
            {
                max = val;
            }
        }

        mMaxHeight = max;
    }
    
    virtual ~TissueSimulationForForceExperimentsShearing()
    {
    }

    /*
     *  Overloaded method.
     */
    void UpdateNodePositions(const std::vector<c_vector<double,2> >& rDrDt)
    {
        double disp = 0.0;
        if(mFirst)
        {
            disp = mShear;
            mFirst = false;
        }
        
        mrTissue.UpdateGhostPositions(mDt);

        for (Tissue<2>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            TissueCell& cell = *cell_iter;
            unsigned index = cell.GetNodeIndex();
            
            ChastePoint<2> new_point(mrTissue.rGetMesh().GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
            
            // all nodes must stay in y>0
            if(new_point.rGetLocation()[1] < 0)
            {
                new_point.rGetLocation()[1] = 0;
            }                
            // all nodes must stay in y<max height
            if(new_point.rGetLocation()[1] > mMaxHeight)
            {
                new_point.rGetLocation()[1] = mMaxHeight;
            }

            // top surface
//            if(index>=317 && index<=327)

            if(index>=330 && index<=344)
            {
                // comment this first line out to allow slip on top surface
                new_point.rGetLocation()[0] = mrTissue.rGetMesh().GetNode(index)->rGetLocation()[0] + disp;
                new_point.rGetLocation()[1] = mMaxHeight;
            }
            
            // bottom surface
//            if(index>=32 && index<=42)

            if(index>=30 && index<=44)
            {
                // comment this first line out to allow slip on bottom surface
                new_point.rGetLocation()[0] = mrTissue.rGetMesh().GetNode(index)->rGetLocation()[0];
                new_point.rGetLocation()[1] = 0;
            }
            
            mrTissue.MoveCell(cell_iter, new_point);
        }
    }


    /*
     *  Calculate the force on the pulled/pushed boundary
     *  
     *  Estimated by looping over all cells on that boundary, and summing the absolute value of
     *  the force in each spring connected those cells
     */
    c_vector<double,2> CalculateTotalForce()
    {
        c_vector<double,2> total_force = zero_vector<double>(2);

        std::set<std::set<unsigned> > node_pairs_checked;
        for (unsigned elem_index = 0; elem_index<mrTissue.rGetMesh().GetNumAllElements(); elem_index++)
        {
            Element<2,2>* p_element = mrTissue.rGetMesh().GetElement(elem_index);
            if (!p_element->IsDeleted())
            {
                for (unsigned k=0; k<2+1; k++)
                {
                    unsigned nodeA = k;
                    
                    for(unsigned l=k+1; l<k+2+1; l++)
                    {
                        unsigned nodeB = l%(2+1);
                    
                        assert(!p_element->GetNode(nodeA)->IsDeleted());
                        assert(!p_element->GetNode(nodeB)->IsDeleted());
                                        
                        unsigned nodeA_global_index = p_element->GetNode(nodeA)->GetIndex();
                        unsigned nodeB_global_index = p_element->GetNode(nodeB)->GetIndex();
                        
                        if(    mrTissue.rGetGhostNodes()[nodeA_global_index] 
                            || mrTissue.rGetGhostNodes()[nodeB_global_index] )
                        {
                            break;
                        }
                        
                        unsigned index = 1;

                        double val_at_A = p_element->GetNode(nodeA)->rGetLocation()[index];
                        double val_at_B = p_element->GetNode(nodeB)->rGetLocation()[index];
    
                        assert(val_at_A <= mMaxHeight );
                        assert(val_at_B <= mMaxHeight );
    
                        if( (fabs(val_at_A-mMaxHeight) < 1e-6) || (fabs(val_at_B-mMaxHeight) < 1e-6) )
                        {
                            // check whether we have already worked out the force between these two...
                            bool is_force_already_calculated = false;
                            
                            std::set<unsigned> current_node_pair;
                            current_node_pair.insert(nodeA_global_index);
                            current_node_pair.insert(nodeB_global_index);
                            
                            // see if the node pair is in the set of node pairs done
                            std::set<std::set<unsigned> >::iterator set_iter = node_pairs_checked.find(current_node_pair);                    
                            if(set_iter!=node_pairs_checked.end())
                            {
                                // node pair found
                                is_force_already_calculated = true;
                            }
                            else
                            {
                                is_force_already_calculated = false;
                                // add the node pair to the list of node pairs
                                node_pairs_checked.insert(current_node_pair);
                            }
                            
                            if(!is_force_already_calculated)
                            {
                                c_vector<double,2> force = CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(nodeA),p_element->GetNodeGlobalIndex(nodeB));
                                if ((!mrTissue.rGetGhostNodes()[nodeA_global_index]) && (!mrTissue.rGetGhostNodes()[nodeB_global_index]))
                                {
                                    for(unsigned i=0; i<2; i++)
                                    {
                                        total_force(i) += fabs(force[i]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return total_force;
    }
};


class TestTissueForceExperiments : public CxxTest::TestSuite
{
private:
    void ApplyRandomDisplacements(Tissue<2>& rTissue)
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = rTissue.rGetMesh();
        
        double max_disp = 0.1;
        
        for(unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            double random_num = RandomNumberGenerator::Instance()->ranf();
            double random_distance = max_disp*random_num;

            random_num = RandomNumberGenerator::Instance()->ranf();
            double random_angle = 2.0*M_PI*random_num;
        
            double u = random_distance*cos(random_angle);
            double v = random_distance*sin(random_angle);
            
            r_mesh.GetNode(i)->rGetModifiableLocation()[0] += u;
            r_mesh.GetNode(i)->rGetModifiableLocation()[1] += v;
        }
        
        rTissue.ReMesh();
    }            


public:
    /*
     *  To run: choose the required stretches and set fix_X_not_Y as required, and set
     *  the run time
     * 
     *   - if squashing: fix in Y direction, and perhaps edit the simulation class above to
     *     allow slip
     *   - if stretching: fix in X direction.
     * 
     *  printed values are the stretch and the force in the relevent direction, and also
     *  the width across the top of the mesh (only relevent if squashing, fixing in Y direction
     *  and allowing slip so that the width can change.
     * 
     */
    void TestMeinekeIncremental() throw(Exception)
    {      
        CancerParameters::Instance()->Reset();

        std::vector<double> stretches;
        std::vector<double> forces;
        std::vector<double> areas;

        // To get this test to pass use a stretch of 1.1, fixXNotY = true 
        // and end time of 0.3
        for(unsigned i=1; i<=12; i++)
        {
            stretches.push_back(1 - 0.05*i);
        }

        for(unsigned i=11; i>=1; i--)
        {
            stretches.push_back(1 - 0.05*i);
        }

        
        bool fix_X_not_Y = false;
        double run_time = 3;


        ////////////////////////////////////////////////////        
        // the main code 
        ////////////////////////////////////////////////////
        int num_cells_depth = 20; // the TissueSimulationForForceExperiments class expects these values!
        int num_cells_width = 11; // the TissueSimulationForForceExperiments class expects these values!
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();        
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        double old_stretch = 1.0;
        for(unsigned i=0; i<stretches.size(); i++)
        {
            double stretch = stretches[i];

            if(fix_X_not_Y)
            {
                p_mesh->Scale(stretch/old_stretch, 1.0);
            }
            else
            {
    	        p_mesh->Scale(1.0, stretch/old_stretch);
            }

            std::stringstream string_stream;
            
            if(fix_X_not_Y)
            {
       	        string_stream << "TissueForceExperimentMeinekeIncrementalX_" << stretch;
            }
            else
            {
        	    string_stream << "TissueForceExperimentMeinekeIncrementalY_" << stretch;
            }
    
            std::string output_directory = string_stream.str();
	        
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
		
            std::vector<TissueCell> cells;
	
            for(unsigned j=0; j<p_mesh->GetNumNodes(); j++)
            {
                TissueCell cell(DIFFERENTIATED, HEALTHY, 0, new FixedCellCycleModel());
                double birth_time = -10;
                cell.SetNodeIndex(j);
	            cell.SetBirthTime(birth_time);
                cells.push_back(cell);
            }

            Tissue<2> tissue(*p_mesh,cells);
            tissue.SetGhostNodes(ghost_node_indices);

            // apply a random displacement to each node. doesn't fix the
            // meineke remeshing of a stretched mesh issue..
            //ApplyRandomDisplacements(tissue);
            
	        TissueSimulationForForceExperiments simulator(tissue, fix_X_not_Y);

            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(run_time);
            tissue.ReMesh();

            simulator.Solve();
            double width = p_mesh->GetNode(327)->rGetLocation()[0] - p_mesh->GetNode(317)->rGetLocation()[0];	
            c_vector<double,2> force = simulator.CalculateTotalForce();
        
            if(fix_X_not_Y)
            {
                forces.push_back( force[0] );
            }
            else
            {
                forces.push_back( force[1] );
            }
            
	        areas.push_back( width );

            tissue.ReMesh();

            std::vector<bool> is_ghost_node = tissue.rGetGhostNodes();
            ghost_node_indices.clear();
            for(unsigned j=0; j<is_ghost_node.size(); j++)
            {
                if(is_ghost_node[j])
                {
                    ghost_node_indices.insert(j);
                }
            }
            old_stretch = stretch;

            SimulationTime::Destroy();
	        RandomNumberGenerator::Destroy();
        }

        std::cout << "\n\nResults - incremental:\n";
        for(unsigned i=0; i<stretches.size(); i++)
        {
            std::cout << stretches[i] << " " << forces[i] << " " << areas[i] << "\n";
        }
        
        // just check nothing has changed. To get this test to pass use
        // a stretch of 1.1, fixXNotY = true and end time of 0.3
        TS_ASSERT_DELTA(forces[0], 28.4007, 1e-3); 
    }

    
    
    
    void xTestMeinekeIncrementalShearing() throw(Exception)
    {            
        CancerParameters::Instance()->Reset();

        std::vector<double> shears;
        std::vector<double> forces;
        std::vector<double> areas;

        //for the test to pass
//        double disp = 0.2;
//        for(unsigned i=1; i<2; i++)
//        {
//            shears.push_back(disp*i);
//        }
//        double run_time = 0.3;

        double disp = 0.2;
        for(unsigned i=1; i<60; i++)
        {
            shears.push_back(disp*i);
        }
        double run_time = 2; // might need to be more..


        ////////////////////////////////////////////////////        
        // the main code 
        ////////////////////////////////////////////////////
//        int num_cells_depth = 20; // the TissueSimulationForForceExperiments class expects these values!
//        int num_cells_width = 11; // the TissueSimulationForForceExperiments class expects these values!
        
//        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
//        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();        


        int num_cells_depth = 20; // the TissueSimulationForForceExperiments class expects these values!
        int num_cells_width = 15; // the TissueSimulationForForceExperiments class expects these values!
 
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
 
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        //double old_stretch = 1.0;
        for(unsigned i=0; i<shears.size(); i++)
        {
            double shear = shears[i];

            std::stringstream string_stream;
            string_stream << "TissueForceExperimentMeinekeShear/" << shear;
            std::string output_directory = string_stream.str();
            
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
        
            std::vector<TissueCell> cells;
    
            for(unsigned j=0; j<p_mesh->GetNumNodes(); j++)
            {
                TissueCell cell(DIFFERENTIATED, HEALTHY, 0, new FixedCellCycleModel());
                double birth_time = -10;
                cell.SetNodeIndex(j);
                cell.SetBirthTime(birth_time);
                cells.push_back(cell);
            }
        
            Tissue<2> tissue(*p_mesh,cells);
            tissue.SetGhostNodes(ghost_node_indices);

            TissueSimulationForForceExperimentsShearing simulator(tissue, disp);

            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(run_time);
            tissue.ReMesh();

            simulator.Solve();

            c_vector<double,2> force = simulator.CalculateTotalForce();
        
            forces.push_back( force[1] );
            
            std::cout << " -> " << shear << " " << force[1] << "\n";
 
            tissue.ReMesh();

            std::vector<bool> is_ghost_node = tissue.rGetGhostNodes();
            ghost_node_indices.clear();
            for(unsigned j=0; j<is_ghost_node.size(); j++)
            {
                if(is_ghost_node[j])
                {
                    ghost_node_indices.insert(j);
                }
            }
 

            SimulationTime::Destroy();
            RandomNumberGenerator::Destroy();
        }

        std::cout << "\n\nResults - shearing:\n";
        for(unsigned i=0; i<shears.size(); i++)
        {
            std::cout << shears[i] << " " << forces[i] << "\n";
        }
        
        // just check nothing has changed. To get this test to pass use
        // disp = shears[0] = 0.2, run_time = 0.3
        TS_ASSERT_DELTA(forces[0], 9.5612, 1e-3); 
    }
};

#endif /*TESTTISSUEFORCEEXPERIMENTS_HPP_*/
