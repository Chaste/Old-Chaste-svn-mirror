#ifndef CRYPTSIMULATION2D_HPP_
#define CRYPTSIMULATION2D_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataWriter.hpp"
#include <cmath>
#include <ctime>
#include <iostream>
#include "TrianglesMeshWriter.cpp"
#include "Exception.hpp"
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "SimulationTime.hpp"
#include "StochasticCellCycleModel.hpp"

#include "ColumnDataWriter.hpp"

#include "MeinekeCryptCellTypes.hpp"
//!!!!!!!!!!!!!!!!!!!!!!!!
/**
 * Solve a 2D crypt simulation based on the Meineke paper.
 *
 * The spring lengths are governed by the equations
 * dr/dt = stem_cycle_time*(mu/eta) sum_j r_hat_i,j*(|r_i,j|-s0)
 *       = alpha sum_j r_hat_i,j*(|r_i,j|-s0)
 *
 * where alpha = stem_cycle_time*(mu/eta) = stem_cycle_time*meineke_lambda.
 *       s0    = natural length of the spring.

 * Length is scaled by natural length.
 * Time is scaled by a stem cell cycle time.
 *
 * meineke_lambda = mu (spring constant) / eta (damping) = 0.01 (from Meineke - note
 * that the value we use for Meineke lambda is completely different because we have
 * nondimensionalised)
 * 
 * The mesh should be surrounded by at least one layer of ghost nodes.  These are nodes which 
 * do not correspond to a cell, but are necessary for remeshing (because the remesher tries to 
 * create a convex hull of the set of nodes) and visualising purposes.  The mesh is passed into
 * the constructor and the class is told about the ghost nodes by using the method SetGhostNodes. 
 */

class CryptSimulation2D
{
private:
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<2,2> &mrMesh;
    
    bool mIncludeRandomBirth;
    bool mIncludeVariableRestLength;
    /**< A boolean saying to fix all four boundaries.*/  
    bool mFixedBoundaries;    
    /**< A boolean saying whether to remesh at each timestep or not (defaults to true).*/    
    bool mReMesh;
     /**< A vector of bools saying whether a node is ghosted-ified or not.*/  
    std::vector <bool> mIsGhostNode; 
    unsigned mMaxCells;
    unsigned mMaxElements ;
    
    std::string mOutputDirectory;
    
    std::vector<MeinekeCryptCell> mCells;
    
    CancerParameters *mpParams;
    
public:

    /** Constructor
     *  @param cells is defaulted to the empty vector, in which case SetIncludeRandomBirth()
     *  should be called for any birth to happen.
     */
    CryptSimulation2D(ConformingTetrahedralMesh<2,2> &rMesh,
                      std::vector<MeinekeCryptCell> cells = std::vector<MeinekeCryptCell>())
            : mrMesh(rMesh),
            mCells(cells)
    {
        mpParams = CancerParameters::Instance();
        mDt = 1.0/(mpParams->GetStemCellCycleTime()*120); // NOTE: hardcoded 120?
        mEndTime = 5.0;
        
        //srandom(time(NULL));
        srandom(0);
        mpParams->SetMeinekeLambda(15.0);
        
        mIncludeRandomBirth = false;
        mIncludeVariableRestLength = false;
        mFixedBoundaries =false;
        mOutputDirectory = "";
        // Set up the ghost nodes bool.  Assume initially that the maximum number of nodes is
        // ten times the mesh size.  Note that more memory is allocated later, if necessary.
        mIsGhostNode.resize(10*mrMesh.GetNumAllNodes()); // Note the hard-coding of 10.
        for (unsigned i=0; i<mIsGhostNode.size(); i++)
        {
            mIsGhostNode[i] = false;
        }
        
        mReMesh = true;
        
        mMaxCells = 10*mrMesh.GetNumNodes();
        mMaxElements = 10*mrMesh.GetNumElements();
    }
    
    void SetDt(double dt)
    {
        assert(dt>0);
        mDt=dt;
    }
    
    void SetEndTime(double endTime)
    {
        assert(endTime>0);
        mEndTime=endTime;
    }
    
//    /**
//     *  Call this before Solve() if no cells have been specified. Randomly adds a new
//     *  node every 1 time unit, starting 0.1
//     */
//    void SetIncludeRandomBirth()
//    {
//        mIncludeRandomBirth = true;
//    }

    void SetOutputDirectory(std::string outputDirectory)
    {
        mOutputDirectory = outputDirectory;
    }
    
    /**
     *  Call this before Solve() to simulate cell growth after cell division.
     *  (will eventually become SetIncludeCellBirth() and then become the default)
     */
    void SetIncludeVariableRestLength()
    {
        mIncludeVariableRestLength = true;
    }
    
    void SetMaxCells(unsigned maxCells)
    {
        mMaxCells = maxCells;
    }
    
    void SetMaxElements(unsigned maxElements)
    {
        mMaxElements = maxElements;
    }
    
    /**
    *  Call this before Solve() to fix the boundary of the mesh
    */
    void SetFixedBoundaries()
    {
        mFixedBoundaries = true;
    }
    
    /**
     * Main Solve method.
     * 
     * Once CryptSimulation object has been set up, call this to run simulation
     */
    void Solve()
    {
        if (mOutputDirectory=="")
        {
            EXCEPTION("OutputDirectory not set");
        }
        
        
         //Creating Column Data Writer Handler
        ColumnDataWriter tabulated_node_writer(mOutputDirectory+"Results", "tabulated_node_results");
        ColumnDataWriter tabulated_element_writer(mOutputDirectory+"Results", "tabulated_element_results");
        
        int time_var_id = tabulated_node_writer.DefineUnlimitedDimension("Time","hours");
        int time_var_id_elem = tabulated_element_writer.DefineUnlimitedDimension("Time","hours");
        
        
        std::vector<int> type_var_ids;
        std::vector<int> x_position_var_ids, y_position_var_ids;
        
        
        type_var_ids.resize(mMaxCells);
        x_position_var_ids.resize(mMaxCells);
        y_position_var_ids.resize(mMaxCells);
        
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // set up columns
        for (unsigned cell=0; cell<mMaxCells; cell++)
        {
            std::stringstream cell_type_var_name;
            std::stringstream cell_x_position_var_name;
            std::stringstream cell_y_position_var_name;
            cell_type_var_name << "cell_type_" << cell;
            cell_x_position_var_name << "cell_x_position_" << cell;
            cell_y_position_var_name << "cell_y_position_" << cell;            
            x_position_var_ids[cell]=tabulated_node_writer.DefineVariable(cell_x_position_var_name.str(),"rest_spring_length");
            y_position_var_ids[cell]=tabulated_node_writer.DefineVariable(cell_y_position_var_name.str(),"rest_spring_length");
            type_var_ids[cell]=tabulated_node_writer.DefineVariable(cell_type_var_name.str(),"dimensionless");
        }
        
        
        tabulated_node_writer.EndDefineMode();
        
        
        // Set up columns for element writer
        
        //std::vector<int> type_var_ids;
        std::vector<int> nodeA_var_ids, nodeB_var_ids, nodeC_var_ids;
        
        
        //type_var_ids.resize(mMaxCells);
        nodeA_var_ids.resize(mMaxElements);
        nodeB_var_ids.resize(mMaxElements);
        nodeC_var_ids.resize(mMaxElements);
        
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // set up columns
        for (unsigned elem_index = 0; elem_index<mMaxElements; elem_index++)
        {
            std::stringstream nodeA_var_name;
            std::stringstream nodeB_var_name;
            std::stringstream nodeC_var_name;
            
            
            nodeA_var_name << "nodeA_" << elem_index;
            nodeB_var_name << "nodeB_" << elem_index;
            nodeC_var_name << "nodeC_" << elem_index;
            
            nodeA_var_ids[elem_index] = tabulated_element_writer.DefineVariable(nodeA_var_name.str(),"dimensionless");
            nodeB_var_ids[elem_index] = tabulated_element_writer.DefineVariable(nodeB_var_name.str(),"dimensionless");
            nodeC_var_ids[elem_index] = tabulated_element_writer.DefineVariable(nodeC_var_name.str(),"dimensionless");
        }
                
        tabulated_element_writer.EndDefineMode();
        
        
        
        NodeMap map(mrMesh.GetNumAllNodes());
        double time = 0.0;
        double time_since_last_birth = 0.9;
        
        int num_births = 0;
        int num_deaths = 0;
        
        
        
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        

        std::vector<double> new_point_position(mrMesh.GetNumAllNodes());
        
        
        std::vector<bool> sloughed_cells(10*mrMesh.GetNumAllNodes());
        for (unsigned i=0; i<sloughed_cells.size(); i++)
        {
            sloughed_cells[i] = false;
        }
        
        static int step_number=0;
        
        // Creating Simple File Handler
        OutputFileHandler output_file_handler(mOutputDirectory);
        out_stream p_node_file = output_file_handler.OpenOutputFile("nodes.dat");
        
        out_stream p_element_file = output_file_handler.OpenOutputFile("elements.dat");



        //std::cout<< "\n mOutputDirectory " << mOutputDirectory ;
        
        while (time < mEndTime)
        {
            //std::cout << "** TIME = " << time << " **\n";
            
            // Cell birth
            if (!mCells.empty())
            {
                for (unsigned i=0; i<mCells.size(); i++)
                {
                    if(mrMesh.GetNodeAt(i)->IsDeleted()) continue; // Skip deleted cells
                    
                    if(mIsGhostNode[i]) continue;
                    // Check for this cell dividing
                    if(mCells[i].ReadyToDivide(time*mpParams->GetStemCellCycleTime()))
                    {
                        // Create new cell
                        MeinekeCryptCell new_cell = mCells[i].Divide(time*mpParams->GetStemCellCycleTime());
                        // Add new node to mesh
                        Node<2> *p_our_node = mrMesh.GetNodeAt(i);
                        
                        //Note: May need to check which side element is put esp. at the ends
                        Element<2,2> *p_element = mrMesh.GetElement(p_our_node->GetNextContainingElementIndex());
                        
                        //unsigned new_node_index = AddNodeToElement(p_element);
                        double new_x_value = (1.0/3.0)*(p_element->GetNode(0)->GetPoint().rGetLocation()[0]
                                                        +  p_element->GetNode(1)->GetPoint().rGetLocation()[0]
                                                        +  p_element->GetNode(2)->GetPoint().rGetLocation()[0] );
                                                        
                        double new_y_value = (1.0/3.0)*(p_element->GetNode(0)->GetPoint().rGetLocation()[1]
                                                        +  p_element->GetNode(1)->GetPoint().rGetLocation()[1]
                                                        +  p_element->GetNode(2)->GetPoint().rGetLocation()[1] );
                                                        
                        Point<2> new_point(new_x_value, new_y_value);
                        //std::cout<< "try1" <<std::endl<< std::flush;
                        unsigned new_node_index = mrMesh.RefineElement(p_element, new_point);
                        //std::cout<< "try2" <<std::endl<< std::flush;
                        // Update cells vector
                        new_cell.SetNodeIndex(new_node_index);
                        if (new_node_index == mCells.size())
                        {
                            mCells.push_back(new_cell);
                        }
                   
                        else
                        {
                            mCells[new_node_index] = new_cell;
                        }
                        
                        // Update size of IsGhostNode if necessary
                        if((int)mrMesh.GetNumNodes() > (int)mIsGhostNode.size())
                        {
                            mIsGhostNode.resize(mrMesh.GetNumNodes());
                            mIsGhostNode[new_node_index] = false;
                        }
                        num_births++;
                        //std::cout<< "num_births=" << num_births <<std::endl<< std::flush;
                    }
                }
            }
            
            // calculate node velocities
            std::vector<std::vector<double> > drdt(mrMesh.GetNumAllNodes());
            for (int i=0; i<mrMesh.GetNumAllNodes(); i++)
            {
                drdt[i].resize(2);
            }
                        
            double rest_length=1.0;
            
            ////////////////////////////////////////////////////////////////////
            // loop over element and for each one loop over it's three edges
            ////////////////////////////////////////////////////////////////////
            for (int elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
            {
                Element<2,2>* p_element = mrMesh.GetElement(elem_index);
                if(!p_element->IsDeleted())
                {
                    for (int k=0; k<3; k++)
                    {
                        int nodeA, nodeB;
                        if(k<2)
                        {
                            nodeA=k;
                            nodeB=k+1;
                        }
                        else
                        {
                            nodeA=2;
                            nodeB=0;
                        }
                        
                        c_vector<double, 2> drdt_contribution;
                        
                        c_vector<double, 2> unit_difference;
                        unit_difference(0)=p_element->GetNodeLocation(nodeB,0)-p_element->GetNodeLocation(nodeA,0);
                        unit_difference(1)=p_element->GetNodeLocation(nodeB,1)-p_element->GetNodeLocation(nodeA,1);
                        double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1));
                        
                        unit_difference = unit_difference/distance_between_nodes;
                        
                        drdt_contribution = mpParams->GetAlpha() *  unit_difference  * (distance_between_nodes - rest_length) ;
                        
                        assert(!p_element->GetNode(nodeA)->IsDeleted());
                        assert(!p_element->GetNode(nodeB)->IsDeleted());

                        // Assume that if both nodes are real, or both are ghosts, then they both
                        // exert forces on each other, but if one is real and one is ghost then
                        // the real node exerts a force on the ghost node, but the ghost node 
                        // does NOT exert a force on the real node.   
                        if(mIsGhostNode[p_element->GetNodeGlobalIndex(nodeA)] == false)
                        {
                            drdt[ p_element->GetNode(nodeB)->GetIndex()][0] -= drdt_contribution(0);
                            drdt[ p_element->GetNode(nodeB)->GetIndex()][1] -= drdt_contribution(1);

                            if(mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)] == false)
                            {
                                drdt[ p_element->GetNode(nodeA)->GetIndex()][0] += drdt_contribution(0);
                                drdt[ p_element->GetNode(nodeA)->GetIndex()][1] += drdt_contribution(1);
                            }
                        }
                        else
                        {
                            drdt[ p_element->GetNode(nodeA)->GetIndex()][0] += drdt_contribution(0);
                            drdt[ p_element->GetNode(nodeA)->GetIndex()][1] += drdt_contribution(1);
 
                            if(mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)] == true)
                            {
                               drdt[ p_element->GetNode(nodeB)->GetIndex()][0] -= drdt_contribution(0);
                               drdt[ p_element->GetNode(nodeB)->GetIndex()][1] -= drdt_contribution(1);
                            }
                        }       
                    }
                }
            }
            
            ////////////////////////////////////////////////////////////////////////////////////////
            // Also loop over boundary edges so that all edges have been looped over exactly twice.
            ////////////////////////////////////////////////////////////////////////////////////////
            ConformingTetrahedralMesh<2,2>::BoundaryElementIterator elem_iter
               = mrMesh.GetBoundaryElementIteratorBegin();
            
            while ( elem_iter != mrMesh.GetBoundaryElementIteratorEnd() )
            {
                BoundaryElement<1,2>* p_edge = *elem_iter;
                if(!p_edge->IsDeleted())
                {
                    c_vector<double, 2> drdt_contribution;
                    
                    c_vector<double, 2> unit_difference;
                    
                    int nodeA = 0;
                    int nodeB = 1;
                    
                    unit_difference(0)=p_edge->GetNodeLocation(nodeB,0)-p_edge->GetNodeLocation(nodeA,0);
                    unit_difference(1)=p_edge->GetNodeLocation(nodeB,1)-p_edge->GetNodeLocation(nodeA,1);
                    double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1));
                    
                    unit_difference=unit_difference/distance_between_nodes;
                    
                    drdt_contribution = mpParams->GetAlpha() *  unit_difference  * (distance_between_nodes - rest_length);
                       
                    assert(!p_edge->GetNode(nodeA)->IsDeleted());
                    assert(!p_edge->GetNode(nodeB)->IsDeleted());
                        
                    // Assume that if both nodes are real, or both are ghosts, then they both
                    // exert forces on each other, but if one is real and one is ghost then
                    // the real node exerts a force on the ghost node, but the ghost node 
                    // does NOT exert a force on the real node.   
                    if(mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeA)] == false)
                    {
                        drdt[ p_edge->GetNode(nodeB)->GetIndex()][0] -= drdt_contribution(0);
                        drdt[ p_edge->GetNode(nodeB)->GetIndex()][1] -= drdt_contribution(1);

                        if(mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)] == false)
                        {
                            drdt[ p_edge->GetNode(nodeA)->GetIndex()][0] += drdt_contribution(0);
                            drdt[ p_edge->GetNode(nodeA)->GetIndex()][1] += drdt_contribution(1);
                        }
                    }
                    else
                    {
                        drdt[ p_edge->GetNode(nodeA)->GetIndex()][0] += drdt_contribution(0);
                        drdt[ p_edge->GetNode(nodeA)->GetIndex()][1] += drdt_contribution(1);
 
                        if(mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)] == true)
                        {
                           drdt[ p_edge->GetNode(nodeB)->GetIndex()][0] -= drdt_contribution(0);
                           drdt[ p_edge->GetNode(nodeB)->GetIndex()][1] -= drdt_contribution(1);
                        }
                    }                            
                }
                elem_iter++;
            }            
            
            ////////////////////////////////////////////////////////////////////////////////////
            // update node positions
            ////////////////////////////////////////////////////////////////////////////////////
            for (int index = 0; index<mrMesh.GetNumAllNodes(); index++)
            {                
                if(mFixedBoundaries)
                {
                    // All Boundaries x=0, x=crypt_width, y=0, y=crypt_length.
                    if(mrMesh.GetNodeAt(index)->rGetPoint()[1]>0)
                    {
                        if(mrMesh.GetNodeAt(index)->rGetPoint()[1]<mpParams->GetCryptLength())
                        {
                            if(mrMesh.GetNodeAt(index)->rGetPoint()[0]>0)
                            {
                                if(mrMesh.GetNodeAt(index)->rGetPoint()[0]<mpParams->GetCryptWidth())
                                {
                                    if(!mrMesh.GetNodeAt(index)->IsDeleted())
                                    {
                                        //  std::cerr<<"Updating index "<<index<<"\n";
                                        Point<2> old_point = mrMesh.GetNodeAt(index)->rGetPoint();
                                        Point<2> new_point;
                                        
                                        // note factor of 0.5 in the update because drdt was twice
                                        // as large as it should be since edges were looped over twice.
                                        new_point.rGetLocation()[0] = old_point[0] + 0.5*mDt*drdt[index][0]; // new_point_position[index];
                                        new_point.rGetLocation()[1] = old_point[1] + 0.5*mDt*drdt[index][1]; // new_point_position[index];
                                        mrMesh.SetNode(index, new_point, false);
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    // assume stem cells are fixed, or if there are no cells, fix node 0
                    if(mrMesh.GetNodeAt(index)->rGetPoint()[1]>0)
                    {
                        if(!mrMesh.GetNodeAt(index)->IsDeleted())
                        {
                            //  std::cerr<<"Updating index "<<index<<"\n";
                            Point<2> old_point = mrMesh.GetNodeAt(index)->rGetPoint();
                            Point<2> new_point;
                            new_point.rGetLocation()[0] = old_point[0] + 0.5*mDt*drdt[index][0]; // new_point_position[index];
                            new_point.rGetLocation()[1] = old_point[1] + 0.5*mDt*drdt[index][1]; // new_point_position[index];
                            mrMesh.SetNode(index, new_point, false);

                        }
                    }
                }
            }
            
            
            /////////////////////////////////////////////////////////////////////////
            // SLOUGHING BY DELETING NODES. When ReMesh() works properly this will
            // probably need to be brought back.
            /////////////////////////////////////////////////////////////////////////
            
            /*
            // Remove nodes that are beyond the crypt
            while (true)
            {
                bool sloughed_node = false;
                ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator it = mrMesh.GetBoundaryNodeIteratorEnd();
                    
                while (it != mrMesh.GetBoundaryNodeIteratorBegin())
                {
                    it--; 
                    //Node<2> *p_node=mrMesh.GetNodeAt(0);
                    Node<2> *p_node = *it;
                    if(!p_node->IsDeleted())
                    {
                     double x = p_node->rGetPoint()[0];
                     double y = p_node->rGetPoint()[1];
                     //unsigned sloughing_node_index=p_node->GetIndex();
                     if((x > 1)||(x<-1)||(y>1)) /// changed
                     {
                         unsigned boundary_element_index=p_node->GetNextBoundaryElementIndex();
                     	BoundaryElement<1,2>* p_boundary_element=mrMesh.GetBoundaryElement(boundary_element_index);
                     	unsigned target_node_index=p_boundary_element->GetNodeGlobalIndex(0);
                     	if(target_node_index==sloughing_node_index)
                     	{
                     		target_node_index=p_boundary_element->GetNodeGlobalIndex(1);
                        	}
                     	
                     	std::cout<<time << "\t"<<sloughing_node_index<<"\t"<<target_node_index<<"\n";
                         // It's fallen off
                         assert(!mrMesh.GetNodeAt(target_node_index)->IsDeleted());
                         mrMesh.SetNode(sloughing_node_index,target_node_index);
            
            mIsGhostNode[p_node->GetIndex()] = true;
            
                         num_deaths++;
                         //std::cout<< "num_deaths=" << num_deaths <<std::endl<< std::flush;
                         sloughed_node = true;
                        
                         break;
                     }
                    }
                }
                if(!sloughed_node) break;
            }
            */
            step_number++;
            
            ///////////////////////////////////////////////////////////////////////////////////
            // Alternate method of sloughing.  Turns boundary nodes into ghost nodes.
            ///////////////////////////////////////////////////////////////////////////////////
            for (int i=0; i<mrMesh.GetNumNodes(); i++)
            {
                Node<2> *p_node = mrMesh.GetNodeAt(i);
                if(!p_node->IsDeleted())
                {
                    double x = p_node->rGetPoint()[0];
                    double y = p_node->rGetPoint()[1];
                    
                    double crypt_length=mpParams->GetCryptLength();
                    double crypt_width=mpParams->GetCryptWidth();
                    
                    if( (x>crypt_width) || (x<0.0) || (y>crypt_length))
                    {
                        mIsGhostNode[p_node->GetIndex()] = true;
                        num_deaths++;
                        //std::cout<< "num_deaths=" << num_deaths <<std::endl<< std::flush;
                    }
                }
            }
            
            if( mReMesh )
            {
                mrMesh.ReMesh(map);
            }
            
            ////////////////////////////////////////////////////////////////////////////////
            // Write results to file
            ////////////////////////////////////////////////////////////////////////////////
            (*p_node_file) << time << "\t";
            (*p_element_file) << time << "\t";
            
            for (int index = 0; index<mrMesh.GetNumAllElements(); index++)
            {
                if (!mrMesh.GetElement(index)->IsDeleted())
                {
                    (*p_element_file) << mrMesh.GetElement(index)->GetNodeGlobalIndex(0)<< " " << mrMesh.GetElement(index)->GetNodeGlobalIndex(1)<< " "<< mrMesh.GetElement(index)->GetNodeGlobalIndex(2)<< " ";
                }
            }
            (*p_element_file) << "\n";
            
            
            
            
            
            int cell=0; // NB this is not the index in mCells, but the index in the mesh!
//            for (int index = 0; index<mrMesh.GetNumAllNodes(); index++)
//            {
//                if (!mrMesh.GetNodeAt(index)->IsDeleted())
//                {
//                    Point<2> point = mrMesh.GetNodeAt(index)->rGetPoint();
//                    (*p_node_file) << point.rGetLocation()[0] << " ";
//                    (*p_node_file) << point.rGetLocation()[1] << " ";
//                    
//                    cell++;
//                }
//            }
//            
//            (*p_node_file) << "\n";
//            
            
            // write results to file
            //(*p_node_file) << time << "\t";
            
            
            
            for (int index = 0; index<mrMesh.GetNumAllNodes(); index++)
            {
                int colour = 0; // all green if no cells have been passed in
                
                if(mIsGhostNode[index]==true)
                {
                    colour = 3; // visualizer treats '4' these as invisible
                }
                else if(mCells.size()>0)
                {
                    if(index < (int) mCells.size())
                    {
                        CryptCellType type = mCells[index].GetCellType();
                        
                        if(type == STEM)
                        {
                            colour = 0;
                        }
                        else if(type == TRANSIT)
                        {
                            colour = 1;
                        }
                        else
                        {
                            colour = 2;
                        }
                    }
                    else
                    {
                        colour = 2; //Fix for segmentation fault
                    }
                    
                }
                
                if(!mrMesh.GetNodeAt(index)->IsDeleted())
                {
                    Point<2> point = mrMesh.GetNodeAt(index)->rGetPoint();
                    (*p_node_file) << point.rGetLocation()[0] << " "<< point.rGetLocation()[1] << " " << colour << " ";
                }
            }
            (*p_node_file) << "\n";
            
            if (time/0.05 - floor(time/0.05) < 1e-6)
            {
                tabulated_node_writer.PutVariable(time_var_id, time);
                tabulated_element_writer.PutVariable(time_var_id_elem, time);
                
                cell=0; // NB this is not the index in mCells, but the index in the mesh!
                for (int index = 0; index<mrMesh.GetNumAllNodes(); index++)
                {
                    if (!mrMesh.GetNodeAt(index)->IsDeleted())
                    {
                        if (mCells.size() > 0)
                        {
                            CryptCellType type  = mCells[index].GetCellType();
                            if (type == STEM)
                            {
                                tabulated_node_writer.PutVariable(type_var_ids[cell], 0);
                            }
                            else if (type == TRANSIT)
                            {
                                tabulated_node_writer.PutVariable(type_var_ids[cell], 1);
                            }
                            else if (type == DIFFERENTIATED)
                            {
                                tabulated_node_writer.PutVariable(type_var_ids[cell], 2);
                            }
                            else
                            {
                                // should be impossible to get here, until cancer cells
                                // are implemented
    #define COVERAGE_IGNORE
                                assert(0);
    #undef COVERAGE_IGNORE
                            }
                        }
                        else
                        {
                            tabulated_node_writer.PutVariable(type_var_ids[cell], -1);
                        }
                        
                        Point<2> point = mrMesh.GetNodeAt(index)->rGetPoint();
                        tabulated_node_writer.PutVariable(x_position_var_ids[cell], point.rGetLocation()[0]);
                        tabulated_node_writer.PutVariable(y_position_var_ids[cell], point.rGetLocation()[1]);
                        cell++;
                    }
                }
                tabulated_node_writer.AdvanceAlongUnlimitedDimension();
                
                for (int elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
                {
                    if (!mrMesh.GetElement(elem_index)->IsDeleted())
                    {
                        tabulated_element_writer.PutVariable(nodeA_var_ids[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0));
                        tabulated_element_writer.PutVariable(nodeB_var_ids[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1));
                        tabulated_element_writer.PutVariable(nodeC_var_ids[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2));
                    }
                }
                tabulated_element_writer.AdvanceAlongUnlimitedDimension();
                
            }
            
           
            time += mDt;
            time_since_last_birth += mDt;
            
        }
        tabulated_node_writer.Close();
        tabulated_element_writer.Close();
    }
    
    
    /**
     * The mesh should be surrounded by at least one layer of ghost nodes.  These are nodes which 
     * do not correspond to a cell, but are necessary for remeshing (because the remesher tries to 
     * create a convex hull of the set of nodes) and visualising purposes.  The mesh is passed into
     * the constructor and the class is told about the ghost nodes by using this method. 
     */     
    void SetGhostNodes(std::vector<int> ghostNodeIndices)
    {
        for (unsigned i = 0; i<ghostNodeIndices.size(); i++)
        {
            mIsGhostNode[ghostNodeIndices[i]]=true;
        }
    } 
    
    /**
     * Get the mesh to be remeshed at every time step.
     */   
    void SetReMeshRule(bool remesh)
    {
        mReMesh = remesh;
    }
    
};

#endif /*CRYPTSIMULATION2D_HPP_*/
