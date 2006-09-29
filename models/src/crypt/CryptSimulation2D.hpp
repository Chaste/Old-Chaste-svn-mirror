#ifndef CRYPTSIMULATION2D_HPP_
#define CRYPTSIMULATION2D_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"
#include "CancerParameters.hpp"
#include <cmath>
#include <ctime>
#include <iostream>
#include "TrianglesMeshWriter.cpp"

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
 */

class CryptSimulation2D
{
private:
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<2,2> &mrMesh;
    
    bool mIncludeRandomBirth;
    bool mIncludeVariableRestLength;
    bool mFixedBoundaries;
    
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
        NodeMap map(mrMesh.GetNumAllNodes());
 		double time = 0.0;
        double time_since_last_birth = 0.9;
        
        int num_births = 0;
        int num_deaths = 0;
        
        std::vector<double> new_point_position(mrMesh.GetNumAllNodes());
        
        OutputFileHandler output_file_handler(mOutputDirectory);
        out_stream p_node_file = output_file_handler.OpenOutputFile("nodes.dat");
        out_stream p_element_file = output_file_handler.OpenOutputFile("elements.dat");
        
        std::vector<bool> sloughed_cells(10*mrMesh.GetNumAllNodes());
        for(unsigned i=0; i<sloughed_cells.size(); i++)
        {
        	sloughed_cells[i] = false;
        }
			
		static int step_number=0;
        
        while (time < mEndTime)
        {
        	//std::cout << "** TIME = " << time << " **\n";
        	
            // Cell birth
            if (mIncludeRandomBirth && time_since_last_birth > 1)
            {
//                AddRandomNode();
//                time_since_last_birth = 0 ;
            }
            else if (!mCells.empty())
            {
                for (unsigned i=0; i<mCells.size(); i++)
                {
                    if (mrMesh.GetNodeAt(i)->IsDeleted()) continue; // Skip deleted cells
                    // Check for this cell dividing
                    if (mCells[i].ReadyToDivide(time*mpParams->GetStemCellCycleTime()))
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
                        
                        unsigned new_node_index = mrMesh.RefineElement(p_element, new_point);
                        
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
                        num_births++;
                        //std::cout<< "num_births=" << num_births <<std::endl<< std::flush;
                    }
                }
            }
            
            // calculate node velocities
            std::vector<std::vector<double> > drdt(mrMesh.GetNumAllNodes());
            for(int i=0; i<mrMesh.GetNumAllNodes(); i++)
            {
	            drdt[i].resize(2);
            }
            
            
            // the following commented section is old 1d code:
            
//            if (mIncludeVariableRestLength && !mCells.empty())
//            {
//                for (int elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
//                {
//                    Element<1,1>* element = mrMesh.GetElement(elem_index);
//                    if (!element->IsDeleted())
//                    {
//                        c_vector<double, 2> drdt_contributions;
//                        
//                        double distance_between_nodes = fabs(element->GetNodeLocation(1,0) - element->GetNodeLocation(0,0));
//                        double unit_vector_backward = -1;
//                        double unit_vector_forward = 1;
//                        
//                        double age1 = mCells[element->GetNode(0)->GetIndex()].GetAge(time*mpParams->GetStemCellCycleTime());
//                        double age2 = mCells[element->GetNode(1)->GetIndex()].GetAge(time*mpParams->GetStemCellCycleTime());
//                        double rest_length=mpParams->GetNaturalSpringLength();
//                        double time_scale = mpParams->GetStemCellCycleTime();
//                        if (age1<1.0/time_scale && age2<1.0/time_scale && fabs(age1-age2)<1e-6)
//                        {
//                            /* Spring Rest Length Increases to 1 from 0.1 over 1 hour
//                             * This doesnt happen at present as when the full line is included the tests fail
//                             */
//                            rest_length=0.9*rest_length;//+0.1*age1*time_scale;
//                            assert(rest_length<=mpParams->GetNaturalSpringLength());
//                        }
//                        
//                        drdt_contributions(0) = mpParams->GetAlpha() *(  unit_vector_forward  * (distance_between_nodes - rest_length) );
//                        drdt_contributions(1) = mpParams->GetAlpha() *(  unit_vector_backward * (distance_between_nodes - rest_length) );
//                        
//                        drdt[ element->GetNode(0)->GetIndex() ] += drdt_contributions(0);
//						drdt[ element->GetNode(1)->GetIndex() ] += drdt_contributions(1);
//                    }
//                }
//            }
//            else
//            {
            
            // loop over element and for each one loop over it's three edges
			for (int elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
            {
            	Element<2,2>* p_element = mrMesh.GetElement(elem_index);
                if (!p_element->IsDeleted())
                {
                	for (int k=0; k<3; k++)
                	{
                		int nodeA, nodeB;
                		if (k<2)
                		{
                			nodeA=k;
                			nodeB=k+1; 
                		}
                		else
                		{
                			nodeA=2;
                			nodeB=0; 
                		}   
                		
                	    c_matrix<double, 2, 2> drdt_contributions;
               	      
               	        c_vector<double, 2> unit_difference;
               	        unit_difference(0)=p_element->GetNodeLocation(nodeB,0)-p_element->GetNodeLocation(nodeA,0);
               	        unit_difference(1)=p_element->GetNodeLocation(nodeB,1)-p_element->GetNodeLocation(nodeA,1);
               	        double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1));
               	      
               	        unit_difference(0)=unit_difference(0)/distance_between_nodes;
               	        unit_difference(1)=unit_difference(1)/distance_between_nodes;

						// if neither node is sloughed include force contribution
		               	if( (sloughed_cells[p_element->GetNodeGlobalIndex(nodeA)] == false) && (sloughed_cells[p_element->GetNodeGlobalIndex(nodeB)] == false))
	                	{
       		       	        // Refactor this to a vector operation
       		       	        double rest_length=mpParams->GetNaturalSpringLength();
            	   	        drdt_contributions(0,0) = mpParams->GetAlpha() *(  unit_difference(0)  * (distance_between_nodes - rest_length) );
               		        drdt_contributions(0,1) = mpParams->GetAlpha() *(  unit_difference(1)  * (distance_between_nodes - rest_length) );
               	    	    drdt_contributions(1,0) = mpParams->GetAlpha() *(  -unit_difference(0)  * (distance_between_nodes - rest_length) );
               	        	drdt_contributions(1,1) = mpParams->GetAlpha() *(  -unit_difference(1)  * (distance_between_nodes - rest_length) );
               	        
	               	        assert(!p_element->GetNode(nodeA)->IsDeleted());
    	           	        assert(!p_element->GetNode(nodeB)->IsDeleted());
        
        	       	        drdt[ p_element->GetNode(nodeA)->GetIndex()][0] += drdt_contributions(0,0);
            	   	        drdt[ p_element->GetNode(nodeA)->GetIndex()][1] += drdt_contributions(0,1);
               		        drdt[ p_element->GetNode(nodeB)->GetIndex()][0] += drdt_contributions(1,0);
               	    	    drdt[ p_element->GetNode(nodeB)->GetIndex()][1] += drdt_contributions(1,1);   
	                	}
                	}
                }
            }

            // Also loop over boundary edges so that all edges have been looped over exactly 
            // twice.
            ConformingTetrahedralMesh<2,2>::BoundaryElementIterator elem_iter 
               = mrMesh.GetBoundaryElementIteratorBegin(); 
            
	        while( elem_iter != mrMesh.GetBoundaryElementIteratorEnd() )
    	    {
    	    	BoundaryElement<1,2>* p_edge = *elem_iter;
    	    	if (!p_edge->IsDeleted())
                {
    	    	    c_matrix<double, 2, 2> drdt_contributions;
               	      
               	    c_vector<double, 2> unit_difference;
               	
               	    int nodeA = 0; 
               	    int nodeB = 1;
               	
               	    unit_difference(0)=p_edge->GetNodeLocation(nodeB,0)-p_edge->GetNodeLocation(nodeA,0);
	               	unit_difference(1)=p_edge->GetNodeLocation(nodeB,1)-p_edge->GetNodeLocation(nodeA,1);
	               	double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1));
	               	     
	               	unit_difference(0)=unit_difference(0)/distance_between_nodes;
	               	unit_difference(1)=unit_difference(1)/distance_between_nodes;
	               	      
	    			// if neither node is sloughed include force contribution
	               	if( (sloughed_cells[p_edge->GetNodeGlobalIndex(nodeA)] == false) && (sloughed_cells[p_edge->GetNodeGlobalIndex(nodeB)] == false))
                	{
		               	// Refactor this to a vector operation
		               	double rest_length=mpParams->GetNaturalSpringLength();
		               	drdt_contributions(0,0) = mpParams->GetAlpha() *(  unit_difference(0)  * (distance_between_nodes - rest_length) );
		               	drdt_contributions(0,1) = mpParams->GetAlpha() *(  unit_difference(1)  * (distance_between_nodes - rest_length) );
	    	           	drdt_contributions(1,0) = mpParams->GetAlpha() *(  -unit_difference(0)  * (distance_between_nodes - rest_length) );
	        	       	drdt_contributions(1,1) = mpParams->GetAlpha() *(  -unit_difference(1)  * (distance_between_nodes - rest_length) );
	               	
	             	  	assert(!p_edge->GetNode(nodeA)->IsDeleted());
	               		assert(!p_edge->GetNode(nodeB)->IsDeleted());
		               	              
		               	drdt[ p_edge->GetNode(nodeA)->GetIndex()][0] += drdt_contributions(0,0);
	    	           	drdt[ p_edge->GetNode(nodeA)->GetIndex()][1] += drdt_contributions(0,1);
	        	       	drdt[ p_edge->GetNode(nodeB)->GetIndex()][0] += drdt_contributions(1,0);
	            	   	drdt[ p_edge->GetNode(nodeB)->GetIndex()][1] += drdt_contributions(1,1);   
                	}
                }
                elem_iter++;
    	    }
           
           

            // update node positions
            for (int index = 0; index<mrMesh.GetNumAllNodes(); index++)
            {
            	// note factor of 0.5 in the update because drdt was twice 
            	// as large as it should be since edges were looped over twice.
            	
            	if(mFixedBoundaries)
            	{
            		// All Bonndaries x=0, x=crypt_width, y=0, y=crypt_length.
	                if (mrMesh.GetNodeAt(index)->rGetPoint()[1]>0)
	                {
	                	if (mrMesh.GetNodeAt(index)->rGetPoint()[1]<mpParams->GetCryptLength())
	                	{
	                		if (mrMesh.GetNodeAt(index)->rGetPoint()[0]>0)
	                		{
	                			if (mrMesh.GetNodeAt(index)->rGetPoint()[0]<mpParams->GetCryptWidth())
	                            {
	                                if (!mrMesh.GetNodeAt(index)->IsDeleted())
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
	                }
            	}
            	else
            	{
            		// assume stem cells are fixed, or if there are no cells, fix node 0
	                if (mrMesh.GetNodeAt(index)->rGetPoint()[1]>0)
	                {
	                    if (!mrMesh.GetNodeAt(index)->IsDeleted())
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
	                    if ((x > 1)||(x<-1)||(y>1)) /// changed
	                    {
	                        unsigned boundary_element_index=p_node->GetNextBoundaryElementIndex();
	                    	BoundaryElement<1,2>* p_boundary_element=mrMesh.GetBoundaryElement(boundary_element_index);
	                    	unsigned target_node_index=p_boundary_element->GetNodeGlobalIndex(0);
	                    	if (target_node_index==sloughing_node_index)
	                    	{
	                    		target_node_index=p_boundary_element->GetNodeGlobalIndex(1);
                        	}
	                    	
	                    	std::cout<<time << "\t"<<sloughing_node_index<<"\t"<<target_node_index<<"\n";
	                        // It's fallen off
	                        assert(!mrMesh.GetNodeAt(target_node_index)->IsDeleted());
	                        mrMesh.SetNode(sloughing_node_index,target_node_index);

							sloughed_cells[p_node->GetIndex()] = true;

	                        num_deaths++;
	                        //std::cout<< "num_deaths=" << num_deaths <<std::endl<< std::flush;
	                        sloughed_node = true;
	                        
	                        break;
	                    }
                    }
                }
                if (!sloughed_node) break;
            }
            */
            step_number++;

            // sloughing by noting which nodes are on boundary and removing their contributions
            // to force
            for(int i=0; i<mrMesh.GetNumNodes(); i++)
            {
                 Node<2> *p_node = mrMesh.GetNodeAt(i);
                 if(!p_node->IsDeleted())
                 {
                    double x = p_node->rGetPoint()[0];
                    double y = p_node->rGetPoint()[1];
                    
                    double crypt_length=mpParams->GetCryptLength();
                    double crypt_width=mpParams->GetCryptWidth();
                    
                    if ( (x>crypt_width) || (x<0.0) || (y>crypt_length))
                    {
						sloughed_cells[p_node->GetIndex()] = true;
                        num_deaths++;
                        //std::cout<< "num_deaths=" << num_deaths <<std::endl<< std::flush;
                    }
                }
            }



           if (mFixedBoundaries)
           {
	            //Re-mesh the mesh
				mrMesh.ReMesh(map);
           }
            
//          std::stringstream time_ss;
//	        time_ss << step_number++;
//			TrianglesMeshWriter<2,2> mesh_writer("","tempmesh"+time_ss.str());
//			mesh_writer.WriteFilesUsingMesh(mrMesh);	    


           
            // write results to file
            (*p_node_file) << time << "\t";
            (*p_element_file) << time << "\t";
            
            for (int index = 0; index<mrMesh.GetNumAllNodes(); index++)
            {
                int colour = 0; // all green if no cells have been passed in
                if (sloughed_cells[index]==true)
                {
                    colour = 4; // visualizer treats these as invisible
                }
                else if(mCells.size()>0)
                {
                    if (index < (int) mCells.size())
                    {
                        CryptCellType type = mCells[index].GetCellType();
            	
                	if (type==STEM)
                	{
                		colour = 0;
                	}
                	else if (type == TRANSIT)
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
                        colour=2; //Fix for segmentation fault
                    }
                    
                }
            	
                if (!mrMesh.GetNodeAt(index)->IsDeleted())
                {
                    Point<2> point = mrMesh.GetNodeAt(index)->rGetPoint();
                    (*p_node_file) << point.rGetLocation()[0] << " "<< point.rGetLocation()[1] << " " << colour << " ";
                }
            }
            (*p_node_file) << "\n";
            
            
            for (int index = 0; index<mrMesh.GetNumAllElements(); index++)
            {
            	if (!mrMesh.GetElement(index)->IsDeleted())
                {
                	(*p_element_file) << mrMesh.GetElement(index)->GetNodeGlobalIndex(0)<< " " << mrMesh.GetElement(index)->GetNodeGlobalIndex(1)<< " "<< mrMesh.GetElement(index)->GetNodeGlobalIndex(2)<< " ";
                }
            }
            (*p_element_file) << "\n";
            
            time += mDt;
            time_since_last_birth += mDt;
            
        }
    }
};

#endif /*CRYPTSIMULATION2D_HPP_*/
