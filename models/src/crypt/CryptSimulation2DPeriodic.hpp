#ifndef CRYPTSIMULATION2DPERIODIC_HPP_
#define CRYPTSIMULATION2DPERIODIC_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataWriter.hpp"
#include <cmath>
#include <ctime>
#include <iostream>
#include "TrianglesMeshWriter.cpp"
#include "Exception.hpp"
#include "SimulationTime.hpp"
#include "StochasticCellCycleModel.hpp"
#include "ColumnDataWriter.hpp"
#include "MeinekeCryptCellTypes.hpp"
#include "WntCellCycleModel.hpp"
#include "WntGradient.hpp"

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

class CryptSimulation2DPeriodic
{
private:
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<2,2> &mrMesh;
    
    bool mIncludeRandomBirth;
    bool mIncludeVariableRestLength;

    /**< A boolean saying to fix all four boundaries.*/  
    bool mFixedBoundaries;    
    
    /**< A boolean saying to run the simulation with no birth (defaults to false). */
    bool mNoBirth;

    /**< A boolean saying whether to remesh at each timestep or not (defaults to true).*/    
    bool mReMesh;
    
    /**< A boolean saying whether Wnt signalling is included or not (defaults to false).*/    
    bool mWntIncluded;
    
    /**< A boolean saying whether the mesh is periodic or not (defaults to false).*/    
    bool mPeriodicSides;

     /**< A vector of bools saying whether a node is ghosted-ified or not.*/  
    std::vector <bool> mIsGhostNode; 
    std::vector <bool> mIsPeriodicNode;
    
    /**< A std::vector of unsigned integers giving the node indexes of nodes on the left boundary. */ 
    std::vector <unsigned> mLeftCryptBoundary;
    std::vector <unsigned> mOldLeftCryptBoundary;
    /**< A std::vector of unsigned integers giving the node indexes of nodes on the right boundary. */
    std::vector <unsigned> mRightCryptBoundary;
    std::vector <unsigned> mOldRightCryptBoundary;
    /**< A std::vector of unsigned integers giving the node indexes of nodes on the boundary. */
    std::vector <unsigned> mCryptBoundary;
    std::vector <unsigned> mOldCryptBoundary;
    
    unsigned mMaxCells;
    unsigned mMaxElements;
    
    std::string mOutputDirectory;
    
    std::vector<MeinekeCryptCell> mCells;

    RandomNumberGenerator *mpRandomNumberGenerator;    
    bool mCreatedRng;
    
    CancerParameters *mpParams;
    
    WntGradientType mWntGradient;
    
public:

    /** Constructor
     *  @param cells is defaulted to the empty vector, in which case SetIncludeRandomBirth()
     *  should be called for any birth to happen.
     */
    CryptSimulation2DPeriodic(ConformingTetrahedralMesh<2,2> &rMesh,
                      std::vector<MeinekeCryptCell> cells = std::vector<MeinekeCryptCell>(),
                      RandomNumberGenerator *pGen = NULL)
            : mrMesh(rMesh),
              mCells(cells)
    {
        if (pGen!=NULL)
        {
            mpRandomNumberGenerator = pGen;
            mCreatedRng = false;
        }
        else
        {
            mpRandomNumberGenerator = new RandomNumberGenerator;
            mCreatedRng = true;
        }
        
        mpParams = CancerParameters::Instance();
    
        mDt = 1.0/(120.0);
        mEndTime = 120.0; //hours
        
        
        //srandom(time(NULL));
        srandom(0);
        mpParams->SetMeinekeLambda(15.0);
        
        mIncludeRandomBirth = false;
        mIncludeVariableRestLength = false;
        mFixedBoundaries = false;
        mOutputDirectory = "";
        
        // Set up the ghost nodes bool.  Assume initially that the maximum number of nodes is
        // ten times the mesh size.  Note that more memory is allocated later, if necessary.
        mIsGhostNode.resize(10*mrMesh.GetNumAllNodes()); // Note the hard-coding of 10.
        mIsPeriodicNode.resize(mIsGhostNode.size());
        for (unsigned i=0; i<mIsGhostNode.size(); i++)
        {
            mIsGhostNode[i] = false;
            mIsPeriodicNode[i] = false;
        }

        // defaults
        mReMesh = true;
        mNoBirth = false;
        mMaxCells = 10*mrMesh.GetNumNodes();
        mMaxElements = 10*mrMesh.GetNumElements();
        mWntIncluded = false;
        mPeriodicSides = true;
        mWntGradient=NONE;
    }
    
    /**
     * Free any memory allocated by the constructor
     */
    ~CryptSimulation2DPeriodic()
    {
        if (mCreatedRng)
        {
            delete mpRandomNumberGenerator;
        }
        SimulationTime::Destroy();
    }
    
    
    void SetDt(double dt)
    {
        assert(dt>0);
        mDt=dt;
    }
    
    /** 
     * Sets the end time and resets the timestep to be endtime/100
     */
    void SetEndTime(double endTime)
    {
        assert(endTime>0);
        mEndTime=endTime;
    }
    
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
     *  Call this before Solve() to set the boundary conditions
     */
    void SetPeriodicSides(bool periodicSides)
    {
        mPeriodicSides = periodicSides;
    }
    
    /** 
     *  Get the cells vector
     */
    std::vector<MeinekeCryptCell> GetCells()
    {
        assert(mCells.size()>0);
        return mCells;
    }
    
    /** A standard vector of booleans saying whether
     *  a node is a ghost or not
     */
    std::vector <bool> GetGhostNodes()
    {
        return mIsGhostNode;
    }
    
    std::vector<unsigned> GetLeftCryptBoundary()
    {
        return mLeftCryptBoundary;
    }
    
    std::vector<unsigned> GetRightCryptBoundary()
    {
        return mRightCryptBoundary;
    }
    
    std::vector<unsigned> GetCryptBoundary()
    {
        return mCryptBoundary;
    }
    
    // This automatically sets this to be a wnt dependent simulation
    // you should supply cells with a wnt cell cycle...
    void SetWntGradient(WntGradientType wntGradient)
    {
    	mWntIncluded = true;
    	mWntGradient = wntGradient;
    }
    /**
     * Main Solve method.
     * 
     * Once CryptSimulation object has been set up, call this to run simulation
     */
    void Solve()
    {
    	WntGradient wnt_gradient(mWntGradient);
        if (mOutputDirectory=="")
        {
            EXCEPTION("OutputDirectory not set");
        }
        ///////////////////////////////////////////////////////////
        // Set up Column Data Writer 
        ///////////////////////////////////////////////////////////
        ColumnDataWriter tabulated_node_writer(mOutputDirectory+"Results", "tabulated_node_results");
        ColumnDataWriter tabulated_element_writer(mOutputDirectory+"Results", "tabulated_element_results");
        
        int time_var_id = tabulated_node_writer.DefineUnlimitedDimension("Time","hours");
        int time_var_id_elem = tabulated_element_writer.DefineUnlimitedDimension("Time","hours");
        
        std::vector<int> type_var_ids;
        std::vector<int> x_position_var_ids, y_position_var_ids;
        
        type_var_ids.resize(mMaxCells);
        x_position_var_ids.resize(mMaxCells);
        y_position_var_ids.resize(mMaxCells);

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
        
        
        

        int num_time_steps = (int)(mEndTime/mDt+0.5);
        
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);        
         
        int num_births = 0;
        int num_deaths = 0;
        std::cout << num_births << "\n" ;
        
        std::vector<double> new_point_position(mrMesh.GetNumAllNodes());
        
        
        static int step_number=0;
        
        // Creating Simple File Handler
        OutputFileHandler output_file_handler(mOutputDirectory);

        out_stream p_node_file = output_file_handler.OpenOutputFile("results.viznodes");
        out_stream p_element_file = output_file_handler.OpenOutputFile("results.vizelements");


        int counter = 0;
        
		int periodic_division_buffer = 0;
        
        if(mPeriodicSides)
        {
        	CalculateCryptBoundary();
        	if(mLeftCryptBoundary.size()<1 || mRightCryptBoundary.size()<1)
        	{
        			EXCEPTION("Periodic Simulation but mesh is not periodic\nIf you want a non-periodic simulation use SetPeriodicSides(false)");
        	}
        }
        /* Age the cells to the correct time (cells set up with negative birth dates 
         * to gives some that are almost ready to divide).
         * 
         * For some strange reason this seems to take about 3 minutes for a realistic Wnt-Crypt.
         * Not sure why - when the same thing was evaluated in a test it seemed almost instant.
         */
        bool temp;
        if (!mCells.empty())
	    {
	    	for(unsigned i=0; i<mCells.size() ; i++)
        	{
		    	if(mIsGhostNode[i]) continue;
		    	//std::cout << "Perparing Cell "<< i << std::endl;
		    	Node<2> *p_our_node = mrMesh.GetNode(i);
		        double y = p_our_node->GetPoint()[1];
		        std::vector<double> cell_cycle_influences;
		        if(mWntIncluded)
		        {
			    	double wnt_stimulus = wnt_gradient.GetWntLevel(y);
			    	cell_cycle_influences.push_back(wnt_stimulus);
		        }
                temp = mCells[i].ReadyToDivide(cell_cycle_influences);		
            }
	    }
	    
        /////////////////////////////////////////////////////////////////////
        // Main time loop
        /////////////////////////////////////////////////////////////////////
        while (p_simulation_time->GetTimeStepsElapsed() < num_time_steps)
        {
		    std::cout << "** TIME = " << p_simulation_time->GetDimensionalisedTime() << " **\n";
            
            ///////////////////////////////////////////////////////
            // Cell birth
            ///////////////////////////////////////////////////////
            periodic_division_buffer--;
            if (mNoBirth==false)
            {
	            if (!mCells.empty())
	            {
	            	
	                for (unsigned i=0; i<mCells.size(); i++)
	                {
	                    bool skip = false;
	                    if(mrMesh.GetNode(i)->IsDeleted()) skip=true; // Skip deleted cells
	                    if(mIsGhostNode[i]) skip=true; // Skip Ghost nodes
						bool periodic_cell=false;
						unsigned periodic_index = 0;
						if(mPeriodicSides)
						{
							for (unsigned j=0 ; j < mRightCryptBoundary.size() ; j++)
		                    {
		                    	if (mRightCryptBoundary[j]==i)
		                    	{	// Only allow one periodic boundary to have divisions...
		                    		skip=true;
		                    	}	
		                    	if (mLeftCryptBoundary[j]==i)
		                    	{
		                    		if(periodic_division_buffer>0)
		                    		{
                                        #define COVERAGE_IGNORE
                                        // Only allow one periodic cell division per timestep so that mesh can catch up with it.
		                    			skip=true;	
                                        #undef COVERAGE_IGNORE
		                    		}
		                    		periodic_cell = true;
									periodic_index = j;
		                    	}
		                    }
						}
	                    if(skip) continue;
						
	                    // Check for this cell dividing
                    	Node<2> *p_our_node = mrMesh.GetNode(i);
	                    double y = p_our_node->GetPoint()[1];
	                    std::vector<double> cell_cycle_influences;
	                    
	                    if(mWntIncluded)
	                    {
		                    double wnt_stimulus = wnt_gradient.GetWntLevel(y);
		                    cell_cycle_influences.push_back(wnt_stimulus);
	                    }
	                    
						//std::cout << "On cell "<< i << std::endl;
						if(mCells[i].ReadyToDivide(cell_cycle_influences))
	                    {
	                    	std::cout << "cell division at node " << i << "\n";
	                        // Create new cell
	                        MeinekeCryptCell new_cell = mCells[i].Divide();
							if(mPeriodicSides && periodic_cell)
	                        {	
                                #define COVERAGE_IGNORE
	                        	std::cout << "Periodic Division\n";
	                        	periodic_division_buffer=3;
	                        	//Make sure the image cell knows it has just divided and aged a generation
	                        	mCells[mRightCryptBoundary[periodic_index]]=mCells[mLeftCryptBoundary[periodic_index]];
                                #undef COVERAGE_IGNORE
	                        }
	                        
	
	                        // Add new node to mesh
	                        double x = p_our_node->GetPoint()[0];
	                        
	                        
	                        // For all cells find an element that isn't a ghost element to put a new cell into.
							unsigned element_number = mpRandomNumberGenerator->randMod( p_our_node->GetNumContainingElements() );
	
	                        Element<2,2>* p_element = mrMesh.GetElement(p_our_node->GetNextContainingElementIndex());
	                        for(unsigned j=0; j<element_number; j++)
	                        {
	                            p_element = mrMesh.GetElement(p_our_node->GetNextContainingElementIndex());
	                        }
	                        
	                        bool is_ghost_element = (    (mIsGhostNode[p_element->GetNodeGlobalIndex(0)]) 
	                                                  || (mIsGhostNode[p_element->GetNodeGlobalIndex(1)]) 
	                                                  ||  (mIsGhostNode[p_element->GetNodeGlobalIndex(2)]) );
							// Make sure only one of the nodes is periodic (the cell that's dividing)
							bool is_periodic_element = ( ((mIsPeriodicNode[p_element->GetNodeGlobalIndex(0)]) 
	                                                  && (mIsPeriodicNode[p_element->GetNodeGlobalIndex(1)])) 
	                                                  || ((mIsPeriodicNode[p_element->GetNodeGlobalIndex(0)]) 
	                                                  && (mIsPeriodicNode[p_element->GetNodeGlobalIndex(2)]))
	                                                  || ((mIsPeriodicNode[p_element->GetNodeGlobalIndex(1)]) 
	                                                  && (mIsPeriodicNode[p_element->GetNodeGlobalIndex(2)])));
							
							
			                unsigned counter = 0;
	                        while (is_ghost_element || is_periodic_element)
	                        {
                                #define COVERAGE_IGNORE 
	                            p_element = mrMesh.GetElement(p_our_node->GetNextContainingElementIndex());
	                            is_ghost_element = (    (mIsGhostNode[p_element->GetNodeGlobalIndex(0)]) 
	                                                 || (mIsGhostNode[p_element->GetNodeGlobalIndex(1)]) 
	                                                 || (mIsGhostNode[p_element->GetNodeGlobalIndex(2)]) );
	   
								is_periodic_element = ( ((mIsPeriodicNode[p_element->GetNodeGlobalIndex(0)]) 
	                                                  && (mIsPeriodicNode[p_element->GetNodeGlobalIndex(1)])) 
	                                                  || ((mIsPeriodicNode[p_element->GetNodeGlobalIndex(0)]) 
	                                                  && (mIsPeriodicNode[p_element->GetNodeGlobalIndex(2)]))
	                                                  || ((mIsPeriodicNode[p_element->GetNodeGlobalIndex(1)]) 
	                                                  && (mIsPeriodicNode[p_element->GetNodeGlobalIndex(2)])));
	            				counter++;
				                if(counter>p_our_node->GetNumContainingElements()+1)
	            				{
	            					if(periodic_cell)
	            					{// Swap to the image node - that might have a 
	            						//non-periodic element to put new cell into.
	            						if(i==mRightCryptBoundary[periodic_index])
	            						{
	            							assert(0);	
	            						}
	            						p_our_node = mrMesh.GetNode(mRightCryptBoundary[periodic_index]);
	                        			x = p_our_node->GetPoint()[0];
										y = p_our_node->GetPoint()[1];
	            					}
	            					else
	            					{
	            				    	// somehow every connecting element is a ghost element. quit to
	            				    	// avoid infinite loop
	            				    	assert(0);
	            					}
	            				}
                                #undef COVERAGE_IGNORE
	                        }
	                        
	                        std::cout << "New cell being intoduced into element with nodes \n";
	                        std::cout << p_element->GetNodeGlobalIndex(0) << "\t" << p_element->GetNodeGlobalIndex(1) << "\t" <<p_element->GetNodeGlobalIndex(2) << "\n";
	                        	
	                        double x_centroid = (1.0/3.0)*(p_element->GetNode(0)->GetPoint().rGetLocation()[0]
	                                                        +  p_element->GetNode(1)->GetPoint().rGetLocation()[0]
	                                                        +  p_element->GetNode(2)->GetPoint().rGetLocation()[0] );
	                                                        
	                        double y_centroid = (1.0/3.0)*(p_element->GetNode(0)->GetPoint().rGetLocation()[1]
	                                                        +  p_element->GetNode(1)->GetPoint().rGetLocation()[1]
	                                                        +  p_element->GetNode(2)->GetPoint().rGetLocation()[1] );
	                                                        
	                        
	                        // check the new point is in the triangle
	                        double distance_from_node_to_centroid =  sqrt(  (x_centroid - x)*(x_centroid - x)
	                                                                      + (y_centroid - y)*(y_centroid - y) );
	                        
	                        // we assume the new cell is a distance 0.1 away from the old.
	                        // however, to avoid crashing in usual situations we check this
	                        // new position is actually in the triangle being refined.
	                        double distance_of_new_cell_from_parent = 0.1;
	                        if(distance_from_node_to_centroid < (2.0/3.0)*0.1)
	                        {
                                #define COVERAGE_IGNORE
	                            distance_of_new_cell_from_parent = (3.0/2.0)*distance_from_node_to_centroid;
                                #undef COVERAGE_IGNORE
	                        }
	                        
	                        double new_x_value = x + distance_of_new_cell_from_parent*(x_centroid-x);
	                        double new_y_value = y + distance_of_new_cell_from_parent*(y_centroid-y);
	                        
	                        std::cout << "Parent node at x = " << x << "  y = " << y << "\n";
	                        std::cout << "Daughter node at x = " << new_x_value << "  y = " << new_y_value << "\n";
	
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
	                            #define COVERAGE_IGNORE
                                mCells[new_node_index] = new_cell;
                                #undef COVERAGE_IGNORE
	                        }
	                        //mCells[new_node_index].SetBirthTime();
	                        
	                        // Update size of IsGhostNode if necessary
	                        if((int)mrMesh.GetNumNodes() > (int)mIsGhostNode.size())
	                        {
                                #define COVERAGE_IGNORE
	                            mIsGhostNode.resize(mrMesh.GetNumNodes());
	                            mIsGhostNode[new_node_index] = false;
                                #undef COVERAGE_IGNORE
	                        }
	                        num_births++;
	                        //std::cout<< "num_births=" << num_births <<std::endl<< std::flush;
	                        if( mReMesh )
	            			{
	                			ReMesh();
	            			}
	                    
						}
	                }
	            }
        	}

			
			
            
            
            //////////////////////////////////////////////////////////////////////////
            //                    calculate node velocities
            //////////////////////////////////////////////////////////////////////////
            std::vector<std::vector<double> > drdt(mrMesh.GetNumAllNodes());
            for (int i=0; i<mrMesh.GetNumAllNodes(); i++)
            {
                drdt[i].resize(2);
            }
            
                        
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
                        
                        double rest_length = 1.0;
                        
                        if( (mCells.size()>0) && (!mIsGhostNode[p_element->GetNodeGlobalIndex(nodeA)]) && (!mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)]) )
                        {
                            double ageA = mCells[p_element->GetNode(nodeA)->GetIndex()].GetAge();
                            double ageB = mCells[p_element->GetNode(nodeB)->GetIndex()].GetAge();
                            if (ageA<1.0 && ageB<1.0 && fabs(ageA-ageB)<1e-6)
                            {
                                // Spring Rest Length Increases to normal rest length from 0.9 to normal rest length, 1.0, over 1 hour
                                rest_length=(0.1+0.9*ageA);
                                assert(rest_length<=1.0);
                               // std::cout<<p_element->GetNode(nodeA)->GetIndex()<<"\t"<<ageA<<"\t"<<p_simulation_time->GetDimensionalisedTime()<<std::endl;
                            }
                         }
                        
                        drdt_contribution = mpParams->GetMeinekeLambda() * unit_difference * (distance_between_nodes - rest_length) ;
                        
                        assert(!p_element->GetNode(nodeA)->IsDeleted());
                        assert(!p_element->GetNode(nodeB)->IsDeleted());
                        
                        bool AandBInLeftEdge = false;
	                    // Lookup whether node A and node B are both in a left edge element...
	                    if(mPeriodicSides)
	                    {
		                    for(unsigned i=0; i< mLeftCryptBoundary.size();i++)
		                    {
		                    	//std::cout << "i = " << i << "\n";
			                    if(mLeftCryptBoundary[i]==(unsigned)p_element->GetNode(nodeA)->GetIndex())
			                    {
			                    	for(unsigned j=0; j<mLeftCryptBoundary.size(); j++)
			                    	{
			                    		if(mLeftCryptBoundary[j]==(unsigned)p_element->GetNode(nodeB)->GetIndex())
			                    		{
			                    			AandBInLeftEdge = true;
			                    			//std::cout << "Left Boundary; Node A = " << mLeftCryptBoundary[i] << "\tNode B = "<< mLeftCryptBoundary[j] << "\n";
			                    			break;
			                    		}	
			                    	}
			                    }
			                }
	                    }

                        // Assume that if both nodes are real, or both are ghosts, then they both
                        // exert forces on each other, but if one is real and one is ghost then
                        // the real node exerts a force on the ghost node, but the ghost node 
                        // does NOT exert a force on the real node.  
                        if(!AandBInLeftEdge) // If A and B are in the left periodic edge ignore them (it will be handled by right edge)
                        {
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
            }
            
            ////////////////////////////////////////////////////////////////////////////////////////
            // Also loop over boundary edges so that all edges have been looped over exactly twice.
            ////////////////////////////////////////////////////////////////////////////////////////
            ConformingTetrahedralMesh<2,2>::BoundaryElementIterator elem_iter
               = mrMesh.GetBoundaryElementIteratorBegin();
            
            
            
            // this iterates over the outer edge elements (i.e. ghost nodes NOT real edge elements)
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

                    double rest_length = 1.0;
                    
                    if( (mCells.size()>0) &&  (!mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeA)]) && (!mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)]) )
                    {
                        double ageA = mCells[p_edge->GetNode(nodeA)->GetIndex()].GetAge();
                        double ageB = mCells[p_edge->GetNode(nodeB)->GetIndex()].GetAge();
                        if (ageA<1.0 && ageB<1.0 && fabs(ageA-ageB)<1e-6)
                        {
                           // Spring Rest Length Increases to normal rest length from 0.9 to normal rest length, 1.0, over 1 hour
                           #define COVERAGE_IGNORE
                           rest_length=(0.1+0.9*ageA);
                           assert(rest_length<=1.0);
                           #undef COVERAGE_IGNORE
                       }                          
                    }

                    drdt_contribution = mpParams->GetMeinekeLambda() * unit_difference * (distance_between_nodes - rest_length);
                       
                    assert(!p_edge->GetNode(nodeA)->IsDeleted());
                    assert(!p_edge->GetNode(nodeB)->IsDeleted());
                    


                    // Assume that if both nodes are real, or both are ghosts, then they both
                    // exert forces on each other, but if one is real and one is ghost then
                    // the real node exerts a force on the ghost node, but the ghost node 
                    // does NOT exert a force on the real node.   
                    if(mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeA)] == false)
                    {
	                    	// Real A force on any B
	                        drdt[ p_edge->GetNode(nodeB)->GetIndex()][0] -= drdt_contribution(0);
	                        drdt[ p_edge->GetNode(nodeB)->GetIndex()][1] -= drdt_contribution(1);
	
							// B exerts a force back if it is real.
	                        if(mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)] == false)
	                        {
	                            drdt[ p_edge->GetNode(nodeA)->GetIndex()][0] += drdt_contribution(0);
	                            drdt[ p_edge->GetNode(nodeA)->GetIndex()][1] += drdt_contribution(1);
	                        }
                    }
                    else
                    {
                    	// Ghost A receives a force
                        drdt[ p_edge->GetNode(nodeA)->GetIndex()][0] += drdt_contribution(0);
                        drdt[ p_edge->GetNode(nodeA)->GetIndex()][1] += drdt_contribution(1);
 
 						// Only a ghost B also receives a force
                        if(mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)] == true)
                        {
                           drdt[ p_edge->GetNode(nodeB)->GetIndex()][0] -= drdt_contribution(0);
                           drdt[ p_edge->GetNode(nodeB)->GetIndex()][1] -= drdt_contribution(1);
                        }
                    }                            
                }
                elem_iter++;
            }            

			//assert(0);            
            ///////////////
            //  sum forces for periodic nodes
            ////////////////
            if(mPeriodicSides)
            {
	            // Add up the forces on paired up nodes
	            // loop through size of left boundaries
	            for (unsigned i = 0; i < mLeftCryptBoundary.size();i++)
	            {
	            	double x_force = drdt[mLeftCryptBoundary[i]][0] + drdt[mRightCryptBoundary[i]][0];
	            	double y_force = drdt[mLeftCryptBoundary[i]][1] + drdt[mRightCryptBoundary[i]][1];
	            	drdt[mLeftCryptBoundary[i]][0] = x_force;
	            	drdt[mLeftCryptBoundary[i]][1] = y_force;
	            	drdt[mRightCryptBoundary[i]][0] = x_force;
	            	drdt[mRightCryptBoundary[i]][1] = y_force;
	            }
            }
            ////////////////////////////////////////////////////////////////////////////////////
            // update node positions
            ////////////////////////////////////////////////////////////////////////////////////

            
            for (int index = 0; index<mrMesh.GetNumAllNodes(); index++)
            {       
            	//std::cout << "Node " << index << "\t x_force = " << drdt[index][0] << "\t y_force = " << drdt[index][1] << "\n";
                     
                if(mFixedBoundaries)
                {
                    // All Boundaries x=0, x=crypt_width, y=0, y=crypt_length.
                    if(mrMesh.GetNode(index)->rGetLocation()[1]>0)
                    {
                        if(mrMesh.GetNode(index)->rGetLocation()[1]<mpParams->GetCryptLength())
                        {
                            if(mrMesh.GetNode(index)->rGetLocation()[0]>0)
                            {
                                if(mrMesh.GetNode(index)->rGetLocation()[0]<mpParams->GetCryptWidth())
                                {
                                    if(!mrMesh.GetNode(index)->IsDeleted())
                                    {
                                        //  std::cerr<<"Updating index "<<index<<"\n";
                                        c_vector<double,2> old_point = mrMesh.GetNode(index)->rGetLocation();
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
                else if(mCells.size()>0)
                {
                    // move any node as long as it is not a stem cell.
                    if(mCells[index].GetCellType()!=STEM)
                    {
                        if(!mrMesh.GetNode(index)->IsDeleted())
                        {
                            //  std::cerr<<"Updating index "<<index<<"\n";
                            c_vector<double,2> old_point = mrMesh.GetNode(index)->rGetLocation();
                            Point<2> new_point;

                            // note factor of 0.5 in the update because drdt was twice
                            // as large as it should be since edges were looped over twice.
                            new_point.rGetLocation()[0] = old_point[0] + 0.5*mDt*drdt[index][0]; // new_point_position[index];
                            new_point.rGetLocation()[1] = old_point[1] + 0.5*mDt*drdt[index][1]; // new_point_position[index];

                            // if a cell wants to move below y<0 (most likely because it was 
                            // just born from a stem cell), stop it doing so
                            if( (new_point.rGetLocation()[1] < 0.0) && (!mIsGhostNode[index]))
                            {
                            	// Here we give the cell a push upwards so that it doesn't get stuck on y=0 for ever.
                            	// it is a bit of a hack to make it work nicely!
                                #define COVERAGE_IGNORE
                                new_point.rGetLocation()[1] = 0.01;
                                #undef COVERAGE_IGNORE
                            }

                            mrMesh.SetNode(index, new_point, false);
                        }
                    }
                }
                else
                {
                    // no cells, just fix any node on line y=0
                    if(mrMesh.GetNode(index)->rGetLocation()[1]>0)
                    {
                        if(!mrMesh.GetNode(index)->IsDeleted())
                        {
                            //  std::cerr<<"Updating index "<<index<<"\n";
                            c_vector<double,2> old_point = mrMesh.GetNode(index)->rGetLocation();
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
            
            // Ensure no errors can creep in and move left nodes to same position as right ones
            if(mPeriodicSides)
            {
	            for (unsigned i = 0; i < mLeftCryptBoundary.size();i++)
	            {
	            	int RightNodeIndex = mRightCryptBoundary[i];
	            	int LeftNodeIndex = mLeftCryptBoundary[i];
	            	c_vector<double,2> right_point = mrMesh.GetNode(RightNodeIndex)->rGetLocation();
	            	Point<2> left_point;
	            	left_point.rGetLocation()[0] = right_point[0]-mpParams->GetCryptWidth();
	            	left_point.rGetLocation()[1] = right_point[1];
	            	mrMesh.SetNode(LeftNodeIndex, left_point, false);
	            	// Also force them to be the same cell
	            	// needed to synchronise cell cycle models (as R periodic cell cycle models are not run)...
	            	mCells[RightNodeIndex]=mCells[LeftNodeIndex];
	            }
            }
            //std::cout << "***************************************************\n";
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
                    //Node<2> *p_node=mrMesh.GetNode(0);
                    Node<2> *p_node = *it;
                    if(!p_node->IsDeleted())
                    {
                     double x = p_node->rGetLocation()[0];
                     double y = p_node->rGetLocation()[1];
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
                     	
                     	std::cout<<p_simulation_time->GetDimensionalisedTime() << "\t"<<sloughing_node_index<<"\t"<<target_node_index<<"\n";
                         // It's fallen off
                         assert(!mrMesh.GetNode(target_node_index)->IsDeleted());
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
            for (unsigned i=0; i<mrMesh.GetNumNodes(); i++)
            {
                Node<2> *p_node = mrMesh.GetNode(i);
                if(!p_node->IsDeleted())
                {
                    double x = p_node->rGetLocation()[0];
                    double y = p_node->rGetLocation()[1];
                    
                    double crypt_length=mpParams->GetCryptLength();
                    double crypt_width=mpParams->GetCryptWidth();
                    
                    if(!mPeriodicSides)
                    {
	                    if( (x>crypt_width) || (x<0.0) || (y>crypt_length))
	                    {
	                        mIsGhostNode[p_node->GetIndex()] = true;
	                        num_deaths++;
	                        //std::cout<< "num_deaths=" << num_deaths <<std::endl<< std::flush;
	                    }
                    }
                    else
                    {
						if(y>crypt_length)
	                    {
	                        mIsGhostNode[p_node->GetIndex()] = true;
	                        num_deaths++;
	                        if(mPeriodicSides)
	                        {
		                        // And delete the periodic image if appropriate.(don't count as a death since it is an image)
		                        for(unsigned j=0 ; j<mLeftCryptBoundary.size() ; j++)
		                        {
		                        	if((unsigned)p_node->GetIndex()==mLeftCryptBoundary[j])
		                        	{
		                        		mIsGhostNode[mRightCryptBoundary[j]]=true;
		                        	}	
		                        	if((unsigned)p_node->GetIndex()==mRightCryptBoundary[j])
		                        	{
		                        		mIsGhostNode[mLeftCryptBoundary[j]]=true;
		                        	}	
		                        }
	                        }
	                        //std::cout<< "num_deaths=" << num_deaths <<std::endl<< std::flush;
	                    }
                    }
                }
            }
            
            /*/////////////////////////////////////////////////////////
             * 
             * Designate cells as proliferating (transit) or
             * quiescent (differentiated) according to protein concentrations
             * 
             * If the betaCatenin level falls below a certain concentration then
             * the cell will (probably) be differentiated - if it later increases
             * it could become a transit cell again. This is just for visualization...
             * 
             *////////////////////////////////////////////////////////
            if(mWntIncluded)
            {	// Cycle through each cell
            	for (unsigned i=0; i<mrMesh.GetNumNodes(); i++)
            	{
            		if(!mIsGhostNode[i])
            		{
            			if(!(mCells[i].GetCellType()==STEM))
            			{
	            			// If we are in here the cell cycle model must be a WntCellCycleModel
	            			WntCellCycleModel *this_Wnt_model = static_cast<WntCellCycleModel*>(mCells[i].GetCellCycleModel());
	            			double betaCateninLevel = this_Wnt_model->GetProteinConcentrations()[6]+this_Wnt_model->GetProteinConcentrations()[7];
	            			//std::cout << "Cell " << i << ", beta-cat = " << betaCateninLevel << "\n" << std::endl;
	            			
	            			CryptCellType cell_type=TRANSIT;
	            			
	            			// For mitogenic stimulus of 6x10^-4 in Wnt equations
            				if(betaCateninLevel < 0.4127)
            				{
								cell_type = DIFFERENTIATED;
            				}
	            			// For mitogenic stimulus of 5x10^-4 in Wnt equations
//       //\todo get parameter right without breaking the build
//     				if(betaCateninLevel < 0.4954)
//            				{
//								cell_type = DIFFERENTIATED;
//            				}

            				
	           				mCells[i].SetCellType(cell_type);
            			}
            		}
            	}
            }
             
            
            if( mReMesh )
			{
    			ReMesh();
			}
            
            std::cout << " " << std::endl;
            
//            for (int i=0; i<mrMesh.GetNumAllNodes(); i++)
//            {
//            	std::cout<<"\t"<<mIsGhostNode[i]<<"\t"<<mIsGhostNode[map.GetNewIndex(i)]<<"\t"<<47*(mIsGhostNode[i]-mIsGhostNode[map.GetNewIndex(i)])<<"\n";
//            }
//            
            
            ////////////////////////////////////////////////////////////////////////////////
            // Write results to file
            ////////////////////////////////////////////////////////////////////////////////
            
            // Increment simulation time here, so results files look sensible
            p_simulation_time->IncrementTimeOneStep();

            double time = p_simulation_time->GetDimensionalisedTime();
            
            (*p_node_file) <<  time << "\t";
            (*p_element_file) <<  time << "\t";


            
            if(counter==0)
            {
                tabulated_node_writer.PutVariable(time_var_id, time);
                tabulated_element_writer.PutVariable(time_var_id_elem, time);
            }    
                
            
            /////////////////////////////////
            // write node files
            /////////////////////////////////
            for (int index = 0; index<mrMesh.GetNumAllNodes(); index++)
            {
                int colour = 0; // all green if no cells have been passed in
                
                if(mIsGhostNode[index]==true)
                {
                    colour = 4; // visualizer treats '4' these as invisible
                }
                else if(mCells.size()>0)
                {
                    if(index < (int) mCells.size())
                    {
                        CryptCellType type = mCells[index].GetCellType();
                        CryptCellMutationState mutation = mCells[index].GetMutationState();
                        
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
                        
                        if(mutation!=HEALTHY)
                        {
                        	colour = 3;	
                        }
                    }
                    else
                    {
                        #define COVERAGE_IGNORE
                        colour = 2; //Fix for segmentation fault
                        #undef COVERAGE_IGNORE
                    }
                    
                }
                
                if(!mrMesh.GetNode(index)->IsDeleted())
                {
                    const c_vector<double,2>& r_node_loc = mrMesh.GetNode(index)->rGetLocation();
                    (*p_node_file) << r_node_loc[0] << " "<< r_node_loc[1] << " " << colour << " ";

                    if(counter==0)
                    {   
                        tabulated_node_writer.PutVariable(x_position_var_ids[index], r_node_loc[0]);
                        tabulated_node_writer.PutVariable(y_position_var_ids[index], r_node_loc[1]);
                        tabulated_node_writer.PutVariable(type_var_ids[index], colour);
                    }
                }
            }
            (*p_node_file) << "\n";
            tabulated_node_writer.AdvanceAlongUnlimitedDimension();
            

            /////////////////////////////////
            // write element data files
            /////////////////////////////////
            for (int elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
            {
                if (!mrMesh.GetElement(elem_index)->IsDeleted())
                {
                    (*p_element_file) << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0)<< " " << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1)<< " "<< mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2)<< " ";
                    
                    if(counter==0)
                    {
                        tabulated_element_writer.PutVariable(nodeA_var_ids[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0));
                        tabulated_element_writer.PutVariable(nodeB_var_ids[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1));
                        tabulated_element_writer.PutVariable(nodeC_var_ids[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2));
                    }
                }
            }
            (*p_element_file) << "\n";
            tabulated_element_writer.AdvanceAlongUnlimitedDimension();

            counter++; 
            if(counter > 80)
            {
                counter = 0;
            }
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
    	// First set all to not be ghost nodes
    	for (unsigned i=0 ; i<mIsGhostNode.size() ; i++)
    	{
    		mIsGhostNode[i] = false;
    	}
    	// then update which ones are.
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
    
    /**
     * Set the simulation to run with no birth.
     */   
    void SetNoBirth(bool nobirth)
    {
        mNoBirth = nobirth;
    }
    
    
    /**
     * Method to calculate the boundary of the crypt within the whole mesh ie the interface 
     * between normal and ghost nodes.
     */    
    void CalculateCryptBoundary()
    {   
    	assert(mPeriodicSides);

		double crypt_width=mpParams->GetCryptWidth();
               
		// Set all nodes as not being on the boundary...
        std::vector<bool> is_nodes_on_boundary(mIsGhostNode.size());
        for(unsigned i = 0 ; i < is_nodes_on_boundary.size() ; i++)
        {
            is_nodes_on_boundary[i]=false;
            mIsPeriodicNode[i]=false;
		}
       
        // Loop over elements and find boundary nodes of crypt
        for (int elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
        {
            Element<2,2>* p_element = mrMesh.GetElement(elem_index);
    
            if (!mIsGhostNode[p_element->GetNode(0)->GetIndex()])
            {
                if((mIsGhostNode[p_element->GetNode(1)->GetIndex()]) || (mIsGhostNode[p_element->GetNode(2)->GetIndex()])   )
                {
                    is_nodes_on_boundary[p_element->GetNode(0)->GetIndex()]= true;
                }
            }
             
            if (!mIsGhostNode[p_element->GetNode(1)->GetIndex()])
            {
                
                if((mIsGhostNode[p_element->GetNode(0)->GetIndex()]) || (mIsGhostNode[p_element->GetNode(2)->GetIndex()])   )
                {
                    is_nodes_on_boundary[p_element->GetNode(1)->GetIndex()]= true;
                }
            }
             
            if (!mIsGhostNode[p_element->GetNode(2)->GetIndex()])
            {
                if((mIsGhostNode[p_element->GetNode(1)->GetIndex()]) || (mIsGhostNode[p_element->GetNode(0)->GetIndex()]))
                {
                    is_nodes_on_boundary[p_element->GetNode(2)->GetIndex()]= true; 
                }
            }
        }
                  
        std::vector<unsigned> nodes_on_boundary;   
        std::vector<unsigned> nodes_on_left_boundary;
        std::vector<unsigned> nodes_on_right_boundary;

        
        for(unsigned i = 0; i < is_nodes_on_boundary.size(); i++) 
        {
        
            if(is_nodes_on_boundary[i])
            {
                nodes_on_boundary.push_back(i);
            }
        }
		
		
		
        for(unsigned i=0; i<nodes_on_boundary.size(); i++)
        {
        	for(unsigned j=0; j<nodes_on_boundary.size(); j++)
        	{
        		// Check y positions are the same
        		 if(fabs(mrMesh.GetNode(nodes_on_boundary[i])->rGetLocation()[1]-mrMesh.GetNode(nodes_on_boundary[j])->rGetLocation()[1])<1e-3)
        	     {
        	     	// Check x positions are crypt width apart.
        	     	if(fabs(mrMesh.GetNode(nodes_on_boundary[j])->rGetLocation()[0]-mrMesh.GetNode(nodes_on_boundary[i])->rGetLocation()[0]-crypt_width)<1e-3)
        	     	{
        	     		nodes_on_left_boundary.push_back(nodes_on_boundary[i]);
        	     		nodes_on_right_boundary.push_back(nodes_on_boundary[j]);
        	     		mIsPeriodicNode[nodes_on_boundary[i]]=true;
        	     		mIsPeriodicNode[nodes_on_boundary[j]]=true;
        	     	}
        	     }
        	}
        }
        
		mLeftCryptBoundary =  nodes_on_left_boundary;
        mRightCryptBoundary =  nodes_on_right_boundary;
        mCryptBoundary =  nodes_on_boundary; 
		//std::cout << "Total boundary nodes = " << nodes_on_boundary.size() << " left boundary nodes = " << nodes_on_left_boundary.size() << "\n";
    }
    
    /**
     * This function detects when a remesh has caused a cell to 
     * join or leave one of the periodic boundaries. It figures out where an image cell
     * will have to be placed to match it on the other side, and which two
     * periodic nodes have been upset.
     */   
    void DetectNaughtyCellsAtPeriodicEdges()
    {
        assert(mPeriodicSides);

		//std::cout << "Nodes on Boundary : \n";
    	// Check for any surplus cells to remove...
		RemoveSurplusCellsFromPeriodicBoundary();
		
		/* For each node see if it has broken onto a periodic boundary
    	 * If it has create an image node at the other side
    	 */
    	unsigned i = 0;
    	while(i<mCryptBoundary.size())
        {
        	std::vector< unsigned > nodes_on_left_boundary = mLeftCryptBoundary;
        	std::vector< unsigned > nodes_on_right_boundary = mRightCryptBoundary;
        	std::vector< unsigned > nodes_on_boundary = mCryptBoundary;
	        unsigned our_node = nodes_on_boundary[i];
	        //std::cout << our_node << "\n";
	        unsigned number_of_left_periodic_neighbours = 0;
	        unsigned number_of_right_periodic_neighbours = 0;
	        bool our_node_periodic = false;
	        
	        std::vector <unsigned > periodic;
	        
	        for (unsigned j = 0 ; j<nodes_on_left_boundary.size() ; j++)
	        {
	        	if(our_node == nodes_on_left_boundary[j] || our_node == nodes_on_right_boundary[j])
	        	{
	        		our_node_periodic = true;
                    break;
	        	}
	        }
	        
	        if(!our_node_periodic)
	        {
	        	for (int elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
		        {
		            Element<2,2>* p_element = mrMesh.GetElement(elem_index);
		    		unsigned node[3];
		    		node[0] = (unsigned)p_element->GetNode(0)->GetIndex();
		    		node[1] = (unsigned)p_element->GetNode(1)->GetIndex();
		    		node[2] = (unsigned)p_element->GetNode(2)->GetIndex();
		    		
		    		unsigned periodic_left_nodes_in_this_element = 0;
		    		unsigned periodic_right_nodes_in_this_element = 0;
		    		
		        	if(node[0]==our_node || node[1]==our_node || node[2]==our_node)
		        	{
		        		if(mIsGhostNode[node[0]] || mIsGhostNode[node[1]] || mIsGhostNode[node[2]])
			    		{	// If this is one of the new nodes (or top or bottom edges)
							for(unsigned j = 0 ; j<nodes_on_left_boundary.size() ; j++)
							{	//Cycle through periodic nodes and check for periodic neighbours
			    				unsigned periodic_node = nodes_on_left_boundary[j];
								if(node[0]==periodic_node || node[1]==periodic_node || node[2]==periodic_node)
								{
									periodic_left_nodes_in_this_element++;
									periodic.push_back(periodic_node);
								}
								periodic_node = nodes_on_right_boundary[j];
								if(node[0]==periodic_node || node[1]==periodic_node || node[2]==periodic_node)
								{
									periodic_right_nodes_in_this_element++;
									periodic.push_back(periodic_node);
								}
			    			}
			    		}	
		        	}
		        	
					//If both periodic nodes are in this element we are at 
					//the top or bottom corner and shouldn't count it.
					if(periodic_left_nodes_in_this_element==2)
					{
                        #define COVERAGE_IGNORE
						periodic_left_nodes_in_this_element=0;	
                        #undef COVERAGE_IGNORE
					}
					if(periodic_right_nodes_in_this_element==2)
					{
                        #define COVERAGE_IGNORE
						periodic_right_nodes_in_this_element=0;	
                        #undef COVERAGE_IGNORE
					}
					number_of_right_periodic_neighbours += periodic_right_nodes_in_this_element;
					number_of_left_periodic_neighbours += periodic_left_nodes_in_this_element;
				}// next element
			
				// Check we have actually broken into an edge that is periodic and 
				// just sensed it from the top edge i.e. the two periodic neighbours are different.
//				if(periodic.size()>2)
//                {
//                    EXCEPTION("Gary says this should not happen");   
//                }
                if(periodic.size()==2)
                {
                    if(periodic[0] == periodic[1])
    				{
						#define COVERAGE_IGNORE
    					number_of_left_periodic_neighbours=1;
    					number_of_right_periodic_neighbours=1;
						#undef COVERAGE_IGNORE
    				}
                }
                
		        //If there are no periodic neighbours then we are at the bottom/top edge
                //If there is one perioic neighbour then we are at the corner?
                if(number_of_left_periodic_neighbours==2)
		        {
                    #define COVERAGE_IGNORE
		        	// We should have a new periodic node
		        	double old_x = mrMesh.GetNode(our_node)->rGetLocation()[0];
		        	double old_y = mrMesh.GetNode(our_node)->rGetLocation()[1];
		        	double crypt_width = mpParams->GetCryptWidth();
		        	std::cout << "LEFT Node "<< our_node << " has broken into the periodic edge between nodes\n";for(unsigned k=0 ; k<periodic.size() ; k++)
					{
						std::cout << periodic[k] << "\t";
					}
					std::cout << "\n";
					AddACellToPeriodicBoundary(our_node,old_x+crypt_width,old_y,periodic);
					//NodeMap map(mrMesh.GetNumAllNodes());
	    			//mrMesh.ReMesh(map);
					CalculateCryptBoundary();
					RemoveSurplusCellsFromPeriodicBoundary();
                    #undef COVERAGE_IGNORE
				}
		        
		        if(number_of_right_periodic_neighbours==2)
		        {
		        	// We should have a new periodic node
		        	double old_x = mrMesh.GetNode(our_node)->rGetLocation()[0];
		        	double old_y = mrMesh.GetNode(our_node)->rGetLocation()[1];
		        	double crypt_width = mpParams->GetCryptWidth();
		        	std::cout << "RIGHT Node "<< our_node << " has broken into the periodic edge between nodes\n";
					for(unsigned k=0 ; k<periodic.size() ; k++)
					{
						std::cout << periodic[k] << "\t";
					}
					std::cout << "\n";
					AddACellToPeriodicBoundary(our_node,old_x-crypt_width,old_y,periodic);
					//NodeMap map(mrMesh.GetNumAllNodes());
	    			//mrMesh.ReMesh(map);	// it is possible that the only ghost node a
	    			// periodic node was attached to is now real and so we need to remesh so that
	    			// the crypt boundary can still be recognised.
					CalculateCryptBoundary();
					RemoveSurplusCellsFromPeriodicBoundary();
				}
			}// end of if(!our_node_periodic)
            
            i++;
		}// next node on boundary.
	}// end of function
	
	/**
	 * For each of the old periodic boundary nodes check they are still 
	 * on the boundary, if not delete their image.
	 */
	void RemoveSurplusCellsFromPeriodicBoundary()
	{
        assert(mPeriodicSides);

		for(unsigned i=0 ; i<mOldLeftCryptBoundary.size() ; i++)
		{
			bool this_left_node_missing = true;
			bool this_right_node_missing = true;
			
			unsigned old_node_on_left_boundary = mOldLeftCryptBoundary[i];
			unsigned old_node_on_right_boundary = mOldRightCryptBoundary[i];
			
			for(unsigned j=0 ; j<mCryptBoundary.size() ; j++)
			{// search through the new crypt boundaries for this node
				if(old_node_on_left_boundary==mCryptBoundary[j])
				{
					this_left_node_missing=false;
				}
				if(old_node_on_right_boundary==mCryptBoundary[j])
				{
					this_right_node_missing=false;
				}
			}
			
			if(this_left_node_missing && (!mIsGhostNode[old_node_on_left_boundary]))
			{	// The left node has been internalised (was periodic and is not a ghost) so the right node should be spooked
				mIsGhostNode[old_node_on_right_boundary]=true;
				std::cout << "Periodic Node "<< old_node_on_left_boundary <<" internalised \n Right Node " << old_node_on_right_boundary << " spooked\n";
				CalculateCryptBoundary();

			}
			if(this_right_node_missing && (!mIsGhostNode[old_node_on_right_boundary]))
			{	// The right node has been internalised (was periodic and is not a ghost) so the right node should be spooked
				mIsGhostNode[old_node_on_left_boundary]=true;
				std::cout << "Periodic Node " << old_node_on_right_boundary << " internalised \n Left Node " << old_node_on_left_boundary << " spooked\n";
				CalculateCryptBoundary();
			}
		}
		// Update the history vectors now we have used them in case we have to use them again this time step.
		mOldLeftCryptBoundary = mLeftCryptBoundary;
		mOldRightCryptBoundary = mRightCryptBoundary;
		mOldCryptBoundary = mCryptBoundary;
	}
	
	
    /**
     * This function adds in a periodic image cell by moving the ghost node 
     * in the element attached to the upset periodic nodes' images to the 
     * correct location and makes it real. (copying cell info from the original 
     * offending node.)
     */
    void AddACellToPeriodicBoundary(unsigned original_node_index, double new_x, double new_y, std::vector< unsigned > periodic)
    {
    	assert(mPeriodicSides);

        unsigned node_index=0;        
		//std::cout << "Periodic node " << original_node_index<< " should have an image created\n";
        
        
        Point<2> new_point;
        new_point.SetCoordinate(0, new_x);
        new_point.SetCoordinate(1, new_y);
        //std::cout << "New point to be introduced at x = " << new_x << " y = " << new_y << "\n";
        
        
        std::vector <unsigned > periodic_nodes;
        periodic_nodes.reserve(2);
        
        //Find periodic nodes in the boundary lists
        for (unsigned i=0 ; i<mLeftCryptBoundary.size() ; i++)
        {
        	for (unsigned j=0 ; j<2 ; j++)
        	{
	        	if (periodic[j]==mLeftCryptBoundary[i])
	        	{
                    #define COVERAGE_IGNORE
        			periodic_nodes[j] = mRightCryptBoundary[i];
                    #undef COVERAGE_IGNORE
	        	}
	        	if (periodic[j]==mRightCryptBoundary[i])
        		{
        			periodic_nodes[j] = mLeftCryptBoundary[i];
        		}
        	}
        }
        
        //std::cout << "Between cells " << periodic_nodes[0] << " and " << periodic_nodes[1] << "\n";
        
        std::vector <unsigned > ghosts_on_node_0;
        std::vector <unsigned > ghosts_on_node_1;
        // Search all the elements connected to periodic_node[0] for ghost nodes
        for (int elim_index = 0 ; elim_index<mrMesh.GetNumAllElements(); elim_index++ )
        {
        	Element<2,2>* p_element = mrMesh.GetElement(elim_index);
        	unsigned node[3];
        	node[0] = (unsigned)p_element->GetNode(0)->GetIndex();
		    node[1] = (unsigned)p_element->GetNode(1)->GetIndex();
		    node[2] = (unsigned)p_element->GetNode(2)->GetIndex();
        	if	(node[0]==periodic_nodes[0] || node[1]==periodic_nodes[0] || node[2]==periodic_nodes[0])
        	{	// if this element contains our first target node
    			for(unsigned i=0 ; i<3 ;i++)
    			{
    				if(mIsGhostNode[node[i]])
    				{
    					bool already_flagged = false;
		    			for(unsigned j=0; j<ghosts_on_node_0.size() ; j++)
		    			{
		    				if(ghosts_on_node_0[j]==node[i])
		    				{
		    					already_flagged = true;
		    				}
		    			}
		    			if(!already_flagged)
		    			{
			    			ghosts_on_node_0.push_back(node[i]);	
		    			}
					}
    			}
        	}
        	if	(node[0]==periodic_nodes[1] || node[1]==periodic_nodes[1] || node[2]==periodic_nodes[1])
        	{	// if this element contains our first target node
    			for(unsigned i=0 ; i<3 ;i++)
    			{
    				if(mIsGhostNode[node[i]])
    				{
    					bool already_flagged = false;
		    			for(unsigned j=0; j<ghosts_on_node_1.size() ; j++)
		    			{
		    				if(ghosts_on_node_1[j]==node[i])
		    				{
		    					already_flagged = true;
		    				}
		    			}
		    			if(!already_flagged)
		    			{
			    			ghosts_on_node_1.push_back(node[i]);	
		    			}
					}
    			}
        	}
        }
        bool image_created = false;
        for(unsigned i=0 ; i< ghosts_on_node_0.size() ; i++)
        {
        	for(unsigned j=0 ; j< ghosts_on_node_1.size() ; j++)	
        	{
        		if(ghosts_on_node_0[i] == ghosts_on_node_1[j])	
        		{
        			//std::cout << "Shared Node is " << ghosts_on_node_0[i] << "\n";	
        			// Move the ghost node to the correct position
        			node_index = ghosts_on_node_0[i];
					image_created = true;
        		}
        	}
        }
        if(image_created==false)
        {	
            // There are no shared ghost nodes - this could happen if two cells have
        	// broken into boundary at the same time. Make nearest ghost node a real one.
        	#define COVERAGE_IGNORE
            unsigned nearest_node = 0;
        	double nearest_distance = 100.0;
        	for (unsigned i=0 ; i<ghosts_on_node_0.size() ; i++)
        	{
        		Node<2> *p_our_node = mrMesh.GetNode(ghosts_on_node_0[i]);
	            double x = p_our_node->GetPoint()[0];
	            double y = p_our_node->GetPoint()[1];
	            double this_ghost_distance = sqrt((x-new_x)*(x-new_x)+(y-new_y)*(y-new_y));
	            if(this_ghost_distance<nearest_distance)
	            {
	            	nearest_distance=this_ghost_distance;
	            	nearest_node = ghosts_on_node_0[i];	
	            }
        	}
        	for (unsigned i=0 ; i<ghosts_on_node_1.size() ; i++)
        	{
        		Node<2> *p_our_node = mrMesh.GetNode(ghosts_on_node_1[i]);
	            double x = p_our_node->GetPoint()[0];
	            double y = p_our_node->GetPoint()[1];
	            double this_ghost_distance = sqrt((x-new_x)*(x-new_x)+(y-new_y)*(y-new_y));
	            if(this_ghost_distance<nearest_distance)
	            {
	            	nearest_distance=this_ghost_distance;
	            	nearest_node = ghosts_on_node_1[i];	
	            }
        	}
        	node_index = nearest_node;
        	image_created=true;
            #define COVERAGE_IGNORE
        }
        
        if(image_created==false)
        {
            #define COVERAGE_IGNORE
        	EXCEPTION("No ghost nodes to form image cell");
            #undef COVERAGE_IGNORE
        }
        
        mrMesh.SetNode(node_index, new_point, false);
		std::cout << "New image cell created from ghost node " << node_index<< "\n ";
		// Stop it being a ghost node
		mIsGhostNode[node_index] = false;
		// copy relevant cell info across...
		mCells[node_index]=mCells[original_node_index];
        
	}
    
    void ReMesh()
    {
    	NodeMap map(mrMesh.GetNumAllNodes());
    	mrMesh.ReMesh(map);

		if(mPeriodicSides)
		{
			mOldLeftCryptBoundary = mLeftCryptBoundary;
			mOldRightCryptBoundary = mRightCryptBoundary;
			mOldCryptBoundary = mCryptBoundary;
		    CalculateCryptBoundary();
    	
			if(mOldCryptBoundary!=mCryptBoundary)
	        {	// There has been a change in the boundary
	        	std::cout << "Boundary has changed\n";
	        	DetectNaughtyCellsAtPeriodicEdges();
			}
		}
	}
};

#endif /*CRYPTSIMULATION2DPERIODIC_HPP_*/
