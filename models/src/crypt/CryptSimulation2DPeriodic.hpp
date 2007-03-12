#ifndef CRYPTSIMULATION2DPERIODIC_HPP_
#define CRYPTSIMULATION2DPERIODIC_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp> // for archiving vectors
#include <boost/serialization/string.hpp>

#include <fstream>

#include "ConformingTetrahedralMesh.cpp"
#include "MeinekeCryptCell.hpp"
#include "CancerParameters.hpp"
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
//TODO: This should become abstract
#include "RandomCellKiller.hpp"
#include "OutputFileHandler.hpp"

/**
 * Structure encapsulating variable identifiers for the node datawriter
 */
typedef struct node_writer_ids_t
{
    unsigned time;               /**< The simulation time */
    std::vector<unsigned> types; /**< Cell types */
    /** Cell positions */
    std::vector<unsigned> x_positions, y_positions;
} node_writer_ids_t;

/**
 * Structure encapsulating variable identifiers for the element datawriter
 */
typedef struct element_writer_ids_t
{
    unsigned time;/**< The simulation time */
    /** Node indices */
    std::vector<unsigned> nodeAs, nodeBs, nodeCs;
} element_writer_ids_t;


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
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2DPeriodic;
private:
    double mDt;
    double mEndTime;
    ConformingTetrahedralMesh<2,2> &mrMesh;
    
    bool mIncludeRandomBirth;
    bool mIncludeVariableRestLength;

    /** Whether to fix all four boundaries (defaults to false).*/  
    bool mFixedBoundaries;    
    
    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;

    /** Whether to remesh at each timestep or not (defaults to true).*/    
    bool mReMesh;
    
    /** Whether the remeshing has made our periodic handlers do anything and
     *  whether it is worth remeshing again (defaults to false).*/    
	bool mNodesMoved;
    
    /** Whether the mesh is periodic or not (defaults to false).*/    
    bool mPeriodicSides;

    /** Whether each node is ghosted-ified or not.*/  
    std::vector <bool> mIsGhostNode; 
    std::vector <bool> mIsPeriodicNode;
    
    /** The node indexes of nodes on the left boundary. */ 
    std::vector <unsigned> mLeftCryptBoundary;
    std::vector <unsigned> mOldLeftCryptBoundary;
    /** The node indexes of nodes on the right boundary. */
    std::vector <unsigned> mRightCryptBoundary;
    std::vector <unsigned> mOldRightCryptBoundary;
    /** The node indexes of nodes on the whole boundary. */
    std::vector <unsigned> mCryptBoundary;
    std::vector <unsigned> mOldCryptBoundary;
    
    /** The maximum number of cells that this simulation will include (for use by datawriter). */
    unsigned mMaxCells;
    /** The maximum number of elements that this simulation will include (for use by datawriter). */
    unsigned mMaxElements;
    
    std::string mOutputDirectory;
    /** Every cell in the simulation*/
    std::vector<MeinekeCryptCell> mCells;

    RandomNumberGenerator *mpRandomNumberGenerator;    
    /** Whether the random number generator was created by our constructor*/
    bool mCreatedRng;
    
    /** The Meineke and cancer parameters */
    CancerParameters *mpParams;
    
    /** Whether Wnt signalling is included or not (defaults to false).*/    
    bool mWntIncluded;
    /** The Wnt gradient, if any */
    WntGradient mWntGradient;
    
    /** Number of remeshes performed in the current time step */
    unsigned mRemeshesThisTimeStep;

    /** Counts the number of births during the simulation */
    unsigned mNumBirths;
    
    /** Counts the number of deaths during the simulation */
    unsigned mNumDeaths;
    
    /**
     * Prevents multiple near-simultaneous cell divisions on the periodic boundary -
     * once one cell has divided, other divisions are postponed for a couple of 
     * timesteps, until this counter reaches 0.  This is to cope with re-meshing 
     * issues; would be nice to get rid of it eventually.
     */
    unsigned mPeriodicDivisionBuffer;
        
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        mpParams = CancerParameters::Instance();
        archive & *mpParams;
        archive & mpParams;
        
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mDt;
        archive & mEndTime; 
        archive & mIncludeRandomBirth;
        archive & mIncludeVariableRestLength;
        archive & mFixedBoundaries;
        archive & mNoBirth;
        archive & mReMesh;
        archive & mNodesMoved;
        archive & mPeriodicSides;
        archive & mIsGhostNode;
        archive & mIsPeriodicNode;
        archive & mLeftCryptBoundary;
        archive & mOldLeftCryptBoundary;
        archive & mRightCryptBoundary;
        archive & mOldRightCryptBoundary;
        archive & mCryptBoundary;
        archive & mOldCryptBoundary;
        archive & mMaxCells;
        archive & mMaxElements;
//        archive & mOutputDirectory;
        archive & mCells;
//        archive & mpRandomNumberGenerator; 
        archive & mCreatedRng;
        archive & mWntIncluded;
        archive & mWntGradient;
        archive & mRemeshesThisTimeStep;
        archive & mNumBirths;
        archive & mNumDeaths;
        archive & mPeriodicDivisionBuffer;
    }
   
    
    
    /** Cell killer */
    //TODO: Should become an abstract cell killer
    RandomCellKiller<2> *mpCellKiller;
    
    /**
     * Define the variable identifiers in the data writer used to write node-based results.
     * 
     * Uses mMaxCells to decide how many variables to define.
     */
    void SetupNodeWriter(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rVarIds)
    {
        rVarIds.time = rNodeWriter.DefineUnlimitedDimension("Time","hours");
        
        rVarIds.types.resize(mMaxCells);
        rVarIds.x_positions.resize(mMaxCells);
        rVarIds.y_positions.resize(mMaxCells);

        // set up per-cell variables
        for (unsigned cell=0; cell<mMaxCells; cell++)
        {
            std::stringstream cell_type_var_name, cell_x_position_var_name, cell_y_position_var_name;
            cell_type_var_name << "cell_type_" << cell;
            cell_x_position_var_name << "cell_x_position_" << cell;
            cell_y_position_var_name << "cell_y_position_" << cell;            
            rVarIds.types[cell]=rNodeWriter.DefineVariable(cell_type_var_name.str(),"dimensionless");
            rVarIds.x_positions[cell]=rNodeWriter.DefineVariable(cell_x_position_var_name.str(),"rest_spring_length");
            rVarIds.y_positions[cell]=rNodeWriter.DefineVariable(cell_y_position_var_name.str(),"rest_spring_length");
        }
        
        rNodeWriter.EndDefineMode();
    }
    
    /**
     * Define the variable identifiers in the data writer used to write element-based results.
     * 
     * Uses mMaxCells to decide how many variables to define.
     */
    void SetupElementWriter(ColumnDataWriter& rElementWriter, element_writer_ids_t& rVarIds)
    {
        rVarIds.time = rElementWriter.DefineUnlimitedDimension("Time","hours");
        
        // Set up columns for element writer      
        rVarIds.nodeAs.resize(mMaxElements);
        rVarIds.nodeBs.resize(mMaxElements);
        rVarIds.nodeCs.resize(mMaxElements);
        
        for (unsigned elem_index = 0; elem_index<mMaxElements; elem_index++)
        {
            std::stringstream nodeA_var_name;
            std::stringstream nodeB_var_name;
            std::stringstream nodeC_var_name;
            
            nodeA_var_name << "nodeA_" << elem_index;
            nodeB_var_name << "nodeB_" << elem_index;
            nodeC_var_name << "nodeC_" << elem_index;
            
            rVarIds.nodeAs[elem_index] = rElementWriter.DefineVariable(nodeA_var_name.str(),"dimensionless");
            rVarIds.nodeBs[elem_index] = rElementWriter.DefineVariable(nodeB_var_name.str(),"dimensionless");
            rVarIds.nodeCs[elem_index] = rElementWriter.DefineVariable(nodeC_var_name.str(),"dimensionless");
        }
        
        rElementWriter.EndDefineMode();
    }
    
    
    void WriteResultsToFiles(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rNodeVarIds,
                             ColumnDataWriter& rElementWriter, element_writer_ids_t& rElementVarIds,
                             std::ofstream& rNodeFile, std::ofstream& rElementFile,
                             bool writeTabulatedResults,
                             bool writeVisualizerResults)
    {
        // Write current simulation time
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        double time = p_simulation_time->GetDimensionalisedTime();
        
        if(writeVisualizerResults)
        {
            rNodeFile <<  time << "\t";
            rElementFile <<  time << "\t";
        }
        
        if (writeTabulatedResults)
        {
            rNodeWriter.PutVariable(rNodeVarIds.time, time);
            rElementWriter.PutVariable(rElementVarIds.time, time);
        }
        
        /////////////////////////////////
        // write node files
        /////////////////////////////////
        for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
        {
            if(index>mMaxCells)
            {
                #define COVERAGE_IGNORE
                EXCEPTION("\nNumber of cells exceeds mMaxCells. Use SetMaxCells(unsigned) to increase it.\n");
                #undef COVERAGE_IGNORE
            }
            unsigned colour = 0; // all green if no cells have been passed in
            
            if(mIsGhostNode[index]==true)
            {
                colour = 4; // visualizer treats '4' these as invisible
            }
            else if(mCells.size()>0)
            {
                if(index < mCells.size())
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
                if(writeVisualizerResults)
                {
                    rNodeFile << r_node_loc[0] << " "<< r_node_loc[1] << " " << colour << " ";
                }
                if(writeTabulatedResults)
                {   
                    rNodeWriter.PutVariable(rNodeVarIds.x_positions[index], r_node_loc[0]);
                    rNodeWriter.PutVariable(rNodeVarIds.y_positions[index], r_node_loc[1]);
                    rNodeWriter.PutVariable(rNodeVarIds.types[index], colour);
                }
            }
        }

        /////////////////////////////////
        // write element data files
        /////////////////////////////////
        for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
        {
            if (elem_index>mMaxElements)
            {
                #define COVERAGE_IGNORE
                EXCEPTION("Maximum number of elements (mMaxElements) exceeded.\nUse SetMaxElements(unsigned) to increase it.\n");
                #undef COVERAGE_IGNORE
            }
            if (!mrMesh.GetElement(elem_index)->IsDeleted())
            {
                if(writeVisualizerResults)
                {
                    rElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0)<< " " << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1)<< " "<< mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2)<< " ";
                }
                if(writeTabulatedResults)
                {
                    rElementWriter.PutVariable(rElementVarIds.nodeAs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0));
                    rElementWriter.PutVariable(rElementVarIds.nodeBs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1));
                    rElementWriter.PutVariable(rElementVarIds.nodeCs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2));
                }
            }
        }

        if(writeVisualizerResults)
        {
            rNodeFile << "\n";
            rElementFile << "\n";
        }
        if(writeTabulatedResults)
        {
            rNodeWriter.AdvanceAlongUnlimitedDimension();
            rElementWriter.AdvanceAlongUnlimitedDimension();
        }   
    }
    
    /**
     * During a simulation time step, process any cell divisions that need to occur.
     * If the simulation includes cell birth, causes (almost) all cells that are ready to divide
     * to produce daughter cells.
     * 
     * Some divisions are postponed to work around issues with remeshing.  Specifically, only
     * one division on the periodic boundary is allowed within any 3 consecutive time steps.
     * This is controlled by the mPeriodicDivisionBuffer attribute.
     * 
     * @return the number of births that occurred.
     */
    unsigned DoCellBirth()
    {
        unsigned num_births = 0;
        if(mPeriodicDivisionBuffer>0)
        {
            mPeriodicDivisionBuffer--;
        }
        if (!mNoBirth && !mCells.empty())
        {
            // Iterate over all cells, seeing if each one can be divided
            for (unsigned i=0; i<mCells.size(); i++)
            {
                bool skip = false; // Whether to not try dividing this cell
                if(mrMesh.GetNode(i)->IsDeleted()) skip=true; // Skip deleted cells
                if(mIsGhostNode[i]) skip=true; // Skip Ghost nodes
                bool periodic_cell = false;
                unsigned periodic_index = 0; // Index of the periodic cell in the boundary, as opposed to in the mesh

                // Check if this cell is on the periodic boundary; there are more conditions if it is
                if(!skip && mPeriodicSides)
                {
                    for (unsigned j=0 ; j < mRightCryptBoundary.size() ; j++)
                    {
                        if (mRightCryptBoundary[j]==i)
                        {   // Only allow one periodic boundary to have divisions...
                            skip=true;
                        }
                        if (mLeftCryptBoundary[j]==i)
                        {
                            if(mPeriodicDivisionBuffer>0)
                            {
                                // Only allow one periodic cell division per 
                                // timestep so that mesh can catch up with it.
                                // it will divide next timestep anyway
                                skip=true;  
                            } 
                            else 
                            {
                                periodic_cell = true;
                                periodic_index = j;
                            }
                        }
                    }
                }
                if(skip) continue;
                
                // Check for this cell dividing
                //std::cout << "On cell "<< i << std::endl;
                // Construct any influences for the cell cycle...
                Node<2> *p_our_node = mrMesh.GetNode(i);
                std::vector<double> cell_cycle_influences;
                if(mWntIncluded)
                {
                    double y = p_our_node->rGetLocation()[1];
                    double wnt_stimulus = mWntGradient.GetWntLevel(y);
                    cell_cycle_influences.push_back(wnt_stimulus);
                }
                
                // CHECK if this cell is ready to divide - if so create a new cell etc.
                if(mCells[i].ReadyToDivide(cell_cycle_influences))
                {
                    // Create new cell
                    MeinekeCryptCell new_cell = mCells[i].Divide();
                    if(mPeriodicSides && periodic_cell)
                    {   
                        std::cout << "Periodic Division\n";
                        mPeriodicDivisionBuffer=3;
                        //Make sure the image cell knows it has just divided and aged a generation
                        mCells[mRightCryptBoundary[periodic_index]]=mCells[mLeftCryptBoundary[periodic_index]];
                    }
                    else
                    {
                        std::cout << "Cell division at node " << i << "\n";
                    }

                    // Add new node to mesh
                    Element<2,2>* p_element = FindElementForBirth(p_our_node, i,
                                                                  periodic_cell, periodic_index);
                    
                    //std::cout << "New cell being intoduced into element with nodes \n";
                    //std::cout << p_element->GetNodeGlobalIndex(0) << "\t" << p_element->GetNodeGlobalIndex(1) << "\t" <<p_element->GetNodeGlobalIndex(2) << "\n";
                    double x = p_our_node->rGetLocation()[0];
                    double y = p_our_node->rGetLocation()[1];
                        
                    double x_centroid = (1.0/3.0)*(p_element->GetNode(0)->rGetLocation()[0]
                                                    +  p_element->GetNode(1)->rGetLocation()[0]
                                                    +  p_element->GetNode(2)->rGetLocation()[0] );
                                                    
                    double y_centroid = (1.0/3.0)*(p_element->GetNode(0)->rGetLocation()[1]
                                                    +  p_element->GetNode(1)->rGetLocation()[1]
                                                    +  p_element->GetNode(2)->rGetLocation()[1] );
                                                    
                    
                    // check the new point is in the triangle
                    double distance_from_node_to_centroid =  sqrt(  (x_centroid - x)*(x_centroid - x)
                                                                  + (y_centroid - y)*(y_centroid - y) );
                    
                    // we assume the new cell is a distance 0.1 away from the old.
                    // however, to avoid crashing in usual situations we check this
                    // new position is actually in the triangle being refined.
                    // TODO: Check this is correct!
                    double distance_of_new_cell_from_parent = 0.1;
                    if(distance_from_node_to_centroid < (2.0/3.0)*0.1)
                    {
                        #define COVERAGE_IGNORE
                        distance_of_new_cell_from_parent = (3.0/2.0)*distance_from_node_to_centroid;
                        #undef COVERAGE_IGNORE
                    }
                    
                    double new_x_value = x + distance_of_new_cell_from_parent*(x_centroid-x);
                    double new_y_value = y + distance_of_new_cell_from_parent*(y_centroid-y);
                    
                    //std::cout << "Parent node at x = " << x << "  y = " << y << "\n";
                    //std::cout << "Daughter node at x = " << new_x_value << "  y = " << new_y_value << "\n";

                    Point<2> new_point(new_x_value, new_y_value);
                    unsigned new_node_index = mrMesh.RefineElement(p_element, new_point);
                    std::cout << "New Cell Index = " << new_node_index << "\n";
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
                    if(mrMesh.GetNumNodes() > mIsGhostNode.size())
                    {
                        #define COVERAGE_IGNORE
                        mIsGhostNode.resize(mrMesh.GetNumNodes());
                        mIsGhostNode[new_node_index] = false;
                        #undef COVERAGE_IGNORE
                    }
                    num_births++;
                    //std::cout<< "num_births=" << num_births <<std::endl<< std::flush;
                    if(mReMesh && periodic_cell)
                    {
                        ReMesh();
                    }
                } // if (ready to divide)
            } // cell iteration loop
        } // if (simulation has cell birth)
        
        return num_births;
    }
    
    /**
     * During a simulation time step, process any cell sloughing or death
     * 
     * At the moment we just slough cells by turning them into ghost nodes
     * 
     * CELL DEATH TO BE ADDED INTO THIS METHOD
     * 
     * @return the number of deaths that occurred.
     */
    unsigned DoCellRemoval()
    {
        
        
        unsigned num_deaths=0;
        
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
        
        if (mpCellKiller)
        {
            mpCellKiller->TestAndLabelCellsForApoptosis();
            mpCellKiller->RemoveDeadCells();
            ReMesh();
        }
        
        return num_deaths;
    }

    /**
     * Find a suitable element for a new cell to be born in, i.e. to be refined with the new node.
     * 
     * @param rpOurNode  the node that is giving birth.  This has to be a reference to a pointer,
     *     since if it is a cell on the periodic boundary that is dividing, we may have to look
     *     at the mirror node to find a suitable element, so we may need to change the pointer
     *     used by calling code.
     * @param cell_index  the index of this cell within the mCells vector
     * @param periodicCell  whether this cell is on the periodic boundary
     * @param periodicIndex  the index of this cell on the periodic boundary; 0 if it isn't there
     * @return  a (pointer to a) suitable element
     */
    Element<2,2>* FindElementForBirth(Node<2>*& rpOurNode, unsigned cellIndex,
                                      const bool periodicCell, const unsigned periodicIndex)
    {
        // Pick a random element to start with
        Element<2,2>* p_element = mrMesh.GetElement(rpOurNode->GetNextContainingElementIndex());
        unsigned element_number = mpRandomNumberGenerator->randMod(rpOurNode->GetNumContainingElements());
        for(unsigned j=0; j<element_number; j++)
        {
            p_element = mrMesh.GetElement(rpOurNode->GetNextContainingElementIndex());
        }
        
        unsigned counter = 0; // how many elements we've checked for suitability
        bool is_ghost_element, is_periodic_element;
        do
        {
            // A ghost element has at least 1 ghost node
            is_ghost_element = (   mIsGhostNode[p_element->GetNodeGlobalIndex(0)] 
                                || mIsGhostNode[p_element->GetNodeGlobalIndex(1)] 
                                || mIsGhostNode[p_element->GetNodeGlobalIndex(2)] );
            // Make sure only one of the nodes is periodic (the cell that's dividing)
            // A periodic element has 2 periodic nodes
            is_periodic_element = (   (mIsPeriodicNode[p_element->GetNodeGlobalIndex(0)] 
                                       && mIsPeriodicNode[p_element->GetNodeGlobalIndex(1)]) 
                                   || (mIsPeriodicNode[p_element->GetNodeGlobalIndex(0)]
                                       && mIsPeriodicNode[p_element->GetNodeGlobalIndex(2)])
                                   || (mIsPeriodicNode[p_element->GetNodeGlobalIndex(1)]
                                       && mIsPeriodicNode[p_element->GetNodeGlobalIndex(2)]));

            if (is_ghost_element || is_periodic_element)
            {
                // This element isn't suitable
                counter++;
                if(counter >= rpOurNode->GetNumContainingElements())
                {
                    if(periodicCell)
                    {
                        // Swap to the image node - that might have a 
                        // non-periodic element to put new cell into.
                        if(cellIndex == mRightCryptBoundary[periodicIndex])
                        {
                            // We already swapped; give up
                            #define COVERAGE_IGNORE
                            assert(0);  
                            #undef COVERAGE_IGNORE
                        }
                        rpOurNode = mrMesh.GetNode(mRightCryptBoundary[periodicIndex]);
                    }
                    else
                    {
                        // somehow every connecting element is a ghost element. quit to
                        // avoid infinite loop
                        #define COVERAGE_IGNORE
                        assert(0);
                        #undef COVERAGE_IGNORE
                    }
                }
                p_element = mrMesh.GetElement(rpOurNode->GetNextContainingElementIndex());
            }
        } while (is_ghost_element || is_periodic_element);

        return p_element;
    }
    
    /**
     * Calculates the forces on each node
     * 
     * @return drdt the x and y force components on each node
     */
    std::vector<std::vector<double> > CalculateForcesOnEachNode()
    {
        std::vector<std::vector<double> > drdt(mrMesh.GetNumAllNodes());
        for (unsigned i=0; i<mrMesh.GetNumAllNodes(); i++)
        {
            drdt[i].resize(2);
        }
        //////////////////////////////////////////////////////////////////
        // loop over element and for each one loop over it's three edges
        ////////////////////////////////////////////////////////////////////
        for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
        {
            Element<2,2>* p_element = mrMesh.GetElement(elem_index);
            if(!p_element->IsDeleted())
            {
                for (unsigned k=0; k<3; k++)
                {
                    unsigned nodeA = k, nodeB = (k+1)%3;
                    
                    assert(!p_element->GetNode(nodeA)->IsDeleted());
                    assert(!p_element->GetNode(nodeB)->IsDeleted());
                    
                    c_vector<double, 2> drdt_contribution = CalculateForceInThisSpring(p_element,nodeA,nodeB);
                              
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
                        if(!mIsGhostNode[p_element->GetNodeGlobalIndex(nodeA)])
                        {
                            drdt[ p_element->GetNode(nodeB)->GetIndex()][0] -= drdt_contribution(0);
                            drdt[ p_element->GetNode(nodeB)->GetIndex()][1] -= drdt_contribution(1);

                            if(!mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)])
                            {
                                drdt[ p_element->GetNode(nodeA)->GetIndex()][0] += drdt_contribution(0);
                                drdt[ p_element->GetNode(nodeA)->GetIndex()][1] += drdt_contribution(1);
                            }
                        }
                        else
                        {
                            drdt[ p_element->GetNode(nodeA)->GetIndex()][0] += drdt_contribution(0);
                            drdt[ p_element->GetNode(nodeA)->GetIndex()][1] += drdt_contribution(1);
 
                            if(mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)])
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
                unsigned nodeA = 0;
                unsigned nodeB = 1;
                
                assert(!p_edge->GetNode(nodeA)->IsDeleted());
                assert(!p_edge->GetNode(nodeB)->IsDeleted());
                
                c_vector<double, 2> drdt_contribution = CalculateForceInThisBoundarySpring(p_edge);
                               
                // Assume that if both nodes are real, or both are ghosts, then they both
                // exert forces on each other, but if one is real and one is ghost then
                // the real node exerts a force on the ghost node, but the ghost node 
                // does NOT exert a force on the real node.   
                if(!mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeA)])
                {
                        // Real A force on any B
                        drdt[ p_edge->GetNode(nodeB)->GetIndex()][0] -= drdt_contribution(0);
                        drdt[ p_edge->GetNode(nodeB)->GetIndex()][1] -= drdt_contribution(1);

                        // B exerts a force back if it is real.
                        if(!mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)])
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
                    if(mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)])
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
        
        // Here we divide all the foces on the nodes by a factor of two because
        // we looped over them all twice to deal with the boundaries above.
        for(unsigned i=0 ; i<mrMesh.GetNumAllNodes(); i++)
        {
            drdt[i][0]=drdt[i][0]/2.0;
            drdt[i][1]=drdt[i][1]/2.0;
        }
        
        return drdt;    
    }
    
    /**
     * @return the x and y forces in this spring
     */
    c_vector<double, 2> CalculateForceInThisSpring(Element<2,2>*& rPElement,const unsigned& rNodeA,const unsigned& rNodeB)
    {
        c_vector<double, 2> drdt_contribution;
        c_vector<double, 2> unit_difference;
        unit_difference(0)=rPElement->GetNodeLocation(rNodeB,0)-rPElement->GetNodeLocation(rNodeA,0);
        unit_difference(1)=rPElement->GetNodeLocation(rNodeB,1)-rPElement->GetNodeLocation(rNodeA,1);
        double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1));
        
        unit_difference=unit_difference/distance_between_nodes;

        double rest_length = 1.0;
        
        if( (mCells.size()>0) &&  (!mIsGhostNode[rPElement->GetNodeGlobalIndex(rNodeA)]) && (!mIsGhostNode[rPElement->GetNodeGlobalIndex(rNodeB)]) )
        {
            double ageA = mCells[rPElement->GetNode(rNodeA)->GetIndex()].GetAge();
            double ageB = mCells[rPElement->GetNode(rNodeB)->GetIndex()].GetAge();
            if (ageA<1.0 && ageB<1.0 && fabs(ageA-ageB)<1e-6)
            {
               // Spring Rest Length Increases to normal rest length from 0.9 to normal rest length, 1.0, over 1 hour
               #define COVERAGE_IGNORE
               rest_length=(0.1+0.9*ageA);
               assert(rest_length<=1.0);
               #undef COVERAGE_IGNORE
           }                          
        }

        return drdt_contribution = mpParams->GetMeinekeLambda() * unit_difference * (distance_between_nodes - rest_length);
    }
    
    /**
     * @return the x and y forces in this boundary spring
     */
    c_vector<double, 2> CalculateForceInThisBoundarySpring(BoundaryElement<1,2>*& rPEdge)
    {
        const unsigned nodeA = 0;
        const unsigned nodeB = 1;
        c_vector<double, 2> drdt_contribution;
        c_vector<double, 2> unit_difference;
        unit_difference(0)=rPEdge->GetNodeLocation(nodeB,0)-rPEdge->GetNodeLocation(nodeA,0);
        unit_difference(1)=rPEdge->GetNodeLocation(nodeB,1)-rPEdge->GetNodeLocation(nodeA,1);
        double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1));
        
        unit_difference=unit_difference/distance_between_nodes;

        double rest_length = 1.0;
        
        if( (mCells.size()>0) &&  (!mIsGhostNode[rPEdge->GetNodeGlobalIndex(nodeA)]) && (!mIsGhostNode[rPEdge->GetNodeGlobalIndex(nodeB)]) )
        {
            double ageA = mCells[rPEdge->GetNode(nodeA)->GetIndex()].GetAge();
            double ageB = mCells[rPEdge->GetNode(nodeB)->GetIndex()].GetAge();
            if (ageA<1.0 && ageB<1.0 && fabs(ageA-ageB)<1e-6)
            {
               // Spring Rest Length Increases to normal rest length from 0.9 to normal rest length, 1.0, over 1 hour
               #define COVERAGE_IGNORE
               rest_length=(0.1+0.9*ageA);
               assert(rest_length<=1.0);
               #undef COVERAGE_IGNORE
           }                          
        }

        return drdt_contribution = mpParams->GetMeinekeLambda() * unit_difference * (distance_between_nodes - rest_length);
    }
    
    
    
    /**
     * Moves each node to a new position for this timestep
     * 
     * @param rDrDt the x and y force components on each node.
     */
    void UpdateNodePositions(const std::vector< std::vector<double> >& rDrDt)
    {
        for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
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
                                    Point<2> new_point = GetNewNodeLocation(index,rDrDt);
                                    mrMesh.SetNode(index, new_point, false);
                                }
                            }
                        }
                    }
                }
            }
            else if(mCells.size()>0)
            {
                // move any node as long as it is not a real stem cell.
                if(mCells[index].GetCellType()!=STEM || mIsGhostNode[index])
                {
                    if(!mrMesh.GetNode(index)->IsDeleted())
                    {
                        Point<2> new_point = GetNewNodeLocation(index,rDrDt);
                        // if a cell wants to move below y<0 (most likely because it was 
                        // just born from a stem cell), stop it doing so
                        if( (new_point.rGetLocation()[1] < 0.0) && (!mIsGhostNode[index]))
                        {
                            // Here we give the cell a push upwards so that it doesn't get stuck on y=0 for ever.
                            // it is a bit of a hack to make it work nicely!
                            new_point.rGetLocation()[1] = 0.01;
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
                        Point<2> new_point = GetNewNodeLocation(index,rDrDt);
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
                unsigned RightNodeIndex = mRightCryptBoundary[i];
                unsigned LeftNodeIndex = mLeftCryptBoundary[i];
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
    }
    
    Point<2> GetNewNodeLocation(const unsigned& rOldNodeIndex, const std::vector< std::vector<double> >& rDrDt)
    {
        c_vector<double,2> old_point = mrMesh.GetNode(rOldNodeIndex)->rGetLocation();
        Point<2> new_point;
                                    
        // Euler style update to node position
        new_point.rGetLocation()[0] = old_point[0] + mDt*rDrDt[rOldNodeIndex][0]; 
        new_point.rGetLocation()[1] = old_point[1] + mDt*rDrDt[rOldNodeIndex][1]; 
                            
        return new_point;   
    }
    
    /**
     * Change the state of cells
     * 
     * At the moment this turns cells to be differentiated 
     * dependent on a protein concentration when using the Wnt model. 
     */
    void UpdateCellTypes()
    {
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
        {   // Cycle through each cell
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
//                  if(betaCateninLevel < 0.4954)
//                          {
//                              cell_type = DIFFERENTIATED;
//                          }

                        
                        mCells[i].SetCellType(cell_type);
                    }
                }
            }
        }
    }
    
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
        mNodesMoved=false;
        mRemeshesThisTimeStep=0;
        mpCellKiller=NULL;
        mNumBirths = 0;
        mNumDeaths = 0;
        mPeriodicDivisionBuffer = 0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        if(!p_simulation_time->IsStartTimeSetUp())
        {
            EXCEPTION("Start time not set in simulation time singleton object");
        }
        
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
    
    /** 
     * Set the timestep of the simulation
     */
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
    
    /**
     * Sets the maximum number of cells that the simulation will contain (for use by the datawriter)
     * default value is set to 10x the initial mesh value by the constructor. 
     */
    void SetMaxCells(unsigned maxCells)
    {
        mMaxCells = maxCells;
        if(maxCells<mrMesh.GetNumAllNodes())
        {
        	EXCEPTION("mMaxCells is less than the number of cells in the mesh.");	
        }
    }
    
    /**
     * Sets the maximum number of elements that the simulation will contain (for use by the datawriter)
     * default value is set to 10x the initial mesh value by the constructor.     
     */
    void SetMaxElements(unsigned maxElements)
    {
        mMaxElements = maxElements;
        if(maxElements<mrMesh.GetNumAllElements())
        {
        	EXCEPTION("mMaxElements is less than the number of elements in the mesh.");	
        }
    }
    
    /**
     * Call this before Solve() to fix the boundary of the mesh.
     * \todo figure out what this does!
     */
    void SetFixedBoundaries()
    {
        mFixedBoundaries = true;
    }
    
    /**
     *  Call this before Solve() to set the boundary conditions
     * i.e. whether the left and right boundaries should be periodic
     */
    void SetPeriodicSides(bool periodicSides)
    {
        mPeriodicSides = periodicSides;
    }
    
    void SetCellKiller(RandomCellKiller<2>* pCellKiller)
    {
        mpCellKiller=pCellKiller;
        mpCellKiller->SetCellsAndMesh(&mCells, &mrMesh);
    }
    
    /** 
     * Get the cells vector
     * N.B. Returns a copy of the cells - any operations on them will not go back into the
     * simulation.
     * \todo change this to return a const reference
     */
    std::vector<MeinekeCryptCell> GetCells()
    {
        assert(mCells.size()>0);
        return mCells;
    }
    
    /** 
     * Whether each node is a ghost or not.
     * \todo change this to return a const reference
     */
    std::vector <bool> GetGhostNodes()
    {
        return mIsGhostNode;
    }
    
    /** 
     * Return the index of each node on the left periodic boundary
     * \todo change this to return a const reference
     */
    std::vector<unsigned> GetLeftCryptBoundary()
    {
        return mLeftCryptBoundary;
    }
    
    /** 
     * Return the index of each node on the right periodic boundary
     * \todo change this to return a const reference
     */
    std::vector<unsigned> GetRightCryptBoundary()
    {
        return mRightCryptBoundary;
    }
    
    /** 
     * Return the index of each node on the whole boundary
     * \todo change this to return a const reference
     */    
    std::vector<unsigned> GetCryptBoundary()
    {
        return mCryptBoundary;
    }
    
    /** 
     * This automatically sets this to be a wnt dependent simulation.
     * You should supply cells with a wnt cell cycle...
     * 
     */
    void SetWntGradient(WntGradientType wntGradientType)
    {
    	mWntIncluded = true;
    	mWntGradient = WntGradient(wntGradientType);
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
        
        ///////////////////////////////////////////////////////////
        // Set up Simulation
        ///////////////////////////////////////////////////////////
        
        // Data writers for tabulated results data, used in tests
        ColumnDataWriter tabulated_node_writer(mOutputDirectory+"Results", "tabulated_node_results");
        ColumnDataWriter tabulated_element_writer(mOutputDirectory+"Results", "tabulated_element_results");
        
        node_writer_ids_t node_writer_ids;
        SetupNodeWriter(tabulated_node_writer, node_writer_ids);
        
        element_writer_ids_t element_writer_ids;
        SetupElementWriter(tabulated_element_writer, element_writer_ids);
        
        // This keeps track of when tabulated results were last output
        unsigned tabulated_output_counter = 0;
        
        // Create output files for the visualizer
        OutputFileHandler output_file_handler(mOutputDirectory);
        out_stream p_node_file = output_file_handler.OpenOutputFile("results.viznodes");
        out_stream p_element_file = output_file_handler.OpenOutputFile("results.vizelements");
        
        
        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double current_time = p_simulation_time->GetDimensionalisedTime();
        std::cout << "Time at start of Solve Method = " << current_time << std::endl;
        
        unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);
        std::cout << "num timesteps = " << num_time_steps << std::endl;
        if(current_time>0)//use the reset function if necessary
        {
        	p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);    
        }
        else
        {
       		p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);    
        }
        
        
        // Check some parameters for a periodic simulation
        if(mPeriodicSides)
        {
        	CalculateCryptBoundary();
        	if(mLeftCryptBoundary.size()<1 || mRightCryptBoundary.size()<1)
        	{
        		EXCEPTION("Periodic Simulation but mesh is not periodic\nIf you want a non-periodic simulation use SetPeriodicSides(false)");
        	}
        	if(!mReMesh)
        	{
        		#define COVERAGE_IGNORE
        		EXCEPTION("A periodic simulation requires active remeshing\n");
        		#undef COVERAGE_IGNORE
        	}
        }
        
        /* Age the cells to the correct time (cells set up with negative birth dates 
         * to give some that are almost ready to divide).
         * 
         * TODO:For some strange reason this seems to take about 3 minutes for a realistic Wnt-Crypt.
         * Not sure why - when the same code was evaluated in a test it seemed almost instant.
         */
        if (!mCells.empty())
	    {
            bool temp;
	    	for(unsigned i=0; i<mCells.size(); i++)
        	{
		    	if(mIsGhostNode[i]) continue;
		    	//std::cout << "Preparing Cell "<< i << std::endl;
		    	Node<2> *p_our_node = mrMesh.GetNode(i);
		        double y = p_our_node->rGetLocation()[1];
		        std::vector<double> cell_cycle_influences;
		        if(mWntIncluded)
		        {
			    	double wnt_stimulus = mWntGradient.GetWntLevel(y);
			    	cell_cycle_influences.push_back(wnt_stimulus);
		        }
                // We don't use the result; this call is just to force the cells to age to time 0,
                // running their cell cycle models to get there.
                temp = mCells[i].ReadyToDivide(cell_cycle_influences);
            }
	    }
	    
        
        /////////////////////////////////////////////////////////////////////
        // Main time loop
        /////////////////////////////////////////////////////////////////////
        
        while (p_simulation_time->GetTimeStepsElapsed() < num_time_steps)
        {
        	mRemeshesThisTimeStep = 0; // To avoid infinite loops
		    std::cout << "** TIME = " << p_simulation_time->GetDimensionalisedTime() << " **" << std::endl;
		                
            // Cell birth
            mNumBirths += DoCellBirth();
            
            //  calculate node velocities
            std::vector<std::vector<double> > drdt = CalculateForcesOnEachNode();
            
            // update node positions
            UpdateNodePositions(drdt);
            
            //////////////////////////////////////////////
            // Cell death should be included in this method
            /////////////////////////////////////////////
            mNumDeaths += DoCellRemoval();
            
            
            // Change the state of some cells
            // Only active for WntCellCycleModel at the moment
            // but mutations etc. could occur in this function
            UpdateCellTypes();


    		ReMesh();
    		
    		// Increment simulation time here, so results files look sensible
            p_simulation_time->IncrementTimeOneStep();            
            
            // Write results to file
            WriteResultsToFiles(tabulated_node_writer, node_writer_ids,
                                tabulated_element_writer, element_writer_ids,
                                *p_node_file, *p_element_file,
                                tabulated_output_counter==0,
                                true);
			
            tabulated_output_counter++; 
            if(tabulated_output_counter > 80) // TODO: make this configurable!
            {
                tabulated_output_counter = 0;
            }
        } // End main time loop

        
        tabulated_node_writer.Close();
        tabulated_element_writer.Close();
    }
    
    
    /**
     * The mesh should be surrounded by at least one layer of ghost nodes.  These are nodes which 
     * do not correspond to a cell, but are necessary for remeshing (because the remesher tries to 
     * create a convex hull of the set of nodes) and visualising purposes.  The mesh is passed into
     * the constructor and the class is told about the ghost nodes by using this method. 
     */     
    void SetGhostNodes(std::vector<unsigned> ghostNodeIndices)
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
        for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
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
		
//		// Work out if there has been a change in the boundary since last time this was run
//		if(mLeftCryptBoundary!=mOldLeftCryptBoundary)
//		{
//			std::cout << "Left Crypt Boundary vector changed.\n";
//			
//		}
//		if (mRightCryptBoundary!=mOldRightCryptBoundary)
//		{
//			std::cout << "Right Crypt Boundary vector changed.\n";	
//		}
//		if(mCryptBoundary!=mOldCryptBoundary)
//		{
//			std::cout << "Total Boundary vector changed.\n";	
//		}
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
		
		/* For each node see if it has broken onto a periodic boundary
    	 * If it has create an image node at the other side
    	 */
    	unsigned i = 0;
    	while(i<mCryptBoundary.size())
        {
        	std::vector<unsigned> nodes_on_left_boundary  =  mLeftCryptBoundary;
        	std::vector<unsigned> nodes_on_right_boundary =  mRightCryptBoundary;
        	std::vector<unsigned> nodes_on_boundary 	  =  mCryptBoundary;
	        unsigned our_node = nodes_on_boundary[i];
	        //std::cout << our_node << "\n";
	        bool our_node_periodic = false;
	        bool left_break=false;
	        bool right_break = false;
	        
	        std::vector<unsigned> periodic;
	        
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
	        	// Cycle through each pair of elements attached to this node
				// are they attached to a shared ghost node?
				// are they each attached to a different periodic node?
				// If so they are the periodic nodes of interest.
				bool culprit_found = false;
				
				//\TODO replace this section with a method that just scans the elements attached to the node of interest.
	        	for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
		        {
		        	Element<2,2>* p_element = mrMesh.GetElement(elem_index);
		    		unsigned element1_node[3];
		    		element1_node[0] = (unsigned)p_element->GetNode(0)->GetIndex();
		    		element1_node[1] = (unsigned)p_element->GetNode(1)->GetIndex();
		    		element1_node[2] = (unsigned)p_element->GetNode(2)->GetIndex();
		    		
		    		if(element1_node[0]==our_node || element1_node[1]==our_node || element1_node[2]==our_node)
		    		{
		    			for (unsigned element2_index = 0; element2_index<mrMesh.GetNumAllElements() ; element2_index++)	
		    			{	
		    				if (elem_index!=element2_index)
			    			{
			    				Element<2,2>* p_element2 = mrMesh.GetElement(element2_index);
					    		unsigned element2_node[3];
						    	element2_node[0] = (unsigned)p_element2->GetNode(0)->GetIndex();
						    	element2_node[1] = (unsigned)p_element2->GetNode(1)->GetIndex();
			    				element2_node[2] = (unsigned)p_element2->GetNode(2)->GetIndex();
			    				

			    				
			    				if(element2_node[0]==our_node || element2_node[1]==our_node || element2_node[2]==our_node)
			    				{
			    					unsigned ghost_node_element1 = 0;//it will never be this because this will be in bottom left corner of ghosts
			    					unsigned ghost_node_element2 = 0;//it will never be this because this will be in bottom left corner of ghosts
			    					
		    						for(unsigned index = 0 ; index <3 ; index++)
		    						{
		    							if(mIsGhostNode[element1_node[index]])
		    							{
		    								ghost_node_element1 = element1_node[index];
		    							}
		    						}
		    						for(unsigned index = 0 ; index <3 ; index++)
		    						{
		    							if(mIsGhostNode[element2_node[index]])
		    							{
		    								ghost_node_element2 = element2_node[index];
		    							}
		    						}
//		    						std::cout << "Examining two elements attached to node " << our_node << ".\n";
		    						// this mucked up periodic births - new cell attached to more than one ghost node.
		    						//if(ghost_node_element1==ghost_node_element2 && ghost_node_element1>0)
		    						// if there is a ghost node in each of the two elements under consideration
		    	 					if(ghost_node_element1>0 && ghost_node_element2>0)
		    						{
		    							//std::cout << "node "<< our_node << " is attached to ghost nodes " << ghost_node_element1 << " and " << ghost_node_element2 << "\n";
		    							//Now check they are both attached to different periodic nodes
		    							unsigned other_node_element1 = 0;
		    							unsigned other_node_element2 = 0;
		    							
		    							bool left_nodes=false;
		    							for(unsigned j = 0 ; j<nodes_on_left_boundary.size() ; j++)
										{	//Cycle through periodic nodes and check for periodic neighbours
						    				unsigned periodic_node = nodes_on_left_boundary[j];
											if(element1_node[0]==periodic_node || element1_node[1]==periodic_node || element1_node[2]==periodic_node)
											{
												other_node_element1 = periodic_node;
												left_nodes=true;
											}
											if(element2_node[0]==periodic_node || element2_node[1]==periodic_node || element2_node[2]==periodic_node)
											{
												other_node_element2 = periodic_node;
												left_nodes=true;
											}
											periodic_node = nodes_on_right_boundary[j];
											if(element1_node[0]==periodic_node || element1_node[1]==periodic_node || element1_node[2]==periodic_node)
											{
												other_node_element1 = periodic_node;
											}
											if(element2_node[0]==periodic_node || element2_node[1]==periodic_node || element2_node[2]==periodic_node)
											{
												other_node_element2 = periodic_node;
											}
											if(other_node_element1!=other_node_element2 && other_node_element1>0 && other_node_element2>0 && !culprit_found)
											{
//												std::cout << "\nElement " << elem_index << " contains nodes " << element1_node[0] << ", " << element1_node[1] << ", " << element1_node[2] << ".\n";
//												std::cout << "Element " << element2_index << " contains nodes " << element2_node[0] << ", " << element2_node[1] << ", " << element2_node[2] << ".\n";
//												std::cout << "We are considering node " << our_node << "\n";
//		    									std::cout << "attached to ghost node " << ghost_node_element1 << "\n";
//		    									std::cout << "attached to ghost node " << ghost_node_element2 << "\n";
//		    									std::cout << "and periodic node " << other_node_element1 << " by element "<< elem_index <<"\n";
//		    									std::cout << "and periodic node " << other_node_element2 << " by element "<< element2_index <<"\n";
//		    									//assert(0);
		    									periodic.push_back(other_node_element1);
		    									periodic.push_back(other_node_element2);
		    									if(left_nodes)
		    									{
		    										left_break = true;	
		    									}
		    									else
		    									{
		    										right_break = true;	
		    									}
		    									culprit_found=true;
											}
						    			}// end of loop through periodic nodes
		    						}// end of if our node and shared ghosts are in these two elements
								}// end of if our node is in both elements
		    				}// end of if elements are different
		    			}// end of element 2 loop
		    		} // end of if our node is in this element
				}// end of element 1 loop
					
				if(left_break)
		        {
                    // We should have a new periodic node
		        	double old_x = mrMesh.GetNode(our_node)->rGetLocation()[0];
		        	double old_y = mrMesh.GetNode(our_node)->rGetLocation()[1];
		        	double crypt_width = mpParams->GetCryptWidth();
//		        	std::cout << "LEFT Node "<< our_node << " has broken into the periodic edge between nodes\n";
//		        	for(unsigned k=0 ; k<periodic.size() ; k++)
//					{
//						std::cout << periodic[k] << "\t";
//					}
//					std::cout << "\n";
					AddACellToPeriodicBoundary(our_node,old_x+crypt_width,old_y,periodic);
					
					ReMesh();
				}
		        
		        if(right_break)
		        {
		        	// We should have a new periodic node
		        	double old_x = mrMesh.GetNode(our_node)->rGetLocation()[0];
		        	double old_y = mrMesh.GetNode(our_node)->rGetLocation()[1];
		        	double crypt_width = mpParams->GetCryptWidth();
//		        	std::cout << "RIGHT Node "<< our_node << " has broken into the periodic edge between nodes\n";
//					for(unsigned k=0 ; k<periodic.size() ; k++)
//					{
//						std::cout << periodic[k] << "\t";
//					}
//					std::cout << "\n";
					AddACellToPeriodicBoundary(our_node,old_x-crypt_width,old_y,periodic);
					
					ReMesh();
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
			
			// Now this_<side>_node_missing is only true if the node has been internalised (or is a ghost already)
			if(this_left_node_missing && (!mIsGhostNode[old_node_on_left_boundary]))
			{	// The left node has been internalised (was periodic and is not a ghost) so the right node should be spooked
				mIsGhostNode[old_node_on_right_boundary]=true;
				std::cout << "Right Node " << old_node_on_right_boundary << " spooked\n";
				//mNodesMoved=true;
				CalculateCryptBoundary();
			}
			if(this_right_node_missing && (!mIsGhostNode[old_node_on_right_boundary]))
			{	// The right node has been internalised (was periodic and is not a ghost) so the right node should be spooked
				mIsGhostNode[old_node_on_left_boundary]=true;
				std::cout << "Left Node " << old_node_on_left_boundary << " spooked\n";
		    	//mNodesMoved=true;
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
		mNodesMoved=true;
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
                	periodic_nodes[j] = mRightCryptBoundary[i];
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
        for (unsigned elim_index = 0 ; elim_index<mrMesh.GetNumAllElements(); elim_index++ )
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
    
    
    /**
     * This method should be called by all other methods that require a re-mesh
     * it looks for changes in the periodic boundaries that could require
     * more remeshes and carries out the correct number.     * 
     */
    void ReMesh()
    {
    	if(mReMesh)
    	{
	    	CallReMesher();
	
			if(mPeriodicSides)
			{
				mOldLeftCryptBoundary = mLeftCryptBoundary;
				mOldRightCryptBoundary = mRightCryptBoundary;
				mOldCryptBoundary = mCryptBoundary;
			    
			    CalculateCryptBoundary();
			    
			    // do this once every time (sometimes mCryptBoundary does not change
			    // but the connections between the nodes have changed)
			    RemoveSurplusCellsFromPeriodicBoundary();
				DetectNaughtyCellsAtPeriodicEdges();  
				
				while(mNodesMoved)
				{
					CallReMesher();
				
					mOldLeftCryptBoundary = mLeftCryptBoundary;
					mOldRightCryptBoundary = mRightCryptBoundary;
					mOldCryptBoundary = mCryptBoundary;
				
					CalculateCryptBoundary();
			    	
		        	RemoveSurplusCellsFromPeriodicBoundary();
		        	DetectNaughtyCellsAtPeriodicEdges();
				}
			}
    	}
	}
	
	/**
	 * This method actually calls the remesh command on the mesh.
	 * It should only be called by the method above (ReMesh) which ensures
	 * that periodic boundaries are handled properly.
	 */
	void CallReMesher()
	{
		std::cout << "Remeshing \n"<< std::flush;
		NodeMap map(mrMesh.GetNumAllNodes());
    	mrMesh.ReMesh(map);
        
        // TODO: These commented out because they caused a segmentation 
        // fault after the Load function has been called.
        // Possibly necessary for cell death - but missing a method to actually
        // make the cells vector smaller.
//        for (unsigned i=0; i<mCells.size(); i++)
//        {
//            
//            unsigned old_index = mCells[i].GetNodeIndex();
//            unsigned new_index = map.GetNewIndex(old_index);
//            mCells[i].SetNodeIndex(new_index);
//        }
        
    	mNodesMoved=false;
    	mRemeshesThisTimeStep++;
    	assert(mRemeshesThisTimeStep < 1000); //to avoid an infinite loop. If this ever throws try increasing it a bit.
    }
    
    /**
     * Saves the whole crypt simulation for restarting later.
     * 
     * Puts it in the folder mOutputDirectory/archive/
     * 
     * First archives simulation time then the simulation itself.
     */
    void Save()
    {
        // todo: remesh, save mesh

        std::string archive_directory = mOutputDirectory + "/archive/";
        
        // create an output file handler in order to get the full path of the 
        // archive directory. Note the false is so the handler doesn't clean
        // the directory
        OutputFileHandler handler(archive_directory, false);
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "crypt_sim_periodic_2d.arch";
        
        std::ofstream ofs(archive_filename.c_str());       
        boost::archive::text_oarchive output_arch(ofs);

        SimulationTime* p_sim_time = SimulationTime::Instance();
        assert(p_sim_time->IsStartTimeSetUp());

        // cast to const.
        const SimulationTime* p_simulation_time = SimulationTime::Instance();
        output_arch << *p_simulation_time;
        output_arch << static_cast<const CryptSimulation2DPeriodic&>(*this);        
    }
    
    /**
     * Loads a saved crypt simulation
     * 
     * @param rArchiveDirectory the name of the simulation to load 
     * (specified originally by simulator.SetOutputDirectory("wherever"); )
     */
    void Load(const std::string& rArchiveDirectory)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        
        OutputFileHandler any_old_handler("",false);
		std::string test_output_directory = any_old_handler.GetTestOutputDirectory();
        
        std::string archive_filename = test_output_directory + rArchiveDirectory + "/archive/crypt_sim_periodic_2d.arch";
        
        // Create an input archive
        std::ifstream ifs(archive_filename.c_str(), std::ios::binary);       
        boost::archive::text_iarchive input_arch(ifs);

        // read the archive
        assert(p_simulation_time->IsStartTimeSetUp());
        input_arch >> *p_simulation_time;        
        input_arch >> *this;
		
		
        mOutputDirectory = rArchiveDirectory+"/load_temp";

        if(mrMesh.GetNumNodes()!=mCells.size())
        {
            EXCEPTION("Number of Nodes is not equal to number of cells. This is very bad.");   
        }
    }
};

#endif /*CRYPTSIMULATION2DPERIODIC_HPP_*/
