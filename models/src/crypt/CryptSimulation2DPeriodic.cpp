#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp> // for archiving vectors
#include <boost/serialization/string.hpp>

#include <fstream>

#include "Cylindrical2dMesh.cpp"
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
#include "CryptSimulation2DPeriodic.hpp"
#include "TrianglesMeshReader.hpp"

/** Constructor
 *  @param cells is defaulted to the empty vector, in which case SetIncludeRandomBirth()
 *  should be called for any birth to happen.
 */
CryptSimulation2DPeriodic::CryptSimulation2DPeriodic(ConformingTetrahedralMesh<2,2> &rMesh,
                                                     std::vector<MeinekeCryptCell> cells)
        : mrMesh(rMesh),
          mCells(cells),
          mCrypt(rMesh, mCells)
          
          
{ 
    
    mpParams = CancerParameters::Instance();
    
    mDt = 1.0/(120.0);
    mEndTime = 120.0; //hours
    
    srandom(0);
    
    mIncludeRandomBirth = false;
    mIncludeVariableRestLength = false;
    mFixedBoundaries = false;
    mOutputDirectory = "";
    
    // Set up the ghost nodes bool.  Assume initially that the maximum number of nodes is
    // ten times the mesh size.  Note that more memory is allocated later, if necessary.
    mIsGhostNode.resize(10*mrMesh.GetNumAllNodes()); // Note the hard-coding of 10.
    
    for (unsigned i=0; i<mIsGhostNode.size(); i++)
    {
        mIsGhostNode[i] = false;
    
    }
    
    // defaults
    mReMesh = true;
    mNoBirth = false;
    mMaxCells = 10*mrMesh.GetNumNodes();
    mMaxElements = 10*mrMesh.GetNumElements();
    mWntIncluded = false;
    mpCellKiller = NULL;
    mNumBirths = 0;
    mNumDeaths = 0;
    
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    if (!p_simulation_time->IsStartTimeSetUp())
    {
        EXCEPTION("Start time not set in simulation time singleton object");
    }
    
}

/**
 * Free any memory allocated by the constructor
 */
CryptSimulation2DPeriodic::~CryptSimulation2DPeriodic()
{
    SimulationTime::Destroy();
}

/**
* Define the variable identifiers in the data writer used to write node-based results.
*
* Uses mMaxCells to decide how many variables to define.
*/
void CryptSimulation2DPeriodic::SetupNodeWriter(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rVarIds)
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
void CryptSimulation2DPeriodic::SetupElementWriter(ColumnDataWriter& rElementWriter, element_writer_ids_t& rVarIds)
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

void CryptSimulation2DPeriodic::WriteVisualizerSetupFile(std::ofstream& rSetupFile)
{
    rSetupFile << "MeshWidth\t" << mrMesh.GetWidth(1u);// get furthest distance between nodes in the x-direciton
    rSetupFile.close();
}

void CryptSimulation2DPeriodic::WriteResultsToFiles(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rNodeVarIds,
                                                    ColumnDataWriter& rElementWriter, element_writer_ids_t& rElementVarIds,
                                                    std::ofstream& rNodeFile, std::ofstream& rElementFile,
                                                    bool writeTabulatedResults,
                                                    bool writeVisualizerResults)
{
    // Write current simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetDimensionalisedTime();
    
    if (writeVisualizerResults)
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
        if (index>mMaxCells)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("\nNumber of cells exceeds mMaxCells. Use SetMaxCells(unsigned) to increase it.\n");
            #undef COVERAGE_IGNORE
        }
        unsigned colour = 0; // all green if no cells have been passed in
        
        if (mIsGhostNode[index]==true)
        {
            colour = 4; // visualizer treats '4' these as invisible
        }
        else if (mCells.size()>0)
        {
            if (index < mCells.size())
            {
                CryptCellType type = mCells[index].GetCellType();
                CryptCellMutationState mutation = mCells[index].GetMutationState();
                
                if (type == STEM)
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
                
                if (mutation!=HEALTHY)
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
        
        if (!mrMesh.GetNode(index)->IsDeleted())
        {
            const c_vector<double,2>& r_node_loc = mrMesh.GetNode(index)->rGetLocation();
            if (writeVisualizerResults)
            {
                rNodeFile << r_node_loc[0] << " "<< r_node_loc[1] << " " << colour << " ";
            }
            if (writeTabulatedResults)
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
            if (writeVisualizerResults)
            {
                rElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0)<< " " << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1)<< " "<< mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2)<< " ";
            }
            if (writeTabulatedResults)
            {
                rElementWriter.PutVariable(rElementVarIds.nodeAs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0));
                rElementWriter.PutVariable(rElementVarIds.nodeBs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1));
                rElementWriter.PutVariable(rElementVarIds.nodeCs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2));
            }
        }
    }
    
    if (writeVisualizerResults)
    {
        rNodeFile << "\n";
        rElementFile << "\n";
    }
    if (writeTabulatedResults)
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
 * @return the number of births that occurred.
 */
unsigned CryptSimulation2DPeriodic::DoCellBirth()
{
    unsigned num_births = 0;
    if (!mNoBirth && !mCells.empty())
    {
        mCrypt.InitialiseCellIterator();
        
        unsigned cell_index=0;
        while(!mCrypt.CellIteratorIsAtEnd())
//        // Iterate over all cells, seeing if each one can be divided

        {
            
            assert(mCrypt.GetCurrentCell()->GetNodeIndex()==cell_index);
            

            unsigned node_index=cell_index ;
            bool skip = false; // Whether to skip this cell
            if (mrMesh.GetNode(node_index)->IsDeleted()) skip=true; // Skip deleted cells
            //if (mrMesh.GetNode(cell_index)->IsDead()) skip=true; // Skip dead cells
            if (mIsGhostNode[cell_index]) skip=true; // Skip Ghost nodes
            
            if (!skip)
            {
                // Check for this cell dividing
                // Construct any influences for the cell cycle...
                Node<2> *p_our_node = mrMesh.GetNode(cell_index);
                std::vector<double> cell_cycle_influences;
                if (mWntIncluded)
                {
                    double y = p_our_node->rGetLocation()[1];
                    double wnt_stimulus = mWntGradient.GetWntLevel(y);
                    cell_cycle_influences.push_back(wnt_stimulus);
                }
                
                // CHECK if this cell is ready to divide - if so create a new cell etc.
                if (mCells[cell_index].ReadyToDivide(cell_cycle_influences))
                {
                    // Create new cell
                    MeinekeCryptCell new_cell = mCells[cell_index].Divide();
                    std::cout << "Cell division at node " << cell_index << "\n";
                
                    // Add a new node to the mesh
                    unsigned parent_node_index = mCells[cell_index].GetNodeIndex();
                    c_vector<double, 2> new_location = CalculateDividingCellCentreLocations(parent_node_index);
                    Node<2>* p_new_node = new Node<2>(mrMesh.GetNumNodes(), new_location, false);   // never on boundary
                    
                    NodeMap map(mrMesh.GetNumNodes());
                    unsigned new_node_index = mrMesh.AddNodeAndReMesh(p_new_node,map);
                    // Go through all the cells and update their node indices according to the map
//                    for (unsigned i=0 ; i<mCells.size(); i++)
//                    {
//                        unsigned old_index = mCells[i].GetNodeIndex();
//                        mCells[i].SetNodeIndex(map.GetNewIndex(old_index));
//                    }
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
                    // Update size of IsGhostNode if necessary
                    if (mrMesh.GetNumNodes() > mIsGhostNode.size())
                    {
                        #define COVERAGE_IGNORE
                        mIsGhostNode.resize(mrMesh.GetNumNodes());
                        mIsGhostNode[new_node_index] = false;
                        #undef COVERAGE_IGNORE
                    }
                    num_births++;
                } // if (ready to divide)
            }
            cell_index++;
            mCrypt.IncrementCellIterator();
        } // cell iteration loop
        assert(cell_index==mCells.size());
    } // if (simulation has cell birth)
    
    return num_births;
}

/**
 * Calculates the new locations of a dividing cell's cell centres.
 * Moves the dividing node a bit and returns co-ordinates for the new node.
 * It does this by picking a random direction (0->2PI) and placing the parent 
 * and daughter in opposing directions on this axis.
 * 
 * @param node_index The parent node index
 * 
 * @return daughter_coords The coordinates for the daughter cell.
 * 
 */
c_vector<double, 2> CryptSimulation2DPeriodic::CalculateDividingCellCentreLocations(unsigned node_index)
{
    double separation = 0.1;
    c_vector<double, 2> parent_coords = mrMesh.GetNode(node_index)->rGetLocation();
    c_vector<double, 2> daughter_coords;
    
    // Make a random direction vector of the required length
    double random_direction = RandomNumberGenerator::Instance()->ranf();
    random_direction = random_direction*2.0*M_PI;
    c_vector<double, 2> random_vector;
    random_vector[0] = 0.5*separation*cos(random_direction);
    random_vector[1] = 0.5*separation*sin(random_direction);
    
    // add it to the daughter and take it from the parent location
    
    daughter_coords[0] = parent_coords[0]+random_vector[0];
    daughter_coords[1] = parent_coords[1]+random_vector[1];
    parent_coords[0] = parent_coords[0]-random_vector[0];
    parent_coords[1] = parent_coords[1]-random_vector[1];
    
    // set the parent to use this location
    mrMesh.SetNode(node_index, parent_coords, false);
    return daughter_coords;
}

/**
 *  Checks that the indices are in sync in the cells vector, ie that
 *  mCells[i].GetNodeIndex() is equal to i
 */
void CryptSimulation2DPeriodic::CheckIndicesAreInSync()
{
    if (!mCells.empty())
    {
        for (unsigned cell_index=0; cell_index<mCells.size(); cell_index++)
        {
            assert(mCells[cell_index].GetNodeIndex()==cell_index);
        }
    }
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
unsigned CryptSimulation2DPeriodic::DoCellRemoval()
{
    unsigned num_deaths=0;
    
    ///////////////////////////////////////////////////////////////////////////////////
    // Alternate method of sloughing.  Turns boundary nodes into ghost nodes.
    ///////////////////////////////////////////////////////////////////////////////////
    for (unsigned i=0; i<mrMesh.GetNumNodes(); i++)
    {
        Node<2> *p_node = mrMesh.GetNode(i);
        if (!p_node->IsDeleted())
        {
            double x = p_node->rGetLocation()[0];
            double y = p_node->rGetLocation()[1];
            double crypt_length=mpParams->GetCryptLength();
            double crypt_width=mpParams->GetCryptWidth();
            unsigned node_index = p_node->GetIndex();
            
            if ( (x>crypt_width) || (x<0.0) || (y>crypt_length))
            {
                mIsGhostNode[node_index] = true;
                num_deaths++;
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
 * Calculates the forces on each node
 *
 * @return drdt the x and y force components on each node
 */
std::vector<c_vector<double, 2> > CryptSimulation2DPeriodic::CalculateVelocitiesOfEachNode()
{
    std::vector<c_vector<double, 2> > drdt(mrMesh.GetNumAllNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i]=zero_vector<double>(2);
    }

    ////////////////////////////////////////////////////////////////////
    // loop over element and for each one loop over its three edges
    ////////////////////////////////////////////////////////////////////
    for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
    {
        Element<2,2>* p_element = mrMesh.GetElement(elem_index);
        if (!p_element->IsDeleted())
        {
            for (unsigned k=0; k<3; k++)
            {
                unsigned nodeA = k, nodeB = (k+1)%3;
                assert(!p_element->GetNode(nodeA)->IsDeleted());
                assert(!p_element->GetNode(nodeB)->IsDeleted());
                
                c_vector<double, 2> force = CalculateForceInThisSpring(p_element,nodeA,nodeB);
                 
                double damping_constantA = mpParams->GetDampingConstantNormal();
                double damping_constantB = mpParams->GetDampingConstantNormal();
                
                if(!mCells.empty())
                {
                    //note: at the moment the index into the mCells vector is the same
                    //as the node index. later this may not be the case, in which case
                    //the following assertion will trip. to deal with this, a map from 
                    //node index to cell will be needed
                    assert( mCells[p_element->GetNodeGlobalIndex(nodeA)].GetNodeIndex()==p_element->GetNodeGlobalIndex(nodeA));
                    assert( mCells[p_element->GetNodeGlobalIndex(nodeB)].GetNodeIndex()==p_element->GetNodeGlobalIndex(nodeB));
                    
                    if(   (mCells[p_element->GetNodeGlobalIndex(nodeA)].GetMutationState()==HEALTHY)
                       || (mCells[p_element->GetNodeGlobalIndex(nodeA)].GetMutationState()==APC_ONE_HIT))
                    {
                        damping_constantA = mpParams->GetDampingConstantNormal();
                    }
                    else
                    {
                        damping_constantA = mpParams->GetDampingConstantMutant();
                    }
                    
                    if(   (mCells[p_element->GetNodeGlobalIndex(nodeB)].GetMutationState()==HEALTHY)
                       || (mCells[p_element->GetNodeGlobalIndex(nodeB)].GetMutationState()==APC_ONE_HIT))
                    {
                        damping_constantB = mpParams->GetDampingConstantNormal();
                    }
                    else
                    {
                        damping_constantB = mpParams->GetDampingConstantMutant();
                    }
                }
                
                if (!mIsGhostNode[p_element->GetNodeGlobalIndex(nodeA)])
                {
                    drdt[ p_element->GetNode(nodeB)->GetIndex()] -= force / damping_constantB;
                    
                    if (!mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)])
                    {
                        drdt[ p_element->GetNode(nodeA)->GetIndex()] += force / damping_constantA;
                    }
                }
                else
                {
                drdt[ p_element->GetNode(nodeA)->GetIndex()] += force / damping_constantA;
                        
                    if (mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)])
                    {
                        drdt[ p_element->GetNode(nodeB)->GetIndex()] -= force / damping_constantB;
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
        if (!p_edge->IsDeleted())
        {
            unsigned nodeA = 0;
            unsigned nodeB = 1;
            
            assert(!p_edge->GetNode(nodeA)->IsDeleted());
            assert(!p_edge->GetNode(nodeB)->IsDeleted());
            
            c_vector<double, 2> force = CalculateForceInThisBoundarySpring(p_edge);
              
            double damping_constantA = mpParams->GetDampingConstantNormal();
            double damping_constantB = mpParams->GetDampingConstantNormal();
            
            if(!mCells.empty())
            {
                //note: at the moment the index into the mCells vector is the same
                //as the node index. later this may not be the case, in which case
                //the following assertion will trip. to deal with this, a map from 
                //node index to cell will be needed
                assert( mCells[p_edge->GetNodeGlobalIndex(nodeA)].GetNodeIndex()==p_edge->GetNodeGlobalIndex(nodeA));
                assert( mCells[p_edge->GetNodeGlobalIndex(nodeB)].GetNodeIndex()==p_edge->GetNodeGlobalIndex(nodeB));
                
                if(   (mCells[p_edge->GetNodeGlobalIndex(nodeA)].GetMutationState()==HEALTHY)
                   || (mCells[p_edge->GetNodeGlobalIndex(nodeA)].GetMutationState()==APC_ONE_HIT))
                {
                    damping_constantA = mpParams->GetDampingConstantNormal();
                }
                else
                {
                    damping_constantA = mpParams->GetDampingConstantMutant();
                }
                
                if(   (mCells[p_edge->GetNodeGlobalIndex(nodeB)].GetMutationState()==HEALTHY)
                   || (mCells[p_edge->GetNodeGlobalIndex(nodeB)].GetMutationState()==APC_ONE_HIT))
                {
                    damping_constantB = mpParams->GetDampingConstantNormal();
                }
                else
                {
                    damping_constantB = mpParams->GetDampingConstantMutant();
                }
            }
                        
            // Assume that if both nodes are real, or both are ghosts, then they both
            // exert forces on each other, but if one is real and one is ghost then
            // the real node exerts a force on the ghost node, but the ghost node
            // does NOT exert a force on the real node.
            if (!mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeA)])
            {
                // Real A force on any B
                drdt[ p_edge->GetNode(nodeB)->GetIndex()] -= force / damping_constantB;
                
                // B exerts a force back if it is real.
                if (!mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)])
                {
                    drdt[ p_edge->GetNode(nodeA)->GetIndex()] += force / damping_constantA;
                }
            }
            else
            {
                // Ghost A receives a force
                drdt[ p_edge->GetNode(nodeA)->GetIndex()] += force / damping_constantA;
                
                // Only a ghost B also receives a force
                if (mIsGhostNode[p_edge->GetNodeGlobalIndex(nodeB)])
                {
                    drdt[ p_edge->GetNode(nodeB)->GetIndex()] -= force / damping_constantB;
                }
            }
        }
        elem_iter++;
    }
    
    
    // Here we divide all the forces on the nodes by a factor of two because
    // we looped over them all twice to deal with the boundaries above.
    for (unsigned i=0 ; i<mrMesh.GetNumAllNodes(); i++)
    {
        drdt[i]=drdt[i]/2.0;
    }
    
    return drdt;
}



/**
 * @return the x and y forces in this spring
 */
c_vector<double, 2> CryptSimulation2DPeriodic::CalculateForceInThisSpring(Element<2,2>*& rPElement,const unsigned& rNodeA,const unsigned& rNodeB)
{
    unsigned node_a_global_index = rPElement->GetNodeGlobalIndex(rNodeA);
    unsigned node_b_global_index = rPElement->GetNodeGlobalIndex(rNodeB);
    return CalculateForceBetweenNodes(node_a_global_index, node_b_global_index);
}

/**
 * @param rPEdge pointer to a boundary element
 * 
 * @return the x and y forces on node 0 of the boundary element
 */
c_vector<double, 2> CryptSimulation2DPeriodic::CalculateForceInThisBoundarySpring(BoundaryElement<1,2>*& rPEdge)
{
    unsigned node_a_global_index = rPEdge->GetNodeGlobalIndex(0);
    unsigned node_b_global_index = rPEdge->GetNodeGlobalIndex(1);
    return CalculateForceBetweenNodes(node_a_global_index, node_b_global_index);
}

/**
 * Calculates the force between two nodes.
 * 
 * Note that this assumes they are connected by a spring and should therefore
 * only be called by CalculateForceInThisSpring or CalculateForceInThisBoundarySpring
 * 
 * @param NodeAGlobalIndex
 * @param NodeBGlobalIndex
 * 
 * @return The force exerted on Node A by Node B.
 */
c_vector<double, 2> CryptSimulation2DPeriodic::CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
{
    c_vector<double, 2> unit_difference;
    c_vector<double, 2> node_a_location = mrMesh.GetNode(rNodeAGlobalIndex)->rGetLocation();
    c_vector<double, 2> node_b_location = mrMesh.GetNode(rNodeBGlobalIndex)->rGetLocation();
    unit_difference = mrMesh.GetVectorFromAtoB(node_a_location, node_b_location);   
    
    double distance_between_nodes = norm_2(unit_difference);
    
    unit_difference /= distance_between_nodes;
    
    double rest_length = 1.0;
    
    if ( (mCells.size()>0) &&  (!mIsGhostNode[rNodeAGlobalIndex])
                           &&  (!mIsGhostNode[rNodeBGlobalIndex]) )
    {
        double ageA = mCells[rNodeAGlobalIndex].GetAge();
        double ageB = mCells[rNodeBGlobalIndex].GetAge();
        if (ageA<1.0 && ageB<1.0 && fabs(ageA-ageB)<1e-6)
        {
            // Spring Rest Length Increases to normal rest length from 0.9 to normal rest length, 1.0, over 1 hour
#define COVERAGE_IGNORE
            rest_length=(0.1+0.9*ageA);
            assert(rest_length<=1.0);
#undef COVERAGE_IGNORE
        }
    }
    return mpParams->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}


/**
 * Moves each node to a new position for this timestep
 *
 * @param rDrDt the x and y force components on each node.
 */
void CryptSimulation2DPeriodic::UpdateNodePositions(const std::vector< c_vector<double, 2> >& rDrDt)
{
    for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
    {
    if (!mrMesh.GetNode(index)->IsDeleted())
        {
            Point<2> new_point = GetNewNodeLocation(index, rDrDt);   
            if (mFixedBoundaries)
            {
                c_vector<double, 2> node_position = mrMesh.GetNode(index)->rGetLocation();
                // All Boundaries x=0, x=crypt_width, y=0, y=crypt_length.
                if (   node_position[1]>0
                    && node_position[1]<mpParams->GetCryptLength()
                    && node_position[0]>0
                    && node_position[0]<mpParams->GetCryptWidth() )
                {
                    mrMesh.SetNode(index, new_point, false);
                }
            }
            else if (mCells.size()>0)
            {
                if (mWntIncluded)
                {   // A new Wnt feature - even stem cells can move as long as they don't go below zero.
                    if ( (new_point.rGetLocation()[1] < 0.0) && !mIsGhostNode[index])
                    {
                        new_point.rGetLocation()[1] = 0.0;
                    }
                    mrMesh.SetNode(index, new_point, false);
                }
                else
                {
                    // THE 'USUAL' SCENARIO move any node as long as it is not a real stem cell.
                    if (mCells[index].GetCellType()!=STEM || mIsGhostNode[index])
                    {   // if a cell wants to move below y<0 (most likely because it was
                        // just born from a stem cell), stop it doing so
                        if ( (new_point.rGetLocation()[1] < 0.0) && (!mIsGhostNode[index]))
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
                if (mrMesh.GetNode(index)->rGetLocation()[1]>0)
                {
                    mrMesh.SetNode(index, new_point, false);
                }
            }
        }
    }
}

Point<2> CryptSimulation2DPeriodic::GetNewNodeLocation(const unsigned& rOldNodeIndex, const std::vector< c_vector<double, 2> >& rDrDt)
{
    Point<2> new_point( mrMesh.GetNode(rOldNodeIndex)->rGetLocation()
                     + mDt*rDrDt[rOldNodeIndex]);
    return new_point;
}

/**
 * Change the state of cells
 *
 * At the moment this turns cells to be differentiated
 * dependent on a protein concentration when using the Wnt model.
 */
void CryptSimulation2DPeriodic::UpdateCellTypes()
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
    if (mWntIncluded)
    {   // Cycle through each cell
        for (unsigned i=0; i<mrMesh.GetNumNodes(); i++)
        {
            if (!mIsGhostNode[i])
            {
                // If we are in here the cell cycle model must be a WntCellCycleModel
                WntCellCycleModel *this_Wnt_model = static_cast<WntCellCycleModel*>(mCells[i].GetCellCycleModel());
                double betaCateninLevel = this_Wnt_model->GetProteinConcentrations()[6]+this_Wnt_model->GetProteinConcentrations()[7];
                // std::cout << "Cell " << i << ", beta-cat = " << betaCateninLevel << "\n" << std::endl;
                
                CryptCellType cell_type=TRANSIT;
                
                // For mitogenic stimulus of 6x10^-4 in Wnt equations
                if (betaCateninLevel < 0.4127)
                {
                    cell_type = DIFFERENTIATED;
                }
                if (betaCateninLevel >0.8)
                {
                    cell_type = STEM;
                }
                    
                mCells[i].SetCellType(cell_type);
            }
        }
    }
}



/**
 * This method actually calls the remesh command on the mesh.
 * It should only be called by the method above (ReMesh) which ensures
 * that periodic boundaries are handled properly.
 */
void CryptSimulation2DPeriodic::ReMesh()
{
    if(mReMesh)
    {
        NodeMap map(mrMesh.GetNumNodes());
        std::cout << "Remeshing \n"<< std::flush;
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
    }
}

/**
 * Set the timestep of the simulation
 */
void CryptSimulation2DPeriodic::SetDt(double dt)
{
    assert(dt>0);
    mDt=dt;
}

/**
 * Sets the end time and resets the timestep to be endtime/100
 */
void CryptSimulation2DPeriodic::SetEndTime(double endTime)
{
    assert(endTime>0);
    mEndTime=endTime;
}

void CryptSimulation2DPeriodic::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
}

/**
 *  Call this before Solve() to simulate cell growth after cell division.
 *  (will eventually become SetIncludeCellBirth() and then become the default)
 */
void CryptSimulation2DPeriodic::SetIncludeVariableRestLength()
{
    mIncludeVariableRestLength = true;
}

/**
 * Sets the maximum number of cells that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
void CryptSimulation2DPeriodic::SetMaxCells(unsigned maxCells)
{
    mMaxCells = maxCells;
    if (maxCells<mrMesh.GetNumAllNodes())
    {
        EXCEPTION("mMaxCells is less than the number of cells in the mesh.");
    }
}

/**
 * Sets the maximum number of elements that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
void CryptSimulation2DPeriodic::SetMaxElements(unsigned maxElements)
{
    mMaxElements = maxElements;
    if (maxElements<mrMesh.GetNumAllElements())
    {
        EXCEPTION("mMaxElements is less than the number of elements in the mesh.");
    }
}

/**
 * Call this before Solve() to fix the boundary of the mesh.
 * \todo figure out what this does!
 */
void CryptSimulation2DPeriodic::SetFixedBoundaries()
{
    mFixedBoundaries = true;
}

/**
 * The mesh should be surrounded by at least one layer of ghost nodes.  These are nodes which
 * do not correspond to a cell, but are necessary for remeshing (because the remesher tries to
 * create a convex hull of the set of nodes) and visualising purposes.  The mesh is passed into
 * the constructor and the class is told about the ghost nodes by using this method.
 */
void CryptSimulation2DPeriodic::SetGhostNodes(std::vector<unsigned> ghostNodeIndices)
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
    
    mCrypt.SetGhostNodes(mIsGhostNode);
}

/**
 * Get the mesh to be remeshed at every time step.
 */
void CryptSimulation2DPeriodic::SetReMeshRule(bool remesh)
{
    mReMesh = remesh;
}

/**
 * Set the simulation to run with no birth.
 */
void CryptSimulation2DPeriodic::SetNoBirth(bool nobirth)
{
    mNoBirth = nobirth;
}

/**
 * This automatically sets this to be a wnt dependent simulation.
 * You should supply cells with a wnt cell cycle...
 */
void CryptSimulation2DPeriodic::SetWntGradient(WntGradientType wntGradientType)
{
    mWntIncluded = true;
    mWntGradient = WntGradient(wntGradientType);
}

/**
 * Set this simulation to use a cell killer
 */
void CryptSimulation2DPeriodic::SetCellKiller(RandomCellKiller<2>* pCellKiller)
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
std::vector<MeinekeCryptCell> CryptSimulation2DPeriodic::GetCells()
{
    assert(mCells.size()>0);
    return mCells;
}

/**
 * Whether each node is a ghost or not.
 * \todo change this to return a const reference
 */
std::vector <bool> CryptSimulation2DPeriodic::GetGhostNodes()
{
    return mIsGhostNode;
}



/**
 * Get a node's location (ONLY FOR TESTING)
 *
 * @param the node index
 * @return the x and y co-ordinates of this node.
 */
std::vector<double> CryptSimulation2DPeriodic::GetNodeLocation(const unsigned& rNodeIndex)
{
    double x = mrMesh.GetNode(rNodeIndex)->rGetLocation()[0];
    double y = mrMesh.GetNode(rNodeIndex)->rGetLocation()[1];
    std::vector<double> location;
    location.push_back(x);
    location.push_back(y);
    return location;
}

/**
 * Main Solve method.
 *
 * Once CryptSimulation object has been set up, call this to run simulation
 */
void CryptSimulation2DPeriodic::Solve()
{ 
    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetDimensionalisedTime();
    std::cout << "Time at start of Solve Method = " << current_time << std::endl;
    
    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);
    std::cout << "num timesteps = " << num_time_steps << std::endl;
    if (current_time>0)//use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    else
    {
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    
    if (mOutputDirectory=="")
    {
        EXCEPTION("OutputDirectory not set");
    }
    
    double time_now = p_simulation_time->GetDimensionalisedTime();
    std::ostringstream time_string;
    time_string << time_now;
    
    std::string results_directory = mOutputDirectory +"/results_from_time_" + time_string.str();
    
    
    ///////////////////////////////////////////////////////////
    // Set up Simulation
    ///////////////////////////////////////////////////////////
    
    // Data writers for tabulated results data, used in tests
    // first construction clears out the folder
    ColumnDataWriter tabulated_node_writer(results_directory+"/tab_results", "tabulated_node_results",true);
    ColumnDataWriter tabulated_element_writer(results_directory+"/tab_results", "tabulated_element_results",false);
    
    node_writer_ids_t node_writer_ids;
    SetupNodeWriter(tabulated_node_writer, node_writer_ids);
    
    element_writer_ids_t element_writer_ids;
    SetupElementWriter(tabulated_element_writer, element_writer_ids);
    
    // This keeps track of when tabulated results were last output
    unsigned tabulated_output_counter = 0;
    
    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/vis_results/",false);
    out_stream p_node_file = output_file_handler.OpenOutputFile("results.viznodes");
    out_stream p_element_file = output_file_handler.OpenOutputFile("results.vizelements");
    out_stream p_setup_file = output_file_handler.OpenOutputFile("results.vizsetup");
    
    
    /* Age the cells to the correct time (cells set up with negative birth dates
     * to give some that are almost ready to divide).
     * 
     * TODO:For some strange reason this seems to take about 3 minutes for a realistic Wnt-Crypt.
     * Not sure why - when the same code was evaluated in a test it seemed almost instant.
     */
    if (!mCells.empty())
    {
        bool temp;
        for (unsigned i=0; i<mCells.size(); i++)
        {
            if (mIsGhostNode[i]) continue;
            //std::cout << "Preparing Cell "<< i << std::endl;
            Node<2> *p_our_node = mrMesh.GetNode(i);
            double y = p_our_node->rGetLocation()[1];
            std::vector<double> cell_cycle_influences;
            if (mWntIncluded)
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
    // Write initial conditions to file for the visualizer.
    WriteVisualizerSetupFile(*p_setup_file);
    WriteResultsToFiles(tabulated_node_writer, node_writer_ids,
                            tabulated_element_writer, element_writer_ids,
                            *p_node_file, *p_element_file,
                            false,
                            true);
    
    while (p_simulation_time->GetTimeStepsElapsed() < num_time_steps)
    {
        CheckIndicesAreInSync();
        
        mRemeshesThisTimeStep = 0; // To avoid infinite loops
        std::cout << "** TIME = " << p_simulation_time->GetDimensionalisedTime() << " **" << std::endl;
        
        mNumBirths += DoCellBirth();

        //  calculate node velocities
        std::vector<c_vector<double, 2> > drdt = CalculateVelocitiesOfEachNode();
        
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
                            tabulated_output_counter%80==0,
                            true);
                            
        tabulated_output_counter++;
    } // End main time loop
    
    // Write end state to tabulated files (not visualizer - this
    // is taken care of in the main loop).
    WriteResultsToFiles(tabulated_node_writer, node_writer_ids,
                        tabulated_element_writer, element_writer_ids,
                        *p_node_file, *p_element_file,
                        true,
                        false);
                        
    tabulated_node_writer.Close();
    tabulated_element_writer.Close();
}



/**
 * Saves the whole crypt simulation for restarting later.
 *
 * Puts it in the folder mOutputDirectory/archive/
 * and the file "2dCrypt_at_time_<SIMULATION TIME>.arch"
 *
 * First archives simulation time then the simulation itself.
 */
void CryptSimulation2DPeriodic::Save()
{
    SimulationTime* p_sim_time = SimulationTime::Instance();
    assert(p_sim_time->IsStartTimeSetUp());
    
    std::string archive_directory = mOutputDirectory + "/archive/";
    
    std::ostringstream time_stamp;
    time_stamp << p_sim_time->GetDimensionalisedTime();
    
    // create an output file handler in order to get the full path of the
    // archive directory. Note the false is so the handler doesn't clean
    // the directory
    OutputFileHandler handler(archive_directory, false);
    std::string archive_filename = handler.GetTestOutputDirectory() + "2dCrypt_at_time_"+time_stamp.str()+".arch";
    std::string mesh_filename = std::string("mesh_") + time_stamp.str();
    
    NodeMap map(mrMesh.GetNumAllNodes());
    mrMesh.ReMesh(map);
    
    // the false is so the directory isn't cleaned
    TrianglesMeshWriter<2,2> mesh_writer(archive_directory, mesh_filename, false);
    mesh_writer.WriteFilesUsingMesh(mrMesh);
    
    std::ofstream ofs(archive_filename.c_str());
    boost::archive::text_oarchive output_arch(ofs);
    
    // cast to const.
    const SimulationTime* p_simulation_time = SimulationTime::Instance();
    output_arch << *p_simulation_time;
    output_arch << static_cast<const CryptSimulation2DPeriodic&>(*this);
}

/**
 * Loads a saved crypt simulation to run further.
 *
 * @param rArchiveDirectory the name of the simulation to load
 * (specified originally by simulator.SetOutputDirectory("wherever"); )
 * @param rTimeStamp the time at which to load the simulation (this must
 * be one of the times at which the simulation.Save() was called)
 */
void CryptSimulation2DPeriodic::Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
{
    std::ostringstream time_stamp;
    time_stamp << rTimeStamp;
    
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    
    OutputFileHandler any_old_handler("",false);
    std::string test_output_directory = any_old_handler.GetTestOutputDirectory();
    
    std::string archive_filename = test_output_directory + rArchiveDirectory + "/archive/2dCrypt_at_time_"+time_stamp.str() +".arch";
    std::string mesh_filename = test_output_directory + rArchiveDirectory + "/archive/mesh_" + time_stamp.str();
    
    mrMesh.Clear();
    TrianglesMeshReader<2,2> mesh_reader(mesh_filename);
    mrMesh.ConstructFromMeshReader(mesh_reader);
    
    // Create an input archive
    std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
    boost::archive::text_iarchive input_arch(ifs);
    
    // read the archive
    assert(p_simulation_time->IsStartTimeSetUp());
    input_arch >> *p_simulation_time;
    input_arch >> *this;
    
    double time_now = p_simulation_time->GetDimensionalisedTime();
    std::ostringstream time_string;
    time_string << time_now;
    
    mOutputDirectory = rArchiveDirectory;
    
    if (mrMesh.GetNumNodes()!=mCells.size())
    {
        EXCEPTION(" Error in Load: number of nodes is not equal to number of cells.");
    }
}
