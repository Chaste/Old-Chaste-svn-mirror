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
    mCylindrical = false;
    mNodesMoved = false;
    mRemeshesThisTimeStep = 0;
    mpCellKiller = NULL;
    mNumBirths = 0;
    mNumDeaths = 0;
    mPeriodicDivisionBuffer = 0;
    
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
 * Some divisions are postponed to work around issues with remeshing.  Specifically, only
 * one division on the periodic boundary is allowed within any 3 consecutive time steps.
 * This is controlled by the mPeriodicDivisionBuffer attribute.
 *
 * @return the number of births that occurred.
 */
unsigned CryptSimulation2DPeriodic::DoCellBirth()
{
    unsigned num_births = 0;
    if (mPeriodicDivisionBuffer>0)
    {
        mPeriodicDivisionBuffer--;
    }
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
            
            // Check if this cell is on the periodic boundary; there are more conditions if it is
            if (!skip && mPeriodicSides)
            {
                if ( (mLeftToRightBoundary.find(cell_index) != mLeftToRightBoundary.end() 
                      && mPeriodicDivisionBuffer > 0)
                    // Only allow one periodic cell division per
                    // timestep so that mesh can catch up with it.
                    // it will divide next timestep anyway 
                    || mRightToLeftBoundary.find(cell_index) != mRightToLeftBoundary.end()
                    // Only allow one periodic boundary to have divisions...
                   )
                {
                    skip=true;
                }
            }
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
                    if (mPeriodicSides && mLeftToRightBoundary.find(cell_index)!=mLeftToRightBoundary.end())
                    {
                        std::cout << "Periodic Division\n";
                        mPeriodicDivisionBuffer=3;
                        //Make sure the image cell knows it has just divided and aged a generation
                        unsigned image_cell_index = mLeftToRightBoundary[cell_index];
                        mCells[image_cell_index]=mCells[cell_index];
                        mCells[image_cell_index].SetNodeIndex(image_cell_index);
                    }
                    else
                    {
                        std::cout << "Cell division at node " << cell_index << "\n";
                    }
                    
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

                    if (mReMesh && mLeftToRightBoundary.find(cell_index)!=mLeftToRightBoundary.end())
                    {
                        ReMesh();
                    }
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
            
            if (!mPeriodicSides)
            {
                if ( (x>crypt_width) || (x<0.0) || (y>crypt_length))
                {
                    mIsGhostNode[node_index] = true;
                    num_deaths++;
                }
            }
            else
            {
                if (y>crypt_length)
                {
                    mIsGhostNode[p_node->GetIndex()] = true;
                    num_deaths++;
                    // And delete the periodic image if appropriate.(don't count as a death since it is an image)
                    if (mLeftToRightBoundary.find(node_index)!=mLeftToRightBoundary.end())
                    {
                        mIsGhostNode[mLeftToRightBoundary[node_index]]=true;
                    }
                    if (mRightToLeftBoundary.find(node_index)!=mRightToLeftBoundary.end())
                    {
                        mIsGhostNode[mRightToLeftBoundary[node_index]]=true;
                    }
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
Element<2,2>* CryptSimulation2DPeriodic::FindElementForBirth(Node<2>*& rpOurNode, unsigned cellIndex)
{
    // Pick a random element to start with
    Element<2,2>* p_element = mrMesh.GetElement(rpOurNode->GetNextContainingElementIndex());
    unsigned element_number = RandomNumberGenerator::Instance()->randMod(rpOurNode->GetNumContainingElements());
    for (unsigned j=0; j<element_number; j++)
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
            if (counter >= rpOurNode->GetNumContainingElements())
            {
                assert(mPeriodicNodes.find(rpOurNode->GetIndex())!=mPeriodicNodes.end());
                    // somehow every connecting element is a ghost element. quit to avoid infinite loop
                assert(mRightToLeftBoundary.find(cellIndex)==mRightToLeftBoundary.end());
                    // We already swapped; give up
                rpOurNode = mrMesh.GetNode(mLeftToRightBoundary[cellIndex]);
            }
            p_element = mrMesh.GetElement(rpOurNode->GetNextContainingElementIndex());
        }
    }
    while (is_ghost_element || is_periodic_element);
    
    return p_element;
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

    //////////////////////////////////////////////////////////////////
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
                unsigned node_a_index = p_element->GetNode(nodeA)->GetIndex();
                unsigned node_b_index = p_element->GetNode(nodeB)->GetIndex();
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
                

                if (!mPeriodicSides  // If A and B are in the left periodic edge ignore them (it will be handled by right edge)
                    || mLeftToRightBoundary.find(node_a_index)==mLeftToRightBoundary.end()
                    || mLeftToRightBoundary.find(node_b_index)==mLeftToRightBoundary.end() )
                {
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
    
    //assert(0);
    ///////////////
    //  sum forces for periodic nodes
    ////////////////
    if (mPeriodicSides)
    {
        // Add up the forces on paired up nodes
        // loop through size of left boundaries
        for (std::map <unsigned, unsigned>::iterator iterator = mLeftToRightBoundary.begin();
         iterator != mLeftToRightBoundary.end();
         iterator++)
        {
            c_vector<double,2> force = drdt[iterator->first] + drdt[iterator->second];
            drdt[iterator->first] = force;
            drdt[iterator->second] = force;
        }
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
                // move any node as long as it is not a real stem cell.
                if (mCells[index].GetCellType()!=STEM || mIsGhostNode[index])
                {
                    
                    // if a cell wants to move below y<0 (most likely because it was
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
    
    // Ensure no errors can creep in and move left nodes to same position as right ones
    if (mPeriodicSides)
    {
        for (unsigned i = 0; i < mLeftCryptBoundary.size();i++)
        {
            unsigned right_node_index = mRightCryptBoundary[i];
            unsigned left_node_index = mLeftCryptBoundary[i];
            c_vector<double,2> right_point = mrMesh.GetNode(right_node_index)->rGetLocation();
            Point<2> left_point;
            left_point.rGetLocation()[0] = right_point[0]-mpParams->GetCryptWidth();
            left_point.rGetLocation()[1] = right_point[1];
            mrMesh.SetNode(left_node_index, left_point, false);
            // Also force them to be the same cell
            // needed to synchronise cell cycle models (as R periodic cell cycle models are not run)...
            mCells[right_node_index]=mCells[left_node_index];
            mCells[right_node_index].SetNodeIndex(right_node_index);
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
                if (!(mCells[i].GetCellType()==STEM))
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
                    // For mitogenic stimulus of 5x10^-4 in Wnt equations
//                  	if(betaCateninLevel < 0.4954)
//                      {
//              	        cell_type = DIFFERENTIATED;
//                      }


                    mCells[i].SetCellType(cell_type);
                }
            }
        }
    }
}

/**
* Method to calculate the boundary of the crypt within the whole mesh ie the interface
* between normal and ghost nodes.
* 
* NB. Currently the crypt boundary is stored in two different data structures:
* (1) using standard vectors,
* (2) using other containers such as maps and sets.
* We plan to get rid of (1) eventually but keep it there while we refactor.
*/
void CryptSimulation2DPeriodic::CalculateCryptBoundary()
{
    assert(mPeriodicSides);
    mBoundary.clear();
    mPeriodicNodes.clear();
    mRightToLeftBoundary.clear();
    mLeftToRightBoundary.clear();
    
    double crypt_width=mpParams->GetCryptWidth();
    
    // Set all nodes as not being on the boundary...
    std::vector<bool> is_nodes_on_boundary(mIsGhostNode.size());
    
    for (unsigned i = 0 ; i < is_nodes_on_boundary.size() ; i++)
    {
        is_nodes_on_boundary[i]=false;
        mIsPeriodicNode[i]=false;
    }
    
    // Loop over elements and find boundary nodes of crypt
    for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
    {
        Element<2,2>* p_element = mrMesh.GetElement(elem_index);
        
        for (unsigned local_index=0; local_index<3; local_index++)
        {
           if (!mIsGhostNode[p_element->GetNode(local_index)->GetIndex()])
            {
                if ( mIsGhostNode[p_element->GetNode(0)->GetIndex()]
                     || mIsGhostNode[p_element->GetNode(1)->GetIndex()]
                     || mIsGhostNode[p_element->GetNode(2)->GetIndex()] )
                {
                    is_nodes_on_boundary[p_element->GetNode(local_index)->GetIndex()]= true;
                    mBoundary.insert(p_element->GetNode(local_index)->GetIndex());
                }
            }
        }
    }
    
    std::vector<unsigned> nodes_on_boundary;
    std::vector<unsigned> nodes_on_left_boundary;
    std::vector<unsigned> nodes_on_right_boundary;
    
    
    for (unsigned i = 0; i < is_nodes_on_boundary.size(); i++)
    {
    
        if (is_nodes_on_boundary[i])
        {
            nodes_on_boundary.push_back(i);
        }
    }
    assert(nodes_on_boundary.size()==mBoundary.size());
    
    
    for (unsigned i=0; i<nodes_on_boundary.size(); i++)
    {
        for (unsigned j=0; j<nodes_on_boundary.size(); j++)
        {
            // Check y positions are the same
            if (fabs(mrMesh.GetNode(nodes_on_boundary[i])->rGetLocation()[1]-mrMesh.GetNode(nodes_on_boundary[j])->rGetLocation()[1])<1e-3)
            {
                // Check x positions are crypt width apart.
                if (fabs(mrMesh.GetNode(nodes_on_boundary[j])->rGetLocation()[0]-mrMesh.GetNode(nodes_on_boundary[i])->rGetLocation()[0]-crypt_width)<1e-3)
                {
                    nodes_on_left_boundary.push_back(nodes_on_boundary[i]);
                    nodes_on_right_boundary.push_back(nodes_on_boundary[j]);
                    mIsPeriodicNode[nodes_on_boundary[i]]=true;
                    mIsPeriodicNode[nodes_on_boundary[j]]=true;
                }
            }
        }
    }
    
    for (std::set <unsigned>::iterator outer_iterator = mBoundary.begin();
         outer_iterator != mBoundary.end();
         outer_iterator++)
    {
        Node<2>* outer_node=mrMesh.GetNode(*outer_iterator);
        for (std::set <unsigned>::iterator inner_iterator = mBoundary.begin();
             inner_iterator != mBoundary.end();
             inner_iterator++)
        {
            Node<2>* inner_node=mrMesh.GetNode(*inner_iterator);
            // Check y positions are the same
            if (fabs(outer_node->rGetLocation()[1] - inner_node->rGetLocation()[1])<1e-3)
            {
                // Check x positions are crypt width apart.
                if (fabs(inner_node->rGetLocation()[0] - outer_node->rGetLocation()[0]-crypt_width)<1e-3)
                {
                    mLeftToRightBoundary[outer_node->GetIndex()]=inner_node->GetIndex();
                    mRightToLeftBoundary[inner_node->GetIndex()]=outer_node->GetIndex();
                    mPeriodicNodes.insert(outer_node->GetIndex());
                    mPeriodicNodes.insert(inner_node->GetIndex());
                }
            }
        } 
    }
    
    assert(nodes_on_boundary.size()==mBoundary.size());
    for (unsigned i=0; i<nodes_on_boundary.size(); i++)
    {
        assert(mBoundary.find(nodes_on_boundary[i])!=mBoundary.end());
        if (mIsPeriodicNode[nodes_on_boundary[i]])
        {
            assert(mPeriodicNodes.find(nodes_on_boundary[i])!=mPeriodicNodes.end());
        }
    }   
    
// DEBUG: write out nodes_on_right_boundary and mRightToLeftBoundary
//    for (unsigned i=0; i<nodes_on_right_boundary.size(); i++)
//    {
//        std::cout << nodes_on_right_boundary[i] << " " << std::flush;
//    }   
//    std::cout<< std::endl ;
//    for (std::map <unsigned, unsigned>::iterator iterator = mRightToLeftBoundary.begin();
//         iterator != mRightToLeftBoundary.end();
//         iterator++)
//    {
//        std::cout << iterator->first << " " << std::flush;
//    }
//    std::cout<< std::endl ;
    assert (nodes_on_left_boundary.size() == mLeftToRightBoundary.size() );
    assert (nodes_on_right_boundary.size() == mRightToLeftBoundary.size() );
    for (unsigned i=0; i<nodes_on_left_boundary.size(); i++)
    {
        assert( mLeftToRightBoundary[nodes_on_left_boundary[i]] == nodes_on_right_boundary[i]);
        assert( mRightToLeftBoundary[nodes_on_right_boundary[i]] == nodes_on_left_boundary[i]);
    }

    mLeftCryptBoundary =  nodes_on_left_boundary;
    mRightCryptBoundary =  nodes_on_right_boundary;
    mCryptBoundary =  nodes_on_boundary;
}

/**
 * This function detects when a remesh has caused a cell to
 * join or leave one of the periodic boundaries. It figures out where an image cell
 * will have to be placed to match it on the other side, and which two
 * periodic nodes have been upset.
 */
void CryptSimulation2DPeriodic::DetectNaughtyCellsAtPeriodicEdges()
{
    assert(mPeriodicSides);
    
    /* For each node see if it has broken onto a periodic boundary
     * If it has create an image node at the other side
     */
    unsigned i = 0;
    while (i<mCryptBoundary.size())
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
        
        if (mLeftToRightBoundary.find(our_node) != mLeftToRightBoundary.end()
            || mRightToLeftBoundary.find(our_node) != mRightToLeftBoundary.end() )
        {
            our_node_periodic=true;
        }
                
        
        if (!our_node_periodic)
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
                
                if (element1_node[0]==our_node || element1_node[1]==our_node || element1_node[2]==our_node)
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
                            
                            
                            
                            if (element2_node[0]==our_node || element2_node[1]==our_node || element2_node[2]==our_node)
                            {
                                unsigned ghost_node_element1 = 0;//it will never be this because this will be in bottom left corner of ghosts
                                unsigned ghost_node_element2 = 0;//it will never be this because this will be in bottom left corner of ghosts
                                
                                for (unsigned index = 0 ; index <3 ; index++)
                                {
                                    if (mIsGhostNode[element1_node[index]])
                                    {
                                        ghost_node_element1 = element1_node[index];
                                    }
                                }
                                for (unsigned index = 0 ; index <3 ; index++)
                                {
                                    if (mIsGhostNode[element2_node[index]])
                                    {
                                        ghost_node_element2 = element2_node[index];
                                    }
                                }
//		    						std::cout << "Examining two elements attached to node " << our_node << ".\n";
                                // this mucked up periodic births - new cell attached to more than one ghost node.
                                //if(ghost_node_element1==ghost_node_element2 && ghost_node_element1>0)
                                // if there is a ghost node in each of the two elements under consideration
                                if (ghost_node_element1>0 && ghost_node_element2>0)
                                {
                                    //std::cout << "node "<< our_node << " is attached to ghost nodes " << ghost_node_element1 << " and " << ghost_node_element2 << "\n";
                                    //Now check they are both attached to different periodic nodes
                                    unsigned other_node_element1 = 0;
                                    unsigned other_node_element2 = 0;
                                    
                                    bool left_nodes=false;
                                    for (unsigned j = 0 ; j<nodes_on_left_boundary.size() ; j++)
                                    {	//Cycle through periodic nodes and check for periodic neighbours
                                        unsigned periodic_node = nodes_on_left_boundary[j];
                                        if (element1_node[0]==periodic_node || element1_node[1]==periodic_node || element1_node[2]==periodic_node)
                                        {
                                            other_node_element1 = periodic_node;
                                            left_nodes=true;
                                        }
                                        if (element2_node[0]==periodic_node || element2_node[1]==periodic_node || element2_node[2]==periodic_node)
                                        {
                                            other_node_element2 = periodic_node;
                                            left_nodes=true;
                                        }
                                        periodic_node = nodes_on_right_boundary[j];
                                        if (element1_node[0]==periodic_node || element1_node[1]==periodic_node || element1_node[2]==periodic_node)
                                        {
                                            other_node_element1 = periodic_node;
                                        }
                                        if (element2_node[0]==periodic_node || element2_node[1]==periodic_node || element2_node[2]==periodic_node)
                                        {
                                            other_node_element2 = periodic_node;
                                        }
                                        if (other_node_element1!=other_node_element2 && other_node_element1>0 && other_node_element2>0 && !culprit_found)
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
                                            if (left_nodes)
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
            
            if (left_break)
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
            
            if (right_break)
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
void CryptSimulation2DPeriodic::RemoveSurplusCellsFromPeriodicBoundary()
{
    assert(mPeriodicSides);
    
    for (unsigned i=0 ; i<mOldLeftCryptBoundary.size() ; i++)
    {
        bool this_left_node_missing = true;
        bool this_right_node_missing = true;
        
        unsigned old_node_on_left_boundary = mOldLeftCryptBoundary[i];
        unsigned old_node_on_right_boundary = mOldRightCryptBoundary[i];
        
        for (unsigned j=0 ; j<mCryptBoundary.size() ; j++)
        {// search through the new crypt boundaries for this node
            if (old_node_on_left_boundary==mCryptBoundary[j])
            {
                this_left_node_missing=false;
            }
            if (old_node_on_right_boundary==mCryptBoundary[j])
            {
                this_right_node_missing=false;
            }
        }
        
        // Now this_<side>_node_missing is only true if the node has been internalised (or is a ghost already)
        if (this_left_node_missing && (!mIsGhostNode[old_node_on_left_boundary]))
        {	// The left node has been internalised (was periodic and is not a ghost) so the right node should be spooked
            mIsGhostNode[old_node_on_right_boundary]=true;
            std::cout << "Right Node " << old_node_on_right_boundary << " spooked\n";
            //mNodesMoved=true;
            CalculateCryptBoundary();
        }
        if (this_right_node_missing && (!mIsGhostNode[old_node_on_right_boundary]))
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
void CryptSimulation2DPeriodic::AddACellToPeriodicBoundary(unsigned original_node_index, double new_x, double new_y, std::vector< unsigned > periodic)
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
            for (unsigned i=0 ; i<3 ;i++)
            {
                if (mIsGhostNode[node[i]])
                {
                    bool already_flagged = false;
                    for (unsigned j=0; j<ghosts_on_node_0.size() ; j++)
                    {
                        if (ghosts_on_node_0[j]==node[i])
                        {
                            already_flagged = true;
                        }
                    }
                    if (!already_flagged)
                    {
                        ghosts_on_node_0.push_back(node[i]);
                    }
                }
            }
        }
        if	(node[0]==periodic_nodes[1] || node[1]==periodic_nodes[1] || node[2]==periodic_nodes[1])
        {	// if this element contains our first target node
            for (unsigned i=0 ; i<3 ;i++)
            {
                if (mIsGhostNode[node[i]])
                {
                    bool already_flagged = false;
                    for (unsigned j=0; j<ghosts_on_node_1.size() ; j++)
                    {
                        if (ghosts_on_node_1[j]==node[i])
                        {
                            already_flagged = true;
                        }
                    }
                    if (!already_flagged)
                    {
                        ghosts_on_node_1.push_back(node[i]);
                    }
                }
            }
        }
    }
    bool image_created = false;
    for (unsigned i=0 ; i< ghosts_on_node_0.size() ; i++)
    {
        for (unsigned j=0 ; j< ghosts_on_node_1.size() ; j++)
        {
            if (ghosts_on_node_0[i] == ghosts_on_node_1[j])
            {
                //std::cout << "Shared Node is " << ghosts_on_node_0[i] << "\n";
                // Move the ghost node to the correct position
                node_index = ghosts_on_node_0[i];
                image_created = true;
            }
        }
    }
    if (image_created==false)
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
            if (this_ghost_distance<nearest_distance)
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
            if (this_ghost_distance<nearest_distance)
            {
                nearest_distance=this_ghost_distance;
                nearest_node = ghosts_on_node_1[i];
            }
        }
        node_index = nearest_node;
        image_created=true;
    }
    
    if (image_created==false)
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
    mCells[node_index].SetNodeIndex(node_index);
    
}


/**
 * This method should be called by all other methods that require a re-mesh
 * it looks for changes in the periodic boundaries that could require
 * more remeshes and carries out the correct number.
 */
void CryptSimulation2DPeriodic::ReMesh()
{
    if (mReMesh)
    {
        CallReMesher();
        
        if (mPeriodicSides)
        {
            mOldLeftCryptBoundary = mLeftCryptBoundary;
            mOldRightCryptBoundary = mRightCryptBoundary;
            mOldCryptBoundary = mCryptBoundary;
            
            CalculateCryptBoundary();
            
            // do this once every time (sometimes mCryptBoundary does not change
            // but the connections between the nodes have changed)
            RemoveSurplusCellsFromPeriodicBoundary();
            DetectNaughtyCellsAtPeriodicEdges();
            
            while (mNodesMoved)
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
void CryptSimulation2DPeriodic::CallReMesher()
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

    mNodesMoved=false;
    mRemeshesThisTimeStep++;
    assert(mRemeshesThisTimeStep < 1000); //to avoid an infinite loop. If this ever throws try increasing it a bit.
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
 *  Call this before Solve() to set the boundary conditions
 * i.e. whether the left and right boundaries should be periodic
 */
void CryptSimulation2DPeriodic::SetPeriodicSides(bool periodicSides)
{
    mPeriodicSides = periodicSides;
}

/**
 * Call this before Solve() to tell the simulator that it is using
 * a Cylindrical2dMesh and not the old style of periodic sides.
 */
void CryptSimulation2DPeriodic::SetCylindrical()
{
    // TODO:
    // Make some kind of cast here to tell the 
    // simulator that the mesh is cylindrical??
    
    // If we are in here the mesh must be cylindrical
    //Cylindrical2dMesh mrMesh = static_cast<Cylindrical2dMesh>(mrMesh);
                    
    mPeriodicSides = false;
    mCylindrical = true;
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
 * Return the index of each node on the left periodic boundary
 * \todo change this to return a const reference
 */
std::vector<unsigned> CryptSimulation2DPeriodic::GetLeftCryptBoundary()
{
    return mLeftCryptBoundary;
}

/**
 * Return the index of each node on the right periodic boundary
 * \todo change this to return a const reference
 */
std::vector<unsigned> CryptSimulation2DPeriodic::GetRightCryptBoundary()
{
    return mRightCryptBoundary;
}

/**
 * Return the index of each node on the whole boundary
 * \todo change this to return a const reference
 */
std::vector<unsigned> CryptSimulation2DPeriodic::GetCryptBoundary()
{
    return mCryptBoundary;
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
    
    // Check some parameters for a periodic simulation
    if (mPeriodicSides)
    {
        CalculateCryptBoundary();
        if (mLeftCryptBoundary.size()<1 || mRightCryptBoundary.size()<1)
        {
            EXCEPTION("Periodic Simulation but mesh is not periodic\nIf you want a non-periodic simulation use SetPeriodicSides(false)");
        }
        if (!mReMesh)
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
