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
template<unsigned DIM> 
CryptSimulation2DPeriodic<DIM>::CryptSimulation2DPeriodic(ConformingTetrahedralMesh<DIM,DIM> &rMesh,
                                                     std::vector<MeinekeCryptCell> cells)
        : mrMesh(rMesh),
          mCells(cells),
          mCrypt(rMesh, mCells)
          
          
{ 
    
    mpParams = CancerParameters::Instance();
    
    mDt = 1.0/120.0;
    mEndTime = 0.0; // hours - this is set later on.
    
    srandom(0);
    mFixedBoundaries = false;
    mOutputDirectory = "";
    
    // Set up the ghost nodes bool list
    mIsGhostNode.resize(mrMesh.GetNumAllNodes());
    for (unsigned i=0; i<mIsGhostNode.size(); i++)
    {
        mIsGhostNode[i] = false;
    }
    mCrypt.SetGhostNodes(mIsGhostNode);
    
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
        #define COVERAGE_IGNORE
        EXCEPTION("Start time not set in simulation time singleton object");
        #undef COVERAGE_IGNORE
    }
    
}

/**
 * Free any memory allocated by the constructor
 */
template<unsigned DIM> 
CryptSimulation2DPeriodic<DIM>::~CryptSimulation2DPeriodic()
{
    SimulationTime::Destroy();
}

/**
* Define the variable identifiers in the data writer used to write node-based results.
*
* Uses mMaxCells to decide how many variables to define.
*/
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetupNodeWriter(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rVarIds)
{
    rVarIds.time = rNodeWriter.DefineUnlimitedDimension("Time","hours");
    
    rVarIds.types.resize(mMaxCells);
    rVarIds.position_id.resize(mMaxCells);
    
    // set up per-cell variables
    for (unsigned cell=0; cell<mMaxCells; cell++)
    {
        std::stringstream cell_type_var_name, cell_x_position_var_name, cell_y_position_var_name, cell_z_position_var_name;
        cell_type_var_name << "cell_type_" << cell;
        
        cell_x_position_var_name << "cell_x_position_" << cell;
        cell_y_position_var_name << "cell_y_position_" << cell;
        cell_z_position_var_name << "cell_z_position_" << cell;

        std::vector<std::string> cell_position_var_name_string;
        cell_position_var_name_string.push_back(cell_x_position_var_name.str());
        cell_position_var_name_string.push_back(cell_y_position_var_name.str());
        cell_position_var_name_string.push_back(cell_z_position_var_name.str());
        
        rVarIds.types[cell]=rNodeWriter.DefineVariable(cell_type_var_name.str(),"dimensionless");
        for(unsigned i=0; i<DIM; i++)
        {
            rVarIds.position_id[cell](i)=rNodeWriter.DefineVariable(cell_position_var_name_string[i],"rest_spring_length");
        }
    }
    
    rNodeWriter.EndDefineMode();
}

/**
 * Define the variable identifiers in the data writer used to write element-based results.
 *
 * Uses mMaxCells to decide how many variables to define.
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetupElementWriter(ColumnDataWriter& rElementWriter, element_writer_ids_t& rVarIds)
{
    rVarIds.time = rElementWriter.DefineUnlimitedDimension("Time","hours");
    
    // Set up columns for element writer
    rVarIds.node_id.resize(mMaxElements);
    
    for (unsigned elem_index = 0; elem_index<mMaxElements; elem_index++)
    {
        std::stringstream nodeA_var_name, nodeB_var_name, nodeC_var_name, nodeD_var_name;
        
        nodeA_var_name << "nodeA_" << elem_index;
        nodeB_var_name << "nodeB_" << elem_index;
        nodeC_var_name << "nodeC_" << elem_index;
        nodeD_var_name << "nodeD_" << elem_index;

        std::vector<std::string> node_var_name_string;
        node_var_name_string.push_back(nodeA_var_name.str());
        node_var_name_string.push_back(nodeB_var_name.str());
        node_var_name_string.push_back(nodeC_var_name.str());
        node_var_name_string.push_back(nodeD_var_name.str());
        
        for(unsigned i=0; i<DIM+1; i++)
        {
            rVarIds.node_id[elem_index](i) = rElementWriter.DefineVariable(node_var_name_string[i],"dimensionless");
        }
    }
    
    rElementWriter.EndDefineMode();
}

template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::WriteVisualizerSetupFile(std::ofstream& rSetupFile)
{
    assert(DIM==2); // this is 2d specific
    rSetupFile << "MeshWidth\t" << mrMesh.GetWidth(1u);// get furthest distance between nodes in the x-direciton
    rSetupFile.close();
}


template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::WriteResultsToFiles(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rNodeVarIds,
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
            const c_vector<double,DIM>& r_node_loc = mrMesh.GetNode(index)->rGetLocation();
            if (writeVisualizerResults)
            {
                for(unsigned i=0; i<DIM; i++)
                {
                    rNodeFile << r_node_loc[i] << " ";
                }
                rNodeFile << colour << " ";
            }
            if (writeTabulatedResults)
            {
                for(unsigned i=0; i<DIM; i++)
                {
                    rNodeWriter.PutVariable(rNodeVarIds.position_id[index](i), r_node_loc[i]);
                }
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
                for(unsigned i=0; i<DIM; i++)
                {
                    rElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i)<< " ";
                }
            }
            if (writeTabulatedResults)
            {
                for(unsigned i=0; i<DIM+1; i++)
                {
                    rElementWriter.PutVariable(rElementVarIds.node_id[elem_index](i), mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i));
                }
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
template<unsigned DIM>  
unsigned CryptSimulation2DPeriodic<DIM>::DoCellBirth()
{
    unsigned num_births_this_step = 0;
    if (!mNoBirth && !mCells.empty())
    {

        // Iterate over all cells, seeing if each one can be divided
        for (typename Crypt<DIM>::Iterator cell_iter = mCrypt.Begin();
             cell_iter != mCrypt.End();
             ++cell_iter)
        {
            MeinekeCryptCell& cell = *cell_iter;
            Node<2>* p_our_node = cell_iter.GetNode();

            assert(cell.GetNodeIndex() == p_our_node->GetIndex());
            
            // Check for this cell dividing
            // Construct any influences for the cell cycle...

            std::vector<double> cell_cycle_influences;
            if (mWntIncluded)
            {
                double y = p_our_node->rGetLocation()[1];
                double wnt_stimulus = mWntGradient.GetWntLevel(y);
                cell_cycle_influences.push_back(wnt_stimulus);
            }
            
            // CHECK if this cell is ready to divide - if so create a new cell etc.
            if (cell.ReadyToDivide(cell_cycle_influences))
            {
                // Create new cell
                MeinekeCryptCell new_cell = cell.Divide();
                std::cout << "Cell division at node " << cell.GetNodeIndex() << "\n";
            
                // Add a new node to the mesh
                unsigned parent_node_index = p_our_node->GetIndex();
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
                    mIsGhostNode.resize(mrMesh.GetNumNodes());
                    mIsGhostNode[new_node_index] = false;
                }
                num_births_this_step++;
            } // if (ready to divide)
        } // cell iteration loop
    } // if (simulation has cell birth)
    
    return num_births_this_step;
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
template<unsigned DIM> 
c_vector<double, DIM> CryptSimulation2DPeriodic<DIM>::CalculateDividingCellCentreLocations(unsigned node_index)
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
    
    
    if  (  (parent_coords[1]-random_vector[1] > 0.0)
        && (parent_coords[1]+random_vector[1] > 0.0))
    {   // We are not too close to the bottom of the crypt
        // add random vector to the daughter and take it from the parent location
        daughter_coords[0] = parent_coords[0]+random_vector[0];
        daughter_coords[1] = parent_coords[1]+random_vector[1];
        parent_coords[0] = parent_coords[0]-random_vector[0];
        parent_coords[1] = parent_coords[1]-random_vector[1];
    }
    else
    {   // Leave the parent where it is and move daughter in a positive direction
        // to ensure new cells are not born below y=0
        if (random_vector[1]>0.0)
        {
            daughter_coords[0] = parent_coords[0]+random_vector[0];
            daughter_coords[1] = parent_coords[1]+random_vector[1];            
        }
        else
        {
            daughter_coords[0] = parent_coords[0]-random_vector[0];
            daughter_coords[1] = parent_coords[1]-random_vector[1]; 
        }
    }
    assert(daughter_coords[1]>=0.0);// to make sure dividing cells stay in the crypt
    assert(parent_coords[1]>=0.0);// to make sure dividing cells stay in the crypt
    // set the parent to use this location
    mrMesh.SetNode(node_index, parent_coords, false);
    return daughter_coords;
}

/**
 *  Checks that the indices are in sync in the cells vector, ie that
 *  mCells[i].GetNodeIndex() is equal to i
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::CheckIndicesAreInSync()
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
template<unsigned DIM> 
unsigned CryptSimulation2DPeriodic<DIM>::DoCellRemoval()
{
    unsigned num_deaths=0;
    
    ///////////////////////////////////////////////////////////////////////////////////
    // Alternate method of sloughing.  Turns boundary nodes into ghost nodes.
    ///////////////////////////////////////////////////////////////////////////////////
    for (typename Crypt<DIM>::Iterator cell_iter = mCrypt.Begin();
         cell_iter != mCrypt.End();
         ++cell_iter)
    {
        
        double x = cell_iter.rGetLocation()[0];
        double y = cell_iter.rGetLocation()[1];
        
        double crypt_length=mpParams->GetCryptLength();
        double crypt_width=mpParams->GetCryptWidth();
        
        if ( (x>crypt_width) || (x<0.0) || (y>crypt_length))
        {
            mIsGhostNode[cell_iter.GetNode()->GetIndex()] = true;
            num_deaths++;
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
template<unsigned DIM>  
std::vector<c_vector<double, DIM> > CryptSimulation2DPeriodic<DIM>::CalculateVelocitiesOfEachNode()
{
    std::vector<c_vector<double, DIM> > drdt(mrMesh.GetNumAllNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i]=zero_vector<double>(DIM);
    }

    std::vector<std::vector<unsigned> > node_pairs_checked;

    ////////////////////////////////////////////////////////////////////
    // loop over element and for each one loop over its three edges
    ////////////////////////////////////////////////////////////////////
    for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
    {
        Element<DIM,DIM>* p_element = mrMesh.GetElement(elem_index);
        if (!p_element->IsDeleted())
        {
            for (unsigned k=0; k<DIM+1; k++)
            {
                unsigned nodeA = k;
                
                for(unsigned l=k+1; l<k+DIM+1; l++)
                {
                    unsigned nodeB = l%(DIM+1);
    
                    assert(!p_element->GetNode(nodeA)->IsDeleted());
                    assert(!p_element->GetNode(nodeB)->IsDeleted());
                                    
                    unsigned nodeA_global_index = p_element->GetNode(nodeA)->GetIndex();
                    unsigned nodeB_global_index = p_element->GetNode(nodeB)->GetIndex();
                    
                    // check whether we have already worked out the force between these two...
                    bool is_force_already_calculated = false;
                    
                    for (unsigned i=0; i<node_pairs_checked.size(); i++)
                    {
                        std::vector<unsigned> node_pair = node_pairs_checked[i];
                        if(node_pair[0]==nodeA_global_index || node_pair[1]==nodeA_global_index)
                        { 
                            // first node is in node_pair
                            if(node_pair[0]==nodeB_global_index || node_pair[1]==nodeB_global_index)
                            {
                                // both are in node_pair
                                is_force_already_calculated = true;
                                break;
                            } 
                        } 
                    } 
                    
                    if(!is_force_already_calculated)
                    {
                        
                        c_vector<double, DIM> force = CalculateForceInThisSpring(p_element,nodeA,nodeB);
                         
                        double damping_constantA = mpParams->GetDampingConstantNormal();
                        double damping_constantB = mpParams->GetDampingConstantNormal();
                        
                        if(!mCells.empty())
                        {
                            //note: at the moment the index into the mCells vector is the same
                            //as the node index. later this may not be the case, in which case
                            //the following assertion will trip. to deal with this, a map from 
                            //node index to cell will be needed
                            assert( mCells[nodeA_global_index].GetNodeIndex()==nodeA_global_index);
                            assert( mCells[nodeB_global_index].GetNodeIndex()==nodeB_global_index);
                            
                            if(   (mCells[nodeA_global_index].GetMutationState()==HEALTHY)
                               || (mCells[nodeA_global_index].GetMutationState()==APC_ONE_HIT))
                            {
                                damping_constantA = mpParams->GetDampingConstantNormal();
                            }
                            else
                            {
                                damping_constantA = mpParams->GetDampingConstantMutant();
                            }
                            
                            if(   (mCells[nodeB_global_index].GetMutationState()==HEALTHY)
                               || (mCells[nodeB_global_index].GetMutationState()==APC_ONE_HIT))
                            {
                                damping_constantB = mpParams->GetDampingConstantNormal();
                            }
                            else
                            {
                                damping_constantB = mpParams->GetDampingConstantMutant();
                            }
                        }
                        
                        if (!mIsGhostNode[nodeA_global_index])
                        {
                            drdt[ p_element->GetNode(nodeB)->GetIndex()] -= force / damping_constantB;
                            
                            if (!mIsGhostNode[nodeB_global_index])
                            {
                                drdt[ p_element->GetNode(nodeA)->GetIndex()] += force / damping_constantA;
                            }
                        }
                        else
                        {
                            drdt[ p_element->GetNode(nodeA)->GetIndex()] += force / damping_constantA;
                                
                            if (mIsGhostNode[nodeB_global_index])
                            {
                                drdt[ p_element->GetNode(nodeB)->GetIndex()] -= force / damping_constantB;
                            }
                        }
    
                        std::vector<unsigned> this_pair;
                        this_pair.push_back(nodeA_global_index);
                        this_pair.push_back(nodeB_global_index);
                        node_pairs_checked.push_back(this_pair);
                    }
                }
            }
        }
    }
  
    return drdt;
}



/**
 * @return the x and y forces in this spring
 */
template<unsigned DIM> 
c_vector<double, DIM> CryptSimulation2DPeriodic<DIM>::CalculateForceInThisSpring(Element<DIM,DIM>*& rPElement,const unsigned& rNodeA,const unsigned& rNodeB)
{
    assert(rNodeA!=rNodeB);
    unsigned node_a_global_index = rPElement->GetNodeGlobalIndex(rNodeA);
    unsigned node_b_global_index = rPElement->GetNodeGlobalIndex(rNodeB);
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
template<unsigned DIM> 
c_vector<double, DIM> CryptSimulation2DPeriodic<DIM>::CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
{
    assert(rNodeAGlobalIndex!=rNodeBGlobalIndex);
    c_vector<double, DIM> unit_difference;
    c_vector<double, DIM> node_a_location = mrMesh.GetNode(rNodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = mrMesh.GetNode(rNodeBGlobalIndex)->rGetLocation();
    
    // there is reason not to substract one position from the other (cyclidrical meshes). clever gary
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
            // Spring Rest Length Increases to normal rest length from 0.1 to normal rest length, 1.0, over 1 hour
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
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rDrDt)
{
    for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
    {
        if (!mrMesh.GetNode(index)->IsDeleted())
        {
            Point<DIM> new_point = GetNewNodeLocation(index, rDrDt);
            
            if(DIM==2)
            {
                // TODO: simplify/remove these 2d cases
                if (mFixedBoundaries)
                {
                    assert(DIM==2);
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
                }
                else
                {
                    // NOTE: once springs in their own tested class all crypt simulations should have cells, 
                    // so delete this:
                    // no cells, just fix any node on line y=0
                    if (mrMesh.GetNode(index)->rGetLocation()[1]>0)
                    {
                        mrMesh.SetNode(index, new_point, false);
                    }
                }
            }
            else
            {
                // 1d or 3d
                mrMesh.SetNode(index, new_point, false);
            }
        }
    }
}

template<unsigned DIM> 
Point<DIM> CryptSimulation2DPeriodic<DIM>::GetNewNodeLocation(const unsigned& rOldNodeIndex, const std::vector< c_vector<double, DIM> >& rDrDt)
{
    Point<DIM> new_point(   mrMesh.GetNode(rOldNodeIndex)->rGetLocation()
                          + mDt*rDrDt[rOldNodeIndex]);
    return new_point;
}

/**
 * Change the state of cells
 *
 * At the moment this turns cells to be differentiated
 * dependent on a protein concentration when using the Wnt model.
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::UpdateCellTypes()
{
    if (!mCells.empty())
    {
        /*
         * 
         * Designate cells as proliferating (transit) or
         * quiescent (differentiated) according to protein concentrations
         * 
         */
        for (typename Crypt<DIM>::Iterator cell_iter = mCrypt.Begin();
             cell_iter != mCrypt.End();
             ++cell_iter)
        {
            (*cell_iter).UpdateCellType();
        }
    }
}



/**
 * This method actually calls the remesh command on the mesh.
 * It should only be called by the method above (ReMesh) which ensures
 * that periodic boundaries are handled properly.
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::ReMesh()
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
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetDt(double dt)
{
    assert(dt>0);
    mDt=dt;
}

/**
 * Sets the end time and resets the timestep to be endtime/100
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetEndTime(double endTime)
{
    assert(endTime>0);
    mEndTime=endTime;
}

template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
}

/**
 * Sets the maximum number of cells that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
template<unsigned DIM>  
void CryptSimulation2DPeriodic<DIM>::SetMaxCells(unsigned maxCells)
{
    mMaxCells = maxCells;
    if (maxCells<mrMesh.GetNumAllNodes())
    {
        #define COVERAGE_IGNORE
        EXCEPTION("mMaxCells is less than the number of cells in the mesh.");
        #undef COVERAGE_IGNORE
    }
}

/**
 * Sets the maximum number of elements that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetMaxElements(unsigned maxElements)
{
    mMaxElements = maxElements;
    if (maxElements<mrMesh.GetNumAllElements())
    {
        #define COVERAGE_IGNORE
        EXCEPTION("mMaxElements is less than the number of elements in the mesh.");
        #undef COVERAGE_IGNORE
    }
}

/**
 * Call this before Solve() to fix the boundary of the mesh.
 * \todo figure out what this does!
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetFixedBoundaries()
{
    mFixedBoundaries = true;    // This is called by a nightly test.
}

/**
 * The mesh should be surrounded by at least one layer of ghost nodes.  These are nodes which
 * do not correspond to a cell, but are necessary for remeshing (because the remesher tries to
 * create a convex hull of the set of nodes) and visualising purposes.  The mesh is passed into
 * the constructor and the class is told about the ghost nodes by using this method.
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetGhostNodes(std::vector<unsigned> ghostNodeIndices)
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
 * Set whether the mesh should be remeshed at every time step.
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetReMeshRule(bool remesh)
{
    mReMesh = remesh;
}

/**
 * Set the simulation to run with no birth.
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetNoBirth(bool nobirth)
{
    mNoBirth = nobirth;
}

/**
 * This automatically sets this to be a wnt dependent simulation.
 * You should supply cells with a wnt cell cycle...
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetWntGradient(WntGradientType wntGradientType)
{
    mWntIncluded = true;
    mWntGradient = WntGradient(wntGradientType);
}

/**
 * Set this simulation to use a cell killer
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::SetCellKiller(RandomCellKiller<DIM>* pCellKiller)
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
template<unsigned DIM> 
std::vector<MeinekeCryptCell> CryptSimulation2DPeriodic<DIM>::GetCells()
{
    assert(mCells.size()>0);
    return mCells;
}

/**
 * Whether each node is a ghost or not.
 * \todo change this to return a const reference
 */
template<unsigned DIM> 
std::vector <bool> CryptSimulation2DPeriodic<DIM>::GetGhostNodes()
{
    return mIsGhostNode;
}



/**
 * Get a node's location (ONLY FOR TESTING)
 *
 * @param the node index
 * @return the co-ordinates of this node.
 */
template<unsigned DIM> 
std::vector<double> CryptSimulation2DPeriodic<DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for(unsigned i=0; i<DIM; i++)
    {
        location.push_back( mrMesh.GetNode(rNodeIndex)->rGetLocation()[i] );
    }
    return location;
}

/**
 * Main Solve method.
 *
 * Once CryptSimulation object has been set up, call this to run simulation
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::Solve()
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
        #define COVERAGE_IGNORE
        EXCEPTION("OutputDirectory not set");
        #undef COVERAGE_IGNORE
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
        for (typename Crypt<DIM>::Iterator cell_iter = mCrypt.Begin();
             cell_iter != mCrypt.End();
             ++cell_iter)
        {
            //std::cout << "Preparing Cell "<< i << std::endl;
            double y = cell_iter.rGetLocation()[1];
            std::vector<double> cell_cycle_influences;
            if (mWntIncluded)
            {
                double wnt_stimulus = mWntGradient.GetWntLevel(y);
                cell_cycle_influences.push_back(wnt_stimulus);
            }
            // We don't use the result; this call is just to force the cells to age to time 0,
            // running their cell cycle models to get there.
            (*cell_iter).ReadyToDivide(cell_cycle_influences);
        }
    }

    UpdateCellTypes();
    
    // Write initial conditions to file for the visualizer.
    WriteVisualizerSetupFile(*p_setup_file);
    WriteResultsToFiles(tabulated_node_writer, node_writer_ids,
                            tabulated_element_writer, element_writer_ids,
                            *p_node_file, *p_element_file,
                            false,
                            true);
    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////
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
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::Save()
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
    TrianglesMeshWriter<DIM,DIM> mesh_writer(archive_directory, mesh_filename, false);
    mesh_writer.WriteFilesUsingMesh(mrMesh);
    
    std::ofstream ofs(archive_filename.c_str());
    boost::archive::text_oarchive output_arch(ofs);
    
    // cast to const.
    const SimulationTime* p_simulation_time = SimulationTime::Instance();
    output_arch << *p_simulation_time;
    output_arch << static_cast<const CryptSimulation2DPeriodic<DIM>&>(*this);
}

/**
 * Loads a saved crypt simulation to run further.
 *
 * @param rArchiveDirectory the name of the simulation to load
 * (specified originally by simulator.SetOutputDirectory("wherever"); )
 * @param rTimeStamp the time at which to load the simulation (this must
 * be one of the times at which the simulation.Save() was called)
 */
template<unsigned DIM> 
void CryptSimulation2DPeriodic<DIM>::Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
{
    std::ostringstream time_stamp;
    time_stamp << rTimeStamp;
    
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    
    OutputFileHandler any_old_handler("",false);
    std::string test_output_directory = any_old_handler.GetTestOutputDirectory();
    
    std::string archive_filename = test_output_directory + rArchiveDirectory + "/archive/2dCrypt_at_time_"+time_stamp.str() +".arch";
    std::string mesh_filename = test_output_directory + rArchiveDirectory + "/archive/mesh_" + time_stamp.str();
    
    mrMesh.Clear();
    TrianglesMeshReader<DIM,DIM> mesh_reader(mesh_filename);
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
