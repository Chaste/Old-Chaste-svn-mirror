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
#include "ColumnDataWriter.hpp"
#include "MeinekeCryptCellTypes.hpp"
//#include "WntCellCycleModel.hpp"
//#include "WntGradient.hpp"
//TODO: This should become abstract
#include "RandomCellKiller.hpp"
#include "OutputFileHandler.hpp"
#include "SpheroidSimulation3D.hpp"
#include "TrianglesMeshReader.hpp"

/** Constructor
 *  @param cells is defaulted to the empty vector, in which case SetIncludeRandomBirth()
 *  should be called for any birth to happen.
 */
SpheroidSimulation3D::SpheroidSimulation3D(ConformingTetrahedralMesh<3,3> &rMesh,
                                                     std::vector<MeinekeCryptCell> cells)
        : mrMesh(rMesh),
          mCells(cells)
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
    mNodesMoved=false;
    mRemeshesThisTimeStep=0;
    mpCellKiller=NULL;
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
SpheroidSimulation3D::~SpheroidSimulation3D()
{
    SimulationTime::Destroy();
}

/**
* Define the variable identifiers in the data writer used to write node-based results.
*
* Uses mMaxCells to decide how many variables to define.
*/
void SpheroidSimulation3D::SetupNodeWriter(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rVarIds)
{
    rVarIds.time = rNodeWriter.DefineUnlimitedDimension("Time","hours");
    
    rVarIds.types.resize(mMaxCells);
    rVarIds.x_positions.resize(mMaxCells);
    rVarIds.y_positions.resize(mMaxCells);
    rVarIds.z_positions.resize(mMaxCells);
    
    // set up per-cell variables
    for (unsigned cell=0; cell<mMaxCells; cell++)
    {
        std::stringstream cell_type_var_name, cell_x_position_var_name, cell_y_position_var_name, cell_z_position_var_name;
        cell_type_var_name << "cell_type_" << cell;
        cell_x_position_var_name << "cell_x_position_" << cell;
        cell_y_position_var_name << "cell_y_position_" << cell;
        cell_z_position_var_name << "cell_z_position_" << cell;
        rVarIds.types[cell]=rNodeWriter.DefineVariable(cell_type_var_name.str(),"dimensionless");
        rVarIds.x_positions[cell]=rNodeWriter.DefineVariable(cell_x_position_var_name.str(),"rest_spring_length");
        rVarIds.y_positions[cell]=rNodeWriter.DefineVariable(cell_y_position_var_name.str(),"rest_spring_length");
        rVarIds.z_positions[cell]=rNodeWriter.DefineVariable(cell_z_position_var_name.str(),"rest_spring_length");
    }
    
    rNodeWriter.EndDefineMode();
}

/**
 * Define the variable identifiers in the data writer used to write element-based results.
 *
 * Uses mMaxCells to decide how many variables to define.
 */
void SpheroidSimulation3D::SetupElementWriter(ColumnDataWriter& rElementWriter, element_writer_ids_t& rVarIds)
{
    rVarIds.time = rElementWriter.DefineUnlimitedDimension("Time","hours");
    
    // Set up columns for element writer
    rVarIds.nodeAs.resize(mMaxElements);
    rVarIds.nodeBs.resize(mMaxElements);
    rVarIds.nodeCs.resize(mMaxElements);
    rVarIds.nodeDs.resize(mMaxElements);
    
    for (unsigned elem_index = 0; elem_index<mMaxElements; elem_index++)
    {
        std::stringstream nodeA_var_name;
        std::stringstream nodeB_var_name;
        std::stringstream nodeC_var_name;
        std::stringstream nodeD_var_name;
        
        nodeA_var_name << "nodeA_" << elem_index;
        nodeB_var_name << "nodeB_" << elem_index;
        nodeC_var_name << "nodeC_" << elem_index;
        nodeD_var_name << "nodeD_" << elem_index;
        
        rVarIds.nodeAs[elem_index] = rElementWriter.DefineVariable(nodeA_var_name.str(),"dimensionless");
        rVarIds.nodeBs[elem_index] = rElementWriter.DefineVariable(nodeB_var_name.str(),"dimensionless");
        rVarIds.nodeCs[elem_index] = rElementWriter.DefineVariable(nodeC_var_name.str(),"dimensionless");
        rVarIds.nodeDs[elem_index] = rElementWriter.DefineVariable(nodeD_var_name.str(),"dimensionless");
    }
    
    rElementWriter.EndDefineMode();
}


void SpheroidSimulation3D::WriteResultsToFiles(ColumnDataWriter& rNodeWriter, node_writer_ids_t& rNodeVarIds,
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
            const c_vector<double,3>& r_node_loc = mrMesh.GetNode(index)->rGetLocation();
            if (writeVisualizerResults)
            {
                rNodeFile << r_node_loc[0] << " "<< r_node_loc[1] << " " << r_node_loc[2] << " " << colour << " ";
            }
            if (writeTabulatedResults)
            {
                rNodeWriter.PutVariable(rNodeVarIds.x_positions[index], r_node_loc[0]);
                rNodeWriter.PutVariable(rNodeVarIds.y_positions[index], r_node_loc[1]);
                rNodeWriter.PutVariable(rNodeVarIds.z_positions[index], r_node_loc[2]);
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
                rElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0)<< " " << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1)<< " "<< mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2)<< " "<< mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(3)<< " ";
            }
            if (writeTabulatedResults)
            {
                rElementWriter.PutVariable(rElementVarIds.nodeAs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(0));
                rElementWriter.PutVariable(rElementVarIds.nodeBs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(1));
                rElementWriter.PutVariable(rElementVarIds.nodeCs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(2));                
                rElementWriter.PutVariable(rElementVarIds.nodeDs[elem_index], mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(3));
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
unsigned SpheroidSimulation3D::DoCellBirth()
{
    unsigned num_births = 0;
    if (mPeriodicDivisionBuffer>0)
    {
        mPeriodicDivisionBuffer--;
    }
    if (!mNoBirth && !mCells.empty())
    {
        // Iterate over all cells, seeing if each one can be divided
        for (unsigned cell_index=0; cell_index<mCells.size(); cell_index++)
        {
            
            //unsigned node_index = mCells[cell_index].GetNodeIndex();
            //assert(node_index==cell_index);
            // To do: Uncomment code above
            unsigned node_index=cell_index ;
            bool skip = false; // Whether to not try dividing this cell
            if (mrMesh.GetNode(node_index)->IsDeleted()) skip=true; // Skip deleted cells
            //if (mrMesh.GetNode(cell_index)->IsDead()) skip=true; // Skip dead cells
            if (mIsGhostNode[cell_index]) skip=true; // Skip Ghost nodes
            
            if (skip) continue;
            
            // Construct any influences for the cell cycle...
            Node<3> *p_our_node = mrMesh.GetNode(cell_index);
            std::vector<double> cell_cycle_influences;
            
            // CHECK if this cell is ready to divide - if so create a new cell etc.
            if (mCells[cell_index].ReadyToDivide(cell_cycle_influences))
            {
                // Create new cell
                MeinekeCryptCell new_cell = mCells[cell_index].Divide();
                                
                // Add new node to mesh
                Element<3,3>* p_element = FindElementForBirth(p_our_node, cell_index);
                                                              
                std::cout << "New cell being intoduced into element with nodes \n";
                std::cout << p_element->GetNodeGlobalIndex(0) << "\t" << p_element->GetNodeGlobalIndex(1) << "\t" <<p_element->GetNodeGlobalIndex(2) << "\t" << p_element->GetNodeGlobalIndex(3) << "\n";
                double x = p_our_node->rGetLocation()[0];
                double y = p_our_node->rGetLocation()[1];
                double z = p_our_node->rGetLocation()[2];
                
                double x_centroid = (1.0/4.0)*(p_element->GetNode(0)->rGetLocation()[0]
                                               +  p_element->GetNode(1)->rGetLocation()[0]
                                               +  p_element->GetNode(2)->rGetLocation()[0]
                                               +  p_element->GetNode(3)->rGetLocation()[0] );
                                               
                double y_centroid = (1.0/4.0)*(p_element->GetNode(0)->rGetLocation()[1]
                                               +  p_element->GetNode(1)->rGetLocation()[1]
                                               +  p_element->GetNode(2)->rGetLocation()[1] 
                                               +  p_element->GetNode(3)->rGetLocation()[1] );
                                               
                double z_centroid = (1.0/4.0)*(p_element->GetNode(0)->rGetLocation()[2]
                                               +  p_element->GetNode(1)->rGetLocation()[2]
                                               +  p_element->GetNode(2)->rGetLocation()[2] 
                                               +  p_element->GetNode(3)->rGetLocation()[2]);
                                               
                                               
                // check the new point is in the tetrahedron
                double distance_from_node_to_centroid =  sqrt(  (x_centroid - x)*(x_centroid - x)
                                                                + (y_centroid - y)*(y_centroid - y)+ (z_centroid - z)*(z_centroid - z) );
                                                                
                // we assume the new cell is a distance 0.1 away from the old.
                // however, to avoid crashing in usual situations we check this
                // new position is actually in the tetrahedron being refined.
                // TODO: Check this is correct!
                double distance_of_new_cell_from_parent = 0.1;
                if (distance_from_node_to_centroid < (2.0/3.0)*0.1) // for tetrahedra is 2/3 right?
                {
                    #define COVERAGE_IGNORE
                    distance_of_new_cell_from_parent = (3.0/2.0)*distance_from_node_to_centroid;
                    #undef COVERAGE_IGNORE
                }
                
                double new_x_value = x + distance_of_new_cell_from_parent*(x_centroid-x);
                double new_y_value = y + distance_of_new_cell_from_parent*(y_centroid-y);
                double new_z_value = z + distance_of_new_cell_from_parent*(z_centroid-z);
                
                //std::cout << "Parent node at x = " << x << "  y = " << y << "\n";
                //std::cout << "Daughter node at x = " << new_x_value << "  y = " << new_y_value << "\n";
                
                Point<3> new_point(new_x_value, new_y_value, new_z_value);
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
                if (mrMesh.GetNumNodes() > mIsGhostNode.size())
                {
                    #define COVERAGE_IGNORE
                    mIsGhostNode.resize(mrMesh.GetNumNodes());
                    mIsGhostNode[new_node_index] = false;
                    #undef COVERAGE_IGNORE
                }
                num_births++;
                //std::cout<< " num_births=" << num_births <<std::endl<< std::flush;
                if (mReMesh)
                {
                    ReMesh();
                }
            } // if (ready to divide)
        } // cell iteration loop
    } // if (simulation has cell birth)
    
    return num_births;
}

/**
 *  Checks that the indices are in sync in the cells vector, ie that
 *  mCells[i].GetNodeIndex() is equal to i
 */
void SpheroidSimulation3D::CheckIndicesAreInSync()
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
unsigned SpheroidSimulation3D::DoCellRemoval()
{
    unsigned num_deaths=0;
    
    ///////////////////////////////////////////////////////////////////////////////////
    // Alternate method of sloughing.  Turns boundary nodes into ghost nodes.
    ///////////////////////////////////////////////////////////////////////////////////
    for (unsigned i=0; i<mrMesh.GetNumNodes(); i++)
    {
        Node<3> *p_node = mrMesh.GetNode(i);
        if (!p_node->IsDeleted())
        {
            double x = p_node->rGetLocation()[0];
            double y = p_node->rGetLocation()[1];
            
            double crypt_length=mpParams->GetCryptLength();
            double crypt_width=mpParams->GetCryptWidth();
            
            if (!mPeriodicSides)
            {
                if ( (x>crypt_width) || (x<0.0) || (y>crypt_length))
                {
                    mIsGhostNode[p_node->GetIndex()] = true;
                    num_deaths++;
                    //std::cout<< "num_deaths=" << num_deaths <<std::endl<< std::flush;
                }
            }
            else
            {
                if (y>crypt_length)
                {
                    mIsGhostNode[p_node->GetIndex()] = true;
                    num_deaths++;
                    if (mPeriodicSides)
                    {
                        // And delete the periodic image if appropriate.(don't count as a death since it is an image)
                        for (unsigned j=0 ; j<mLeftCryptBoundary.size() ; j++)
                        {
                            if ((unsigned)p_node->GetIndex()==mLeftCryptBoundary[j])
                            {
                                mIsGhostNode[mRightCryptBoundary[j]]=true;
                            }
                            if ((unsigned)p_node->GetIndex()==mRightCryptBoundary[j])
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
 * @return  a (pointer to a) suitable element
 */
Element<3,3>* SpheroidSimulation3D::FindElementForBirth(Node<3>*& rpOurNode, unsigned cellIndex)
{
    // Pick a random element to start with
    Element<3,3>* p_element = mrMesh.GetElement(rpOurNode->GetNextContainingElementIndex());
    unsigned element_number = RandomNumberGenerator::Instance()->randMod(rpOurNode->GetNumContainingElements());
    for (unsigned j=0; j<element_number; j++)
    {
        p_element = mrMesh.GetElement(rpOurNode->GetNextContainingElementIndex());
    }
    
    unsigned counter = 0; // how many elements we've checked for suitability
    bool is_ghost_element;
    do
    {
        // A ghost element has at least 1 ghost node
        is_ghost_element = (   mIsGhostNode[p_element->GetNodeGlobalIndex(0)]
                               || mIsGhostNode[p_element->GetNodeGlobalIndex(1)]
                               || mIsGhostNode[p_element->GetNodeGlobalIndex(2)] );
        
        if (is_ghost_element)
        {
            // This element isn't suitable
            counter++;
            if (counter >= rpOurNode->GetNumContainingElements())
            {
                
                // somehow every connecting element is a ghost element. quit to
                // avoid infinite loop
                #define COVERAGE_IGNORE
                assert(0);
                #undef COVERAGE_IGNORE
                
            }
            p_element = mrMesh.GetElement(rpOurNode->GetNextContainingElementIndex());
        }
    }
    while (is_ghost_element);
    
    return p_element;
}

/**
 *Calculates the forces on each node
 *
 * @return drdt the x, y and z force components on each node
*/
std::vector<std::vector<double> > SpheroidSimulation3D::CalculateVelocitiesOfEachNode()
{
    std::vector<std::vector<double> > drdt(mrMesh.GetNumAllNodes());
    for (unsigned i=0; i<mrMesh.GetNumAllNodes(); i++)
    {
        drdt[i].resize(3);
    }
    
    std::vector<std::vector<unsigned> > node_pairs_checked;
    //////////////////////////////////////////////////////////////////
    // loop over element and for each one loop over its three edges
    ////////////////////////////////////////////////////////////////////
    for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
    {
        Element<3,3>* p_element = mrMesh.GetElement(elem_index);
        
        for (unsigned k=0; k<4; k++)
        {
            unsigned nodeA = k;
            
            for (unsigned l=k+1; l<k+4; l++)
            {
                unsigned nodeB = l%4;
                unsigned nodeA_global_index = p_element->GetNode(nodeA)->GetIndex();
                unsigned nodeB_global_index = p_element->GetNode(nodeB)->GetIndex();
            
            
                // check whether we have already worked out the force between these two...
                bool is_force_already_calculated = false;
                
                for (unsigned i=0 ; i<node_pairs_checked.size() ; i++)
                {
                   std::vector<unsigned> node_pair = node_pairs_checked[i];
                   if(node_pair[0]==nodeA_global_index || node_pair[1]==nodeA_global_index)
                   { // first node is in node_pair
                        if(node_pair[0]==nodeB_global_index || node_pair[1]==nodeB_global_index)
                        { // both are in node_pair
                            is_force_already_calculated = true;
                        } 
                   } 
                }                
            
                if (!is_force_already_calculated)
                {
                    c_vector<double, 3> force = CalculateForceInThisSpring(p_element,nodeA,nodeB);
                    
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
                    
                    // Assume that if both nodes are real, or both are ghosts, then they both
                    // exert forces on each other, but if one is real and one is ghost then
                    // the real node exerts a force on the ghost node, but the ghost node
                    // does NOT exert a force on the real node.
                    if (!mIsGhostNode[p_element->GetNodeGlobalIndex(nodeA)])
                    {
                        drdt[nodeB_global_index][0] -= force(0) / damping_constantB;
                        drdt[nodeB_global_index][1] -= force(1) / damping_constantB;
                        drdt[nodeB_global_index][2] -= force(2) / damping_constantB;
                        
                        if (!mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)])
                        {
                            drdt[nodeA_global_index][0] += force(0) / damping_constantA;
                            drdt[nodeA_global_index][1] += force(1) / damping_constantA;
                            drdt[nodeA_global_index][2] += force(2) / damping_constantA;
                        }
                    }
                    else
                    {
                        drdt[nodeA_global_index][0] += force(0) / damping_constantA;
                        drdt[nodeA_global_index][1] += force(1) / damping_constantA;
                        drdt[nodeA_global_index][1] += force(2) / damping_constantA;
                        
                        if (mIsGhostNode[p_element->GetNodeGlobalIndex(nodeB)])
                        {
                            drdt[nodeB_global_index][0] -= force(0) / damping_constantB;
                            drdt[nodeB_global_index][1] -= force(1) / damping_constantB;
                            drdt[nodeB_global_index][2] -= force(2) / damping_constantB;
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

    return drdt;
}

/**
 * @return the x and y forces in this spring
 */
c_vector<double, 3> SpheroidSimulation3D::CalculateForceInThisSpring(Element<3,3>*& rPElement,const unsigned& rNodeA,const unsigned& rNodeB)
{
    c_vector<double, 3> force;
    c_vector<double, 3> unit_difference;
    unit_difference(0)=rPElement->GetNodeLocation(rNodeB,0)-rPElement->GetNodeLocation(rNodeA,0);
    unit_difference(1)=rPElement->GetNodeLocation(rNodeB,1)-rPElement->GetNodeLocation(rNodeA,1);
    unit_difference(2)=rPElement->GetNodeLocation(rNodeB,2)-rPElement->GetNodeLocation(rNodeA,2);
    
    double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1)+unit_difference(2)*unit_difference(2));
    unit_difference=unit_difference/distance_between_nodes;
    double rest_length = 1.0;
    
    if ( (mCells.size()>0) &&  (!mIsGhostNode[rPElement->GetNodeGlobalIndex(rNodeA)]) && (!mIsGhostNode[rPElement->GetNodeGlobalIndex(rNodeB)]) )
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
    return force = mpParams->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}



/**
 * @param rPEdge pointer to a boundary element
 * 
 * @return the x and y forces on node 0 of the boundary element
 */
c_vector<double, 3> SpheroidSimulation3D::CalculateForceInThisBoundarySpring(BoundaryElement<2,3>*& rPEdge)
{
    const unsigned nodeA = 0;
    const unsigned nodeB = 1;
    c_vector<double, 3> force;
    c_vector<double, 3> unit_difference;
    unit_difference(0)=rPEdge->GetNodeLocation(nodeB,0)-rPEdge->GetNodeLocation(nodeA,0);
    unit_difference(1)=rPEdge->GetNodeLocation(nodeB,1)-rPEdge->GetNodeLocation(nodeA,1);
    unit_difference(2)=rPEdge->GetNodeLocation(nodeB,2)-rPEdge->GetNodeLocation(nodeA,2);
    double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1)+unit_difference(2)*unit_difference(2));
    
    unit_difference=unit_difference/distance_between_nodes;
    
    double rest_length = 1.0;
    
    if ( (mCells.size()>0) &&  (!mIsGhostNode[rPEdge->GetNodeGlobalIndex(nodeA)]) && (!mIsGhostNode[rPEdge->GetNodeGlobalIndex(nodeB)]) )
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
    
    return force = mpParams->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}

    

/**
 * Moves each node to a new position for this timestep
 *
 * @param rDrDt the x and y force components on each node.
 */
void SpheroidSimulation3D::UpdateNodePositions(const std::vector< std::vector<double> >& rDrDt)
{
    for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
    {
        if (!mrMesh.GetNode(index)->IsDeleted())
        {
            Point<3> new_point = GetNewNodeLocation(index,rDrDt);
            mrMesh.SetNode(index, new_point, false);
        }         
    }
}

Point<3> SpheroidSimulation3D::GetNewNodeLocation(const unsigned& rOldNodeIndex, const std::vector< std::vector<double> >& rDrDt)
{
    c_vector<double,3> old_point = mrMesh.GetNode(rOldNodeIndex)->rGetLocation();
    Point<3> new_point;
    
    // Euler style update to node position
    new_point.rGetLocation()[0] = old_point[0] + mDt*rDrDt[rOldNodeIndex][0];
    new_point.rGetLocation()[1] = old_point[1] + mDt*rDrDt[rOldNodeIndex][1];
    new_point.rGetLocation()[2] = old_point[2] + mDt*rDrDt[rOldNodeIndex][2];
    
    return new_point;
}

/**
 * Change the state of cells
 *
 * At the moment this turns cells to be differentiated
 * dependent on a protein concentration when using the Wnt model.
 */
void SpheroidSimulation3D::UpdateCellTypes()
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
//                if (!(mCells[i].GetCellType()==STEM))
//                {
//                    // If we are in here the cell cycle model must be a WntCellCycleModel
//                    WntCellCycleModel *this_Wnt_model = static_cast<WntCellCycleModel*>(mCells[i].GetCellCycleModel());
//                    double betaCateninLevel = this_Wnt_model->GetProteinConcentrations()[6]+this_Wnt_model->GetProteinConcentrations()[7];
//                    // std::cout << "Cell " << i << ", beta-cat = " << betaCateninLevel << "\n" << std::endl;
//                    
//                    CryptCellType cell_type=TRANSIT;
//                    
//                    // For mitogenic stimulus of 6x10^-4 in Wnt equations
//                    if (betaCateninLevel < 0.4127)
//                    {
//                        cell_type = DIFFERENTIATED;
//                    }
//                    // For mitogenic stimulus of 5x10^-4 in Wnt equations
////                  	if(betaCateninLevel < 0.4954)
////                      {
////              	        cell_type = DIFFERENTIATED;
////                      }
//
//
//                    mCells[i].SetCellType(cell_type);
//                }
            }
        }
    }
}

/**
 * This method should be called by all other methods that require a re-mesh
 * it looks for changes in the periodic boundaries that could require
 * more remeshes and carries out the correct number.
 */
void SpheroidSimulation3D::ReMesh()
{
    if (mReMesh)
    {
        CallReMesher();
    }
}

/**
 * This method actually calls the remesh command on the mesh.
 * It should only be called by the method above (ReMesh) which ensures
 * that periodic boundaries are handled properly.
 */
void SpheroidSimulation3D::CallReMesher()
{    
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
 * Set the timestep of the simulation
 */
void SpheroidSimulation3D::SetDt(double dt)
{
    assert(dt>0);
    mDt=dt;
}

/**
 * Sets the end time and resets the timestep to be endtime/100
 */
void SpheroidSimulation3D::SetEndTime(double endTime)
{
    assert(endTime>0);
    mEndTime=endTime;
}

void SpheroidSimulation3D::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
}

/**
 *  Call this before Solve() to simulate cell growth after cell division.
 *  (will eventually become SetIncludeCellBirth() and then become the default)
 */
void SpheroidSimulation3D::SetIncludeVariableRestLength()
{
    mIncludeVariableRestLength = true;
}

/**
 * Sets the maximum number of cells that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
void SpheroidSimulation3D::SetMaxCells(unsigned maxCells)
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
void SpheroidSimulation3D::SetMaxElements(unsigned maxElements)
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
void SpheroidSimulation3D::SetFixedBoundaries()
{
    mFixedBoundaries = true;
}

/**
 *  Call this before Solve() to set the boundary conditions
 * i.e. whether the left and right boundaries should be periodic
 */
void SpheroidSimulation3D::SetPeriodicSides(bool periodicSides)
{
    mPeriodicSides = periodicSides;
}

/**
 * The mesh should be surrounded by at least one layer of ghost nodes.  These are nodes which
 * do not correspond to a cell, but are necessary for remeshing (because the remesher tries to
 * create a convex hull of the set of nodes) and visualising purposes.  The mesh is passed into
 * the constructor and the class is told about the ghost nodes by using this method.
 */
void SpheroidSimulation3D::SetGhostNodes(std::vector<unsigned> ghostNodeIndices)
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
void SpheroidSimulation3D::SetReMeshRule(bool remesh)
{
    mReMesh = remesh;
}

/**
 * Set the simulation to run with no birth.
 */
void SpheroidSimulation3D::SetNoBirth(bool nobirth)
{
    mNoBirth = nobirth;
}

/**
 * Set this simulation to use a cell killer
 */
void SpheroidSimulation3D::SetCellKiller(RandomCellKiller<3>* pCellKiller)
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
std::vector<MeinekeCryptCell> SpheroidSimulation3D::GetCells()
{
    assert(mCells.size()>0);
    return mCells;
}

/**
 * Whether each node is a ghost or not.
 * \todo change this to return a const reference
 */
std::vector <bool> SpheroidSimulation3D::GetGhostNodes()
{
    return mIsGhostNode;
}

/**
 * Return the index of each node on the left periodic boundary
 * \todo change this to return a const reference
 */
std::vector<unsigned> SpheroidSimulation3D::GetLeftCryptBoundary()
{
    return mLeftCryptBoundary;
}

/**
 * Return the index of each node on the right periodic boundary
 * \todo change this to return a const reference
 */
std::vector<unsigned> SpheroidSimulation3D::GetRightCryptBoundary()
{
    return mRightCryptBoundary;
}

/**
 * Return the index of each node on the whole boundary
 * \todo change this to return a const reference
 */
std::vector<unsigned> SpheroidSimulation3D::GetCryptBoundary()
{
    return mCryptBoundary;
}

/**
 * Get a node's location (ONLY FOR TESTING)
 *
 * @param the node index
 * @return the x and y co-ordinates of this node.
 */
std::vector<double> SpheroidSimulation3D::GetNodeLocation(const unsigned& rNodeIndex)
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
void SpheroidSimulation3D::Solve()
{ 
    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetDimensionalisedTime();
    std::cout << "Time at start of Solve Method = \n" << current_time << std::flush;
    
    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);
    std::cout << "num timesteps = " << num_time_steps << "\n" << std::flush;
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
    OutputFileHandler output_file_handler(results_directory+"/vis_results/",true);
    out_stream p_node_file = output_file_handler.OpenOutputFile("results.viznodes");
    out_stream p_element_file = output_file_handler.OpenOutputFile("results.vizelements");
    
    
    
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
            //Node<3> *p_our_node = mrMesh.GetNode(i);
         //   double y = p_our_node->rGetLocation()[1];
            std::vector<double> cell_cycle_influences;
//            if (mWntIncluded)
//            {
//                double wnt_stimulus = mWntGradient.GetWntLevel(y);
//                cell_cycle_influences.push_back(wnt_stimulus);
//            }
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
        CheckIndicesAreInSync();
        
        mRemeshesThisTimeStep = 0; // To avoid infinite loops
 //       std::cout << "** TIME = " << p_simulation_time->GetDimensionalisedTime() << " **" << std::endl;
        
        // Cell birth
        mNumBirths += DoCellBirth();
        
        //  calculate node velocities
        std::vector<std::vector<double> > drdt = CalculateVelocitiesOfEachNode();
        
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
        if (tabulated_output_counter > 80) // TODO: make this configurable!
        {
            tabulated_output_counter = 0;
        }
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
