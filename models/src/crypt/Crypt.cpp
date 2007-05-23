#ifndef CRYPT_CPP
#define CRYPT_CPP

#include "Crypt.hpp"


///\todo: make this constructor take in ghost nodes, and validate the three objects
// are in sync ie num cells + num ghost nodes = num_nodes ? this would mean all ghosts
// *cannot* be cells, making it more difficult to construct the cells.
// also check cell.GetNodeIndices() is in the mesh, and covers the mesh, etc.
template<unsigned DIM>
Crypt<DIM>::Crypt(ConformingTetrahedralMesh<DIM, DIM>& rMesh,
                  std::vector<MeinekeCryptCell>& rCells)
             : mrMesh(rMesh),
               mrCells(rCells)
{
    mSelfSetGhostNodes=true;
    mpGhostNodes=new std::vector<bool>(mrMesh.GetNumNodes(), false);
    
    mMaxCells = 10*mrMesh.GetNumNodes();
    mMaxElements = 10*mrMesh.GetNumElements();
}

template<unsigned DIM>
Crypt<DIM>::~Crypt()
{
    if (mSelfSetGhostNodes)
    {
        delete mpGhostNodes;
    }
}

template<unsigned DIM>
ConformingTetrahedralMesh<DIM, DIM>& Crypt<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
std::vector<MeinekeCryptCell>& Crypt<DIM>::rGetCells()
{
    return mrCells;
}

template<unsigned DIM>
std::vector<bool>& Crypt<DIM>::rGetGhostNodes()
{
    return *mpGhostNodes;
}

template<unsigned DIM>
void Crypt<DIM>::SetGhostNodes(std::vector<bool>& rGhostNodes)
{
    if (mSelfSetGhostNodes)
    {
        delete mpGhostNodes;
    }
    mSelfSetGhostNodes=false;
    mpGhostNodes = &rGhostNodes;
}

template<unsigned DIM>
void Crypt<DIM>::RemoveDeadCells()
{
    std::vector< MeinekeCryptCell > living_cells;
    for (unsigned i=0; i<mrCells.size(); i++)
    {
        MeinekeCryptCell* p_cell=&(mrCells[i]);
        if (p_cell->IsDead())
        {
            mrMesh.DeleteNode(p_cell->GetNodeIndex());
        }
        else
        {
            living_cells.push_back(*p_cell);
        }
    }
    
    mrCells=living_cells;
    //Remesh and re-index (is moved to caller)     
}

template<unsigned DIM>
MeinekeCryptCell& Crypt<DIM>::Iterator::operator*()
{
    assert((*this) != mrCrypt.End());
    return mrCrypt.rGetCells()[mCellIndex];
}


template<unsigned DIM>
MeinekeCryptCell* Crypt<DIM>::Iterator::operator->()
{
    assert((*this) != mrCrypt.End());
    return &(mrCrypt.rGetCells()[mCellIndex]);
}


template<unsigned DIM>
Node<DIM>* Crypt<DIM>::Iterator::GetNode()
{
    assert((*this) != mrCrypt.End());
    return mrCrypt.rGetMesh().GetNode(mNodeIndex);
}

template<unsigned DIM>
const c_vector<double, DIM>& Crypt<DIM>::Iterator::rGetLocation()
{
    return GetNode()->rGetLocation();
}

template<unsigned DIM>
bool Crypt<DIM>::Iterator::operator!=(const Crypt<DIM>::Iterator& other)
{
    return mCellIndex != other.mCellIndex;   
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator& Crypt<DIM>::Iterator::operator++()
{
    do
    {
        mCellIndex++;
        if((*this) != mrCrypt.End())
        {
            mNodeIndex = mrCrypt.rGetCells()[mCellIndex].GetNodeIndex();
        }
    }
    while ((*this) != mrCrypt.End() && !IsRealCell());
    
      
    return (*this);
}

template<unsigned DIM>
bool Crypt<DIM>::Iterator::IsRealCell()
{
    assert(mrCrypt.rGetGhostNodes().size() == mrCrypt.rGetMesh().GetNumAllNodes() );
    return !(mrCrypt.rGetGhostNodes()[mNodeIndex] || GetNode()->IsDeleted());
}

template<unsigned DIM>
Crypt<DIM>::Iterator::Iterator(Crypt& rCrypt, unsigned cellIndex)
    : mrCrypt(rCrypt),
      mCellIndex(cellIndex)
{
    mNodeIndex=0;//\todo
    // Make sure the crypt isn't empty
    assert(mrCrypt.rGetCells().size() > 0);
    // Make sure we start at a real cell
    if (mCellIndex == 0 && !IsRealCell())
    {
        ++(*this);
    }
    if (mCellIndex != mrCrypt.rGetCells().size())
    {
        mNodeIndex = mrCrypt.rGetCells()[mCellIndex].GetNodeIndex();
    }
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator Crypt<DIM>::Begin()
{
    return Iterator(*this, 0);
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator Crypt<DIM>::End()
{
    return Iterator(*this, mrCells.size());
}

template<unsigned DIM>
void Crypt<DIM>::UpdateGhostPositions(const std::vector< c_vector<double, DIM> >& rDrDt, double dt)
{
    for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
    {
        if ((!mrMesh.GetNode(index)->IsDeleted()) && (*mpGhostNodes)[index])
        {
            Point<DIM> new_point(mrMesh.GetNode(index)->rGetLocation() + dt*rDrDt[index]);
            mrMesh.SetNode(index, new_point, false);
        }
    }
}

template<unsigned DIM>
void Crypt<DIM>::MoveCell(Iterator iter, Point<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    mrMesh.SetNode(index, rNewLocation, false);
}

template<unsigned DIM>  
void Crypt<DIM>::AddCell(MeinekeCryptCell newCell, c_vector<double,DIM> newLocation)
{
    Node<DIM>* p_new_node = new Node<DIM>(mrMesh.GetNumNodes(), newLocation, false);   // never on boundary
                
    NodeMap map(mrMesh.GetNumNodes());
    unsigned new_node_index = mrMesh.AddNodeAndReMesh(p_new_node,map);

    newCell.SetNodeIndex(new_node_index);
    if (new_node_index == mrCells.size())
    {
        mrCells.push_back(newCell);
    }
    else
    {
        #define COVERAGE_IGNORE
        mrCells[new_node_index] = newCell;
        #undef COVERAGE_IGNORE
    }

    // Update size of IsGhostNode if necessary
    if (mrMesh.GetNumNodes() > mpGhostNodes->size())
    {
        mpGhostNodes->resize(mrMesh.GetNumNodes());
        (*mpGhostNodes)[new_node_index] = false;
    }   
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//                             output methods                               // 
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>  
void Crypt<DIM>::SetMaxCells(unsigned maxCells)
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
void Crypt<DIM>::SetMaxElements(unsigned maxElements)
{
    mMaxElements = maxElements;
    if (maxElements<mrMesh.GetNumAllElements())
    {
        #define COVERAGE_IGNORE
        EXCEPTION("mMaxElements is less than the number of elements in the mesh.");
        #undef COVERAGE_IGNORE
    }
}

template<unsigned DIM>  
void Crypt<DIM>::SetupTabulatedWriters(ColumnDataWriter& rNodeWriter, ColumnDataWriter& rElementWriter)
{   
    ////////////////////////
    // set up node writer
    ////////////////////////
    mNodeVarIds.time = rNodeWriter.DefineUnlimitedDimension("Time","hours");
    
    mNodeVarIds.types.resize(mMaxCells);
    mNodeVarIds.position_id.resize(mMaxCells);
    
    // set up per-cell variables
    for (unsigned cell=0; cell<mMaxCells; cell++)
    {
        std::stringstream cell_type_var_name, cell_x_position_var_name, cell_y_position_var_name, cell_z_position_var_name;
        cell_type_var_name << "cell_type_" << cell;
        
        cell_x_position_var_name << "cell_x_position_" << cell;
        cell_y_position_var_name << "cell_y_position_" << cell;
        cell_z_position_var_name << "cell_z_position_" << cell; // not used in 2d

        std::vector<std::string> cell_position_var_name_string;
        cell_position_var_name_string.push_back(cell_x_position_var_name.str());
        cell_position_var_name_string.push_back(cell_y_position_var_name.str());
        cell_position_var_name_string.push_back(cell_z_position_var_name.str()); // not used in 2d
        
        mNodeVarIds.types[cell]=rNodeWriter.DefineVariable(cell_type_var_name.str(),"dimensionless");
        for(unsigned i=0; i<DIM; i++)
        {
            mNodeVarIds.position_id[cell](i)=rNodeWriter.DefineVariable(cell_position_var_name_string[i],"rest_spring_length");
        }
    }
    
    rNodeWriter.EndDefineMode();

    //////////////////////////
    // set up element writer
    //////////////////////////
    mElemVarIds.time = rElementWriter.DefineUnlimitedDimension("Time","hours");
    
    // Set up columns for element writer
    mElemVarIds.node_id.resize(mMaxElements);
    
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
            mElemVarIds.node_id[elem_index](i) = rElementWriter.DefineVariable(node_var_name_string[i],"dimensionless");
        }
    }
    
    rElementWriter.EndDefineMode();
}


template<unsigned DIM>  
void Crypt<DIM>::WriteResultsToFiles(ColumnDataWriter& rNodeWriter,
                                     ColumnDataWriter& rElementWriter,
                                     std::ofstream& rNodeFile, 
                                     std::ofstream& rElementFile,
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
        rNodeWriter.PutVariable(mNodeVarIds.time, time);
        rElementWriter.PutVariable(mElemVarIds.time, time);
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
        
        if ((*mpGhostNodes)[index]==true)
        {
            colour = 4; // visualizer treats '4' these as invisible
        }
// remove this else - facade eventually shouldn't be able to have empty cells vector
        else if (mrCells.size()>0)
        {
            if (index < mrCells.size())
            {
                CryptCellType type = mrCells[index].GetCellType();
                CryptCellMutationState mutation = mrCells[index].GetMutationState();
                
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
            const c_vector<double,DIM>& position = mrMesh.GetNode(index)->rGetLocation();
            if (writeVisualizerResults)
            {
                for(unsigned i=0; i<DIM; i++)
                {
                    rNodeFile << position[i] << " ";
                }
                rNodeFile << colour << " ";
            }
            if (writeTabulatedResults)
            {
                for(unsigned i=0; i<DIM; i++)
                {
                    rNodeWriter.PutVariable(mNodeVarIds.position_id[index](i), position[i]);
                }
                rNodeWriter.PutVariable(mNodeVarIds.types[index], colour);
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
                for(unsigned i=0; i<DIM+1; i++)
                {
                    rElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i)<< " ";
                }
            }
            if (writeTabulatedResults)
            {
                for(unsigned i=0; i<DIM+1; i++)
                {
                    rElementWriter.PutVariable(mElemVarIds.node_id[elem_index](i), mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i));
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



#endif //CRYPT_CPP

