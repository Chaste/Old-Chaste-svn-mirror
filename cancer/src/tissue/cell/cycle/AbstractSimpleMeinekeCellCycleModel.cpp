#include "AbstractSimpleMeinekeCellCycleModel.hpp"

void AbstractSimpleMeinekeCellCycleModel::ResetModel()
{ 
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    SetG1Duration();
}

std::vector<CellType> AbstractSimpleMeinekeCellCycleModel::GetNewCellTypes()
{
    CellType cell_type = mpCell->GetCellType();
    std::vector<CellType> new_cell_types(2);
    CancerParameters *p_params = CancerParameters::Instance();
    
    if (cell_type == STEM)
    {
        new_cell_types[0] = STEM;
        new_cell_types[1] = TRANSIT;
    }       
    else
    {
        if ((mGeneration-1u) < p_params->GetMaxTransitGenerations())
        {
            new_cell_types[0] = cell_type;
            new_cell_types[1] = TRANSIT;
        }
        else
        {
            new_cell_types[0] = DIFFERENTIATED;
            new_cell_types[1] = DIFFERENTIATED;
        }        
    }
               
    return new_cell_types;
}

void AbstractSimpleMeinekeCellCycleModel::SetMotherGeneration()
{
    CellType cell_type = mpCell->GetCellType();
    if (cell_type == STEM)
    {
        mGeneration = 0u;
    }
}
