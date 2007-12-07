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
            
    if (cell_type != STEM)
    {
        //std::cout << "Generation  " << mGeneration << " max trans " << p_params->GetMaxTransitGenerations() << "\n";
        if (cell_type == HEPA_ONE)
        {
            new_cell_types[0] = cell_type;
            new_cell_types[1] = cell_type;
        }
        else if ((mGeneration-1u) < p_params->GetMaxTransitGenerations())
        {
            //std::cout << "Generation T " << mGeneration << "\n";
            new_cell_types[0] = cell_type;
            new_cell_types[1] = TRANSIT;
        }
        else
        {
            //std::cout << "Generation D " << mGeneration << "\n";
            new_cell_types[0] = DIFFERENTIATED;
            new_cell_types[1] = DIFFERENTIATED;
        }
    }
    else
    {
        new_cell_types[0] = STEM;
        new_cell_types[1] = TRANSIT;
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
