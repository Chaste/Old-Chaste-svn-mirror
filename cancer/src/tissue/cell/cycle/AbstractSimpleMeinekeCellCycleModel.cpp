#include "AbstractSimpleMeinekeCellCycleModel.hpp"

void AbstractSimpleMeinekeCellCycleModel::ResetModel()
{ 
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    SetG1Duration();
}

std::vector<CellType> AbstractSimpleMeinekeCellCycleModel::GetNewCellTypes(CellType cellType)
{
   
    std::vector<CellType> cell_types(2);
    CancerParameters *p_params = CancerParameters::Instance();
    
    if (cellType != STEM)
    {
        //std::cout << "Generation  " << mGeneration << " max trans " << p_params->GetMaxTransitGenerations() << "\n";
        if (cellType == HEPA_ONE)
        {
            assert(GetGeneration()==mGeneration);  
            cell_types[0] = cellType;
            cell_types[1] = cellType;
        }
        else if ((mGeneration-1u) < p_params->GetMaxTransitGenerations())
        {
            //std::cout << "Generation T " << mGeneration << "\n";
            assert(GetGeneration()==mGeneration);     
            cell_types[0] = cellType;
            cell_types[1] = TRANSIT;
        }
        else
        {
            //std::cout << "Generation D " << mGeneration << "\n";
            assert(GetGeneration()==mGeneration);
            cell_types[0] = DIFFERENTIATED;
            cell_types[1] = DIFFERENTIATED;
        }
    }
    else
    {
        //SetGeneration(0u);                      
        cell_types[0] = cellType;
        cell_types[1] = TRANSIT;
        
    }
    
    return cell_types;
}

void AbstractSimpleMeinekeCellCycleModel::SetMotherGeneration(CellType cellType)
{
    if (cellType == STEM)
    {
        mGeneration--;
    }
    
}
