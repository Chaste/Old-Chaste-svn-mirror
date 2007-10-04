#ifndef CELLSGENERATOR_HPP_
#define CELLSGENERATOR_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "WntGradient.hpp"



// Possible types of Cell Cycle Model (just for GenerateForCrypt method)
typedef enum CellCycleType_
{
    FIXED,
    STOCHASTIC,
    WNT,
    STOCHASTIC_WNT,
    TYSONNOVAK
} CellCycleType;


/*
 * A helper class for generating a vector of cells for a given mesh
 */
template<unsigned DIM>
class CellsGenerator
{
public :
    static void GenerateBasic(std::vector<TissueCell>& rCells, 
                               ConformingTetrahedralMesh<DIM,DIM>& rMesh)
    {
        rCells.clear();
        rCells.reserve(rMesh.GetNumNodes());
        for(unsigned i=0; i<rMesh.GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            rCells.push_back(cell);
        }
    }

    static void GenerateForCrypt(std::vector<TissueCell>& rCells, 
                                 ConformingTetrahedralMesh<2,2>& rMesh, 
                                 CellCycleType cycleType,
                                 bool randomBirthTimes,
                                 double y0 = 0.3,
                                 double y1 = 2.0,
                                 double y2 = 3.0,
                                 double y3 = 4.0)
    {
        assert(DIM==2);
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        unsigned num_cells = rMesh.GetNumNodes();

        AbstractCellCycleModel* p_cell_cycle_model = NULL;
        double typical_transit_cycle_time;
        double typical_stem_cycle_time;
        
        CancerParameters* p_params = CancerParameters::Instance();
        
        rCells.reserve(num_cells);
        
        for (unsigned i=0; i<num_cells; i++)
        {
            CellType cell_type;
            unsigned generation;

            double y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (cycleType==FIXED)
            {
                p_cell_cycle_model = new FixedCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==STOCHASTIC)
            {
                p_cell_cycle_model = new StochasticCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==WNT)
            {
                p_cell_cycle_model = new WntCellCycleModel();
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==STOCHASTIC_WNT)
            {
                p_cell_cycle_model = new StochasticWntCellCycleModel();
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==TYSONNOVAK)
            {
                p_cell_cycle_model = new TysonNovakCellCycleModel();
                typical_transit_cycle_time = 1.25;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else
            {
                EXCEPTION("Cell Cycle Type is not recognised");   
            }
            

            double birth_time = -2; 
            
            if (y <= y0)
            {
                cell_type = STEM;
                generation = 0;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_stem_cycle_time; // hours
                }
            }
            else if (y < y1)
            {
                cell_type = TRANSIT;
                generation = 1;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y2)
            {
                cell_type = TRANSIT;
                generation = 2;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y3)
            {
                cell_type = TRANSIT;
                generation = 3;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else
            {
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
                if(cycleType==WNT || cycleType==TYSONNOVAK)
                {
                    // There are no fully differentiated cells!
                    cell_type = TRANSIT;
                    
                }
                else
                {
                    cell_type = DIFFERENTIATED;
                }                
                generation = 4;
            }

            TissueCell cell(cell_type, HEALTHY, generation, p_cell_cycle_model);
            
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            rCells.push_back(cell);
        }
    }
};

#endif /*CELLSGENERATOR_HPP_*/
