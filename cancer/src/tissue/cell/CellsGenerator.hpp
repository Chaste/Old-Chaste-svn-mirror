#ifndef CELLSGENERATOR_HPP_
#define CELLSGENERATOR_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "IngeWntSwatCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"

// Possible types of Cell Cycle Model (just for GenerateForCrypt method)
typedef enum CellCycleType_
{
    FIXED,
    STOCHASTIC,
    SIMPLE_WNT,
    WNT,
    INGE_WNT_SWAT_HYPOTHESIS_ONE,
    INGE_WNT_SWAT_HYPOTHESIS_TWO,
    STOCHASTIC_WNT,
    TYSONNOVAK
} CellCycleType;


/**
 * A helper class for generating a vector of cells for a given mesh
 */
template<unsigned DIM>
class CellsGenerator
{
public :
    /**
     * Fills a vector of cells with a Fixed cell cycle model, to match 
     * a given mesh. Gives them birth times of 0 for node 0,
     * -1 for node 1, -2 for node 2 etc...
     * 
     * @param rCells  An empty vector of cells to fill up.
     * @param rMesh  The mesh the cells should be associated with.  
     */
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

    /**
     * Generates cells of a specified cell cycle type under the correct 
     * crypt conditions and gives random ages if required, 
     * or gives them an age of 0.0 - creates least work for solver startup.
     * 
     * @param rCells  An empty cells vector for this function to fill up
     * @param rMesh  The crypt mesh (can be cylindrical)
     * @param cycleType (As specified in the enumeration at the top of this file)
     * @param randomBirthTimes  Whether to assign the cells random birth times (this can be expensive computationally with ODE models)
     * @param y0  below this line cells are generation 0 (defaults to 0.3)
     * @param y1  below this line cells are generation 1 (defaults to 2.0)
     * @param y2  below this line cells are generation 2 (defaults to 3.0)
     * @param y3  below this line cells are generation 3 (defaults to 4.0)
     * 
     * \todo Only give generation information to relevant models (#509)
     */
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
        
        rCells.clear();
        rCells.reserve(num_cells);
        
        if(cycleType!=FIXED && cycleType!=STOCHASTIC)
        {   // Only these two models use a fixed number of transit generations.
            CancerParameters::Instance()->SetMaxTransitGenerations(UINT_MAX);                  
        }
        
        for (unsigned i=0; i<num_cells; i++)
        {
            CellType cell_type;
            unsigned generation;

            double y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (cycleType==FIXED)
            {
                p_cell_cycle_model = new FixedCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellG1Duration()
                                            + p_params->GetSG2MDuration();
                typical_stem_cycle_time = p_params->GetStemCellG1Duration()
                                            + p_params->GetSG2MDuration();
            }
            else if (cycleType==STOCHASTIC)
            {
                p_cell_cycle_model = new StochasticCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellG1Duration()
                                                + p_params->GetSG2MDuration();
                typical_stem_cycle_time = p_params->GetStemCellG1Duration()
                                            + p_params->GetSG2MDuration();
            }
            else if (cycleType==SIMPLE_WNT)
            {
                p_cell_cycle_model = new SimpleWntCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellG1Duration()
                                               + p_params->GetSG2MDuration();
                typical_stem_cycle_time = p_params->GetStemCellG1Duration()
                                            + p_params->GetSG2MDuration();
            }
            else if (cycleType==WNT)
            {
                p_cell_cycle_model = new WntCellCycleModel();
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==INGE_WNT_SWAT_HYPOTHESIS_ONE)
            {
                p_cell_cycle_model = new IngeWntSwatCellCycleModel(1u);
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==INGE_WNT_SWAT_HYPOTHESIS_TWO)
            {
                p_cell_cycle_model = new IngeWntSwatCellCycleModel(2u);
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
            { // Useful to provide a warning if any new cell cycle models are being used.
                #define COVERAGE_IGNORE
                EXCEPTION("Cell Cycle Type is not recognised");
                #undef COVERAGE_IGNORE   
            }

            double birth_time = 0.0;
            
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
                if(cycleType!=FIXED && cycleType!=STOCHASTIC)
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
