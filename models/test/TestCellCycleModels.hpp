#ifndef TESTCELLCYCLEMODELS_HPP_
#define TESTCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"

class TestCellCycleModels : public CxxTest::TestSuite
{
public:
    void TestFixedCellCycleModel(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        FixedCellCycleModel our_fixed_cell;
        
        our_fixed_cell.SetCellType(TRANSIT);
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()-0.01));
        TS_ASSERT(our_fixed_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()));
        TS_ASSERT(our_fixed_cell.ReadyToDivide(p_params->GetTransitCellCycleTime()+0.01));
        
        our_fixed_cell.SetCellType(STEM);
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(p_params->GetStemCellCycleTime()-0.01));
        TS_ASSERT(our_fixed_cell.ReadyToDivide(p_params->GetStemCellCycleTime()));
        TS_ASSERT(our_fixed_cell.ReadyToDivide(p_params->GetStemCellCycleTime()+0.01));
        
        our_fixed_cell.SetCellType(DIFFERENTIATED);
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(1.0));
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(1e10));
        TS_ASSERT(!our_fixed_cell.ReadyToDivide(1e100));
    }
    
    void TestStochasticCellCycleModel(void) throw(Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerators rand_gen;
        StochasticCellCycleModel cell_model(&rand_gen);
        
        cell_model.SetCellType(STEM);
        TS_ASSERT(!cell_model.ReadyToDivide(p_params->GetStemCellCycleTime()-0.01));
        TS_ASSERT(cell_model.ReadyToDivide(p_params->GetStemCellCycleTime()));
        TS_ASSERT(cell_model.ReadyToDivide(p_params->GetStemCellCycleTime()+0.01));
        cell_model.SetCellType(DIFFERENTIATED);
        TS_ASSERT(!cell_model.ReadyToDivide(1.0));
        TS_ASSERT(!cell_model.ReadyToDivide(1e10));
        TS_ASSERT(!cell_model.ReadyToDivide(1e100));
        
        // Testing a random generator is hard...
        cell_model.SetCellType(TRANSIT);
        const int TESTS = 100;
        int ready_count = 0;
        for (int i=0; i<TESTS; i++)
        {
            if (cell_model.ReadyToDivide(p_params->GetTransitCellCycleTime()-0.1))
            {
                ready_count++;
            }
        }
        TS_ASSERT(ready_count>0);
    }
    
};

#endif /*TESTCELLCYCLEMODELS_HPP_*/
