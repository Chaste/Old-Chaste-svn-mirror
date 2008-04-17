#ifndef CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_
#define CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>

#include "AbstractCellCycleModel.hpp"

void CheckReadyToDivideAndPhaseIsUpdated(AbstractCellCycleModel* pModel, double g1Duration)
{   
    double age = pModel->GetAge();
    CancerParameters* p_params = CancerParameters::Instance();
    
    if (pModel->GetCell()->GetCellType()==DIFFERENTIATED)
    {
        TS_ASSERT(!pModel->ReadyToDivide());
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(),G_ZERO_PHASE);  
    }
    else if (age < p_params->GetMDuration())
    {   // if in M phase
        TS_ASSERT(!pModel->ReadyToDivide());
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(),M_PHASE);
    }
    else if (age < p_params->GetMDuration() + g1Duration)
    {   // if in G1 phase
        TS_ASSERT(!pModel->ReadyToDivide());
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(),G_ONE_PHASE);
    }
    else if (age < p_params->GetMDuration() + g1Duration + p_params->GetSDuration())
    {   // if in S phase
        TS_ASSERT(!pModel->ReadyToDivide());
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(),S_PHASE);
    }
    else if (age < p_params->GetMDuration() + g1Duration + p_params->GetSDuration() + p_params->GetG2Duration() )
    {   // if in G2 phase
        TS_ASSERT(!pModel->ReadyToDivide());
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(),G_TWO_PHASE);
    }
    else
    {
        TS_ASSERT(pModel->ReadyToDivide());
    }
}

#endif /*CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_*/
