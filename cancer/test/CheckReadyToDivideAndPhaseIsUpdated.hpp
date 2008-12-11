/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_
#define CHECKREADYTODIVIDEANDPHASEISUPDATED_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>

#include "AbstractCellCycleModel.hpp"

void CheckReadyToDivideAndPhaseIsUpdated(AbstractCellCycleModel* pModel,
                                         double g1Duration,
                                         double g2Duration=CancerParameters::Instance()->GetG2Duration())
{
    double age = pModel->GetAge();
    CancerParameters* p_params = CancerParameters::Instance();
    
    const double G1TOL = 1e-5; // How accurate the expected G1 duration is
    
    if ((pModel->GetCell()->GetCellType() != DIFFERENTIATED) &&
        (age >= p_params->GetMDuration()) &&
        (pModel->GetG1Duration() != DOUBLE_UNSET) &&
        (fabs(pModel->GetG1Duration() - g1Duration) > G1TOL))
    {
        std::cout << "G1 duration mismatch: actual = " << pModel->GetG1Duration()
                  << ", expected = " << g1Duration << std::endl;
    }

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
    else if (age < p_params->GetMDuration() + g1Duration - G1TOL)
    {   // if in G1 phase
        TS_ASSERT(!pModel->ReadyToDivide());
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(),G_ONE_PHASE);
        if (pModel->GetCurrentCellCyclePhase() != G_ONE_PHASE)
        {
            std::cout << "Expected G1: " << g1Duration << "; actual: " << pModel->GetG1Duration()
                      << "; age = " << age << "; G1-S transition = " << p_params->GetMDuration() + g1Duration << std::endl;
        }
    }
    else if (age < p_params->GetMDuration() + g1Duration + p_params->GetSDuration() - G1TOL)
    {   // if in S phase
        TS_ASSERT(!pModel->ReadyToDivide());
        TS_ASSERT_EQUALS(pModel->GetCurrentCellCyclePhase(),S_PHASE);
    }
    else if (age < p_params->GetMDuration() + g1Duration + p_params->GetSDuration() + g2Duration  - G1TOL)
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
