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


#ifndef TESTHEARTCONFIG_HPP_
#define TESTHEARTCONFIG_HPP_

#include <cxxtest/TestSuite.h>
#include "HeartConfig.hpp"

class TestHeartConfig : public CxxTest::TestSuite
{
public :
    void TestHeartConfigBasic()
    {
        HeartConfig::Instance()->SetParametersFile("ChasteParameters.xml");
        
        double chi = HeartConfig::Instance()->Parameters()->Physiological().SurfaceAreaToVolumeRatio();
        TS_ASSERT_EQUALS(chi, 1400);
        
        double capacitance = HeartConfig::Instance()->Parameters()->Physiological().Capacitance();
        TS_ASSERT_EQUALS(capacitance, 1.0);

        double conductivity_1 = HeartConfig::Instance()->Parameters()->Physiological().IntracellularConductivities().longi();
        double conductivity_2 = HeartConfig::Instance()->Parameters()->Physiological().ExtracellularConductivities().longi();

        TS_ASSERT_EQUALS(conductivity_1, 1.75);
        TS_ASSERT_EQUALS(conductivity_2, 7.0);

        HeartConfig::Instance()->Destroy();
    }
};

#endif /*TESTHEARTCONFIG_HPP_*/
