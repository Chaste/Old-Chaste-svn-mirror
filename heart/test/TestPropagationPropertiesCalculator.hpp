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


#ifndef _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>

#include "PropagationPropertiesCalculator.hpp"


class TestPropagationPropertiesCalculator : public CxxTest::TestSuite
{
public:

    void TestConductionVelocity1D(void) throw (Exception)
    {
        Hdf5DataReader simulation_data("heart/test/data/Monodomain1d",
                                       "MonodomainLR91_1d", false);

        PropagationPropertiesCalculator ppc(&simulation_data);

        double velocity = ppc.CalculateConductionVelocity(5,15,0.1);
        TS_ASSERT_DELTA(velocity, 0.0499, 0.001);

        // Should throw because AP does not propagate far enough in simulation time
        TS_ASSERT_THROWS_ANYTHING(ppc.CalculateConductionVelocity(5,95,0.9));

        // Should throw because AP does not propagate far enough in simulation time
        TS_ASSERT_THROWS_ANYTHING(ppc.CalculateConductionVelocity(90,100,0.1));

        TS_ASSERT_DELTA(ppc.CalculateMaximumUpstrokeVelocity(1),343.9429,0.001);

        TS_ASSERT_DELTA(ppc.CalculatePeakMembranePotential(5),23.6271,0.001);

        TS_ASSERT_DELTA(ppc.CalculateActionPotentialDuration(50,5),0,0.001);
    }
    
    void TestConductionBidomain3D() throw (Exception)
    {
        //Note: these data files (from notforrelease/test/TestCardiacFastSlowProblem3D.hpp) are incomplete.
        unsigned middle_index = 14895U;
        unsigned rhs_index = 14910U;
   
        
        Hdf5DataReader simulation_data_fs("heart/test/data/BidomainFastSlow3D",
                                       "res", false);
        PropagationPropertiesCalculator properties_fs(&simulation_data_fs);

        /*
         * From Raf's Matlab/Octave script
        r4256_Bi_FastSlow_3d
        Mid APD90: 229.7
        Mid dv/dt_max: 180.28
        Mid dv/dt time: 3.2
        RHS dv/dt time: 5.8
        CV: 0.057692
        *
        */
        TS_ASSERT_DELTA(properties_fs.CalculateActionPotentialDuration(90, middle_index), 229.7, 0.25);
        TS_ASSERT_DELTA(properties_fs.CalculateConductionVelocity(middle_index, rhs_index, 0.15), 0.057692, 0.001);
        TS_ASSERT_DELTA(properties_fs.CalculateMaximumUpstrokeVelocity(middle_index), 180.28, 0.01);
    

        Hdf5DataReader simulation_data_bw("heart/test/data/BidomainBackwardToCompareWithFastSlow3D",
                                       "res", false);
        PropagationPropertiesCalculator properties_bw(&simulation_data_bw);
        /*
         * From Raf's Matlab/Octave script
        r4256_Bi_Back_3d
        Mid APD90: 229.2
        Mid dv/dt_max: 173.02
        Mid dv/dt time: 3.3
        RHS dv/dt time: 6
        CV: 0.055556
        *
        */        
        TS_ASSERT_DELTA(properties_bw.CalculateActionPotentialDuration(90, middle_index), 229.2, 0.25);
        TS_ASSERT_DELTA(properties_bw.CalculateConductionVelocity(middle_index, rhs_index, 0.15), 0.055556, 0.001);
        TS_ASSERT_DELTA(properties_bw.CalculateMaximumUpstrokeVelocity(middle_index), 173.02, 0.01);
        
        
        
    }
    
    void TestUpstrokeBidomain3D()
    {
        Hdf5DataReader mono_fs_reader("heart/test/data/MonodomainFastSlow3D", "res", false);
        //std::vector<double> voltage_fast_slow=mono_fs_reader.GetVariableOverTime("V", mMiddleIndex);
        //DumpVector(voltage_fast_slow, "middle_node_mfs.dat"); 
        
        PropagationPropertiesCalculator ppc_fs(&mono_fs_reader);
        
        TS_ASSERT_DELTA(ppc_fs.CalculateMaximumUpstrokeVelocity(14895U), 174.9, 4.0);
        TS_ASSERT_DELTA(ppc_fs.CalculateActionPotentialDuration(90, 14895U), 230.25, 1.0);
        
        Hdf5DataReader mono_bw_reader("heart/test/data/MonodomainBackwardToCompareWithFastSlow3D", "res", false);
        //std::vector<double> voltage_fast_slow=mono_fs_reader.GetVariableOverTime("V", mMiddleIndex);
        //DumpVector(voltage_fast_slow, "middle_node_mfs.dat"); 
        
        PropagationPropertiesCalculator ppc_bw(&mono_bw_reader);
        
        TS_ASSERT_DELTA(ppc_bw.CalculateMaximumUpstrokeVelocity(14895U), 174.9, 4.0);
        TS_ASSERT_DELTA(ppc_bw.CalculateActionPotentialDuration(90, 14895U), 230.25, 1.0);
     }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
