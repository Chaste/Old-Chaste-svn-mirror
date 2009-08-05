/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef TESTCARDIACFASTSLOWPDE_HPP_
#define TESTCARDIACFASTSLOWPDE_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleStimulus.hpp"
#include "MonodomainFastSlowPde.hpp"
#include "BidomainFastSlowPde.hpp"
#include "FastSlowLuoRudyIModel1991.hpp"

class MyCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    boost::shared_ptr<SimpleStimulus> mpHalfStimulus;
public:
    MyCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-600.0, 0.5)),
          mpHalfStimulus(new SimpleStimulus(-300, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];

        if (fabs(x)<1e-6)
        {
            return new FastSlowLuoRudyIModel1991(mpSolver, mpZeroStimulus); // state unset at the moment
        }
        else if(fabs(x-0.5)<1e-6)
        {
            return new FastSlowLuoRudyIModel1991(mpSolver, mpHalfStimulus); // state unset at the moment
        }
        else
        {
            assert(fabs(x-1.0)<1e-6);
            return new FastSlowLuoRudyIModel1991(mpSolver, mpStimulus); // state unset at the moment
        }
    }

};



class TestCardiacFastSlowPde : public CxxTest::TestSuite
{
public:
    void TestMonodomainFastSlowPde() throw (Exception)
    {
        MixedTetrahedralMesh<2,2> mixed_mesh;
        mixed_mesh.ConstructRectangularMeshes(1.0,1.0,1,2);

        MyCellFactory cell_factory;
        cell_factory.SetMesh(mixed_mesh.GetFineMesh());

        MonodomainFastSlowPde<2> monodomain_fast_slow_pde( &cell_factory, mixed_mesh, 1.0 );

        /////////////////////////////////////////////////////////////////
        // check the cells have been correctly assigned as fast or slow
        /////////////////////////////////////////////////////////////////
        bool correct_cell_state[] = {false, true, false, true, true, true, false, true, false};

        std::vector<FastSlowLuoRudyIModel1991*> cells;
        for (unsigned cell_index = 0;
             cell_index < mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLocalOwnership();
             cell_index++)             
        {
            unsigned global_index = cell_index + mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow(); 
            
            cells.push_back(static_cast<FastSlowLuoRudyIModel1991*>(monodomain_fast_slow_pde.mCellsDistributed[cell_index]));
            TS_ASSERT_EQUALS( cells[cell_index]->IsFastOnly(), correct_cell_state[global_index] );
        }

        Vec voltage = PetscTools::CreateVec(9, -84.5);

        ///////////////////////////////////////////////////////////////////
        // check coarse-slow ODEs solved when they should be
        ///////////////////////////////////////////////////////////////////
        monodomain_fast_slow_pde.SolveCellSystems(voltage, 0, 0.5);
        TS_ASSERT_EQUALS(monodomain_fast_slow_pde.mLastSlowCurrentSolveTime, 0.0);
        TS_ASSERT_EQUALS(monodomain_fast_slow_pde.mNextSlowCurrentSolveTime, 1.0);

        monodomain_fast_slow_pde.SolveCellSystems(voltage, 0.5, 0.9);
        TS_ASSERT_EQUALS(monodomain_fast_slow_pde.mLastSlowCurrentSolveTime, 0.0);
        TS_ASSERT_EQUALS(monodomain_fast_slow_pde.mNextSlowCurrentSolveTime, 1.0);

        monodomain_fast_slow_pde.SolveCellSystems(voltage, 0.9, 1.5);
        TS_ASSERT_EQUALS(monodomain_fast_slow_pde.mLastSlowCurrentSolveTime, 1.0);
        TS_ASSERT_EQUALS(monodomain_fast_slow_pde.mNextSlowCurrentSolveTime, 2.0);

        monodomain_fast_slow_pde.SolveCellSystems(voltage, 1.2, 2.0);
        // don't test whether mNextTime = 3.0, as it may or may not, depending on floating point issues..

        monodomain_fast_slow_pde.SolveCellSystems(voltage, 2.0, 2.5);
        TS_ASSERT_EQUALS(monodomain_fast_slow_pde.mLastSlowCurrentSolveTime, 2.0);
        TS_ASSERT_EQUALS(monodomain_fast_slow_pde.mNextSlowCurrentSolveTime, 3.0);


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // test the values at the end.. (taken for TestMixedOdes - without testing voltage as ComputeExcVolt is used
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ReplicatableVector cell_voltage(9);
        ReplicatableVector gating_var(9);
        ReplicatableVector slow_value_0(9);
        ReplicatableVector slow_value_1(9);


        for (unsigned cell_index = 0;
             cell_index < mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLocalOwnership();
             cell_index++)             
        {
            unsigned global_index = cell_index + mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow();             
            
            cell_voltage[global_index] = cells[cell_index]->GetVoltage();
            gating_var[global_index] = cells[cell_index]->rGetStateVariables()[0];

            std::vector<double> slow_values(2);
            if (!cells[cell_index]->IsFastOnly())
            {
                cells[cell_index]->GetSlowValues(slow_values);
                slow_value_0[global_index] = slow_values[0];
                slow_value_1[global_index] = slow_values[1];
            }
        }

        cell_voltage.Replicate(mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow(), mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetHigh());
        gating_var.Replicate(mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow(), mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetHigh());
        slow_value_0.Replicate(mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow(), mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetHigh());
        slow_value_1.Replicate(mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow(), mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetHigh());

        // test cells for which x=0: correspond to fine nodes 0 (coarse-slow cell), 3 (fine-fast cell), 6 (coarse-slow cell)
        TS_ASSERT_DELTA( cell_voltage[0], cell_voltage[3], 1.0); // check coarse and fine agree in voltage
        TS_ASSERT_DELTA( gating_var[0], gating_var[3], 0.01); // check coarse and fine agree in a gating var
        TS_ASSERT_DELTA(slow_value_0[0], slow_value_0[3], 0.01); // check slow values match
        TS_ASSERT_DELTA(slow_value_1[0], slow_value_1[3], 0.01); // check slow values match

        // test cells for which x=0.5: correspond to fine nodes 1, 4, 7 (all fine-fast cells)
        TS_ASSERT_DELTA( cell_voltage[1], cell_voltage[4], 1.0);

        // test cells for which x=1: correspond to fine nodes 2 (coarse-slow), 5 (fine-fast), 8 (coarse-slow)
        TS_ASSERT_DELTA( cell_voltage[2], cell_voltage[5], 1.0); // check coarse and fine agree in voltage
        TS_ASSERT_DELTA( gating_var[2], gating_var[5], 0.01); // check coarse and fine agree in a gating var

        TS_ASSERT_LESS_THAN(0, slow_value_0[2]);
        TS_ASSERT_DELTA(slow_value_0[2], slow_value_0[5], 0.01); // check slow values match
        TS_ASSERT_DELTA(slow_value_1[2], slow_value_1[5], 0.01); // check slow values match

        VecDestroy(voltage);
    }


    // similar test to TestMonodomainFastSlowPde, but smaller, as testing the same functionality.
    void TestBidomainFastSlowPde() throw (Exception)
    {
        MixedTetrahedralMesh<2,2> mixed_mesh;
        mixed_mesh.ConstructRectangularMeshes(1.0,1.0,1,2);

        MyCellFactory cell_factory;
        cell_factory.SetMesh(mixed_mesh.GetFineMesh());

        BidomainFastSlowPde<2> bidomain_fast_slow_pde( &cell_factory, mixed_mesh, 1.0 );

        /////////////////////////////////////////////////////////////////
        // check the cells have been correctly assigned as fast or slow
        /////////////////////////////////////////////////////////////////
        bool correct_cell_state[] = {false, true, false, true, true, true, false, true, false};

        std::vector<FastSlowLuoRudyIModel1991*> cells;
        for (unsigned cell_index = 0;
             cell_index < mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLocalOwnership();
             cell_index++)             
        {
            unsigned global_index = cell_index + mixed_mesh.GetFineMesh()->GetDistributedVectorFactory()->GetLow();             
            
            cells.push_back(static_cast<FastSlowLuoRudyIModel1991*>(bidomain_fast_slow_pde.mCellsDistributed[cell_index]));
            TS_ASSERT_EQUALS( cells[cell_index]->IsFastOnly(), correct_cell_state[global_index] );
        }

        Vec voltage = PetscTools::CreateVec(9, -84.5);

        ///////////////////////////////////////////////////////////////////
        // check coarse-slow ODEs solved when they should be
        ///////////////////////////////////////////////////////////////////
        bidomain_fast_slow_pde.SolveCellSystems(voltage, 0, 0.5);
        TS_ASSERT_EQUALS(bidomain_fast_slow_pde.mLastSlowCurrentSolveTime, 0.0);
        TS_ASSERT_EQUALS(bidomain_fast_slow_pde.mNextSlowCurrentSolveTime, 1.0);

        VecDestroy(voltage);
    }
};

#endif /*TESTCARDIACFASTSLOWPDE_HPP_*/


