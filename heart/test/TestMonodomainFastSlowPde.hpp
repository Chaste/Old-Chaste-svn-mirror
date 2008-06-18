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


#ifndef TESTMONODOMAINFASTSLOWPDE_HPP_
#define TESTMONODOMAINFASTSLOWPDE_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleStimulus.hpp"
#include "MonodomainFastSlowPde.hpp"
#include "FastSlowLuoRudyIModel1991.hpp"

class MyCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    SimpleStimulus* mpStimulus;
    SimpleStimulus* mpHalfStimulus;
public:
    MyCellFactory() : AbstractCardiacCellFactory<2>(0.01)
    {
        mpStimulus = new SimpleStimulus(-600.0, 0.5);
        mpHalfStimulus = new SimpleStimulus(-300, 0.5);
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned nodeIndex)
    {
        double x = mpMesh->GetNode(nodeIndex)->rGetLocation()[0];

        if (fabs(x)<1e-6)
        {
            return new FastSlowLuoRudyIModel1991(mpSolver, mTimeStep, mpZeroStimulus); // state unset at the moment
        }
        else if(fabs(x-0.5)<1e-6)
        {
            return new FastSlowLuoRudyIModel1991(mpSolver, mTimeStep, mpHalfStimulus); // state unset at the moment
        }
        else
        {
            assert(fabs(x-1.0)<1e-6);
            return new FastSlowLuoRudyIModel1991(mpSolver, mTimeStep, mpStimulus); // state unset at the moment
        }
    }

    ~MyCellFactory(void)
    {
        delete mpStimulus;
        delete mpHalfStimulus;
    }
};



class TestMonodomainFastSlowPde : public CxxTest::TestSuite
{
public:
    void TestSimpleMonodomainFastSlowPde() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        MixedTetrahedralMesh<2,2> mixed_mesh;
        mixed_mesh.ConstructRectangularMeshes(1.0,1.0,1,2);

        MyCellFactory cell_factory;
        cell_factory.SetMesh(mixed_mesh.GetFineMesh());

        MonodomainFastSlowPde<2> monodomain_fast_slow_pde( &cell_factory, mixed_mesh, 0, 1.0 );

        std::vector<FastSlowLuoRudyIModel1991*> cells(mixed_mesh.GetFineMesh()->GetNumNodes());
        for(unsigned i=0; i<mixed_mesh.GetFineMesh()->GetNumNodes(); i++)
        {
            cells[i] = static_cast<FastSlowLuoRudyIModel1991*>(monodomain_fast_slow_pde.mCellsDistributed[i]);
        }

        /////////////////////////////////////////////////////////////////
        // check the cells have been correctly assigned as fast or slow
        /////////////////////////////////////////////////////////////////
        TS_ASSERT_EQUALS( cells[0]->IsFast(), false );
        TS_ASSERT_EQUALS( cells[1]->IsFast(), true  );
        TS_ASSERT_EQUALS( cells[2]->IsFast(), false );
        TS_ASSERT_EQUALS( cells[3]->IsFast(), true  );
        TS_ASSERT_EQUALS( cells[4]->IsFast(), true  );
        TS_ASSERT_EQUALS( cells[5]->IsFast(), true  );
        TS_ASSERT_EQUALS( cells[6]->IsFast(), false );
        TS_ASSERT_EQUALS( cells[7]->IsFast(), true  );
        TS_ASSERT_EQUALS( cells[8]->IsFast(), false );

        std::vector<double> voltage_data(9, -84.5);
        Vec voltage = PetscTools::CreateVec(voltage_data);

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

        std::vector<double> slow_values(2);

        // test cells for which x=0: correspond to fine nodes 0 (coarse-slow cell), 3 (fine-fast cell), 6 (coarse-slow cell)
        TS_ASSERT_DELTA( cells[0]->GetVoltage(), cells[3]->GetVoltage(), 1.0); // check coarse and fine agree in voltage
        TS_ASSERT_DELTA( cells[0]->rGetStateVariables()[0], cells[3]->rGetStateVariables()[0], 0.01); // check coarse and fine agree in a gating var
        cells[0]->GetSlowValues(slow_values);
        TS_ASSERT_DELTA(slow_values[0], cells[3]->mSlowValues[0], 0.01); // check slow values match
        TS_ASSERT_DELTA(slow_values[1], cells[3]->mSlowValues[1], 0.01); // check slow values match

        // test cells for which x=0.5: correspond to fine nodes 1, 4, 7 (all fine-fast cells)
        TS_ASSERT_DELTA( cells[1]->GetVoltage(), cells[4]->GetVoltage(), 1.0);

        // test cells for which x=1: correspond to fine nodes 2 (coarse-slow), 5 (fine-fast), 8 (coarse-slow)
        TS_ASSERT_DELTA( cells[2]->GetVoltage(), cells[5]->GetVoltage(), 1.0); // check coarse and fine agree in voltage
        TS_ASSERT_DELTA( cells[2]->rGetStateVariables()[0], cells[5]->rGetStateVariables()[0], 0.01); // check coarse and fine agree in a gating var
        cells[2]->GetSlowValues(slow_values);
        TS_ASSERT_LESS_THAN(0, slow_values[0]);
        TS_ASSERT_DELTA(slow_values[0], cells[5]->mSlowValues[0], 0.01); // check slow values match
        TS_ASSERT_DELTA(slow_values[1], cells[5]->mSlowValues[1], 0.01); // check slow values match

        VecDestroy(voltage);
    }
};

#endif /*TESTMONODOMAINFASTSLOWPDE_HPP_*/


