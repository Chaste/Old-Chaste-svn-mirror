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


#ifndef TESTBIDOMAINPDE_HPP_
#define TESTBIDOMAINPDE_HPP_


#include <iostream>
#include <vector>

#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
#include "BidomainPde.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "TetrahedralMesh.hpp"
#include <petsc.h>
#include <cxxtest/TestSuite.h>


// cell factory for creating 2 cells with both intra and extracellular stimuli
class MyCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    AbstractStimulusFunction* mpStimulus;
    AbstractStimulusFunction* mpExtracellularStimulus1;
    AbstractStimulusFunction* mpExtracellularStimulus2;
public:

    MyCardiacCellFactory() : AbstractCardiacCellFactory<1>()
    {
        mpStimulus = new SimpleStimulus(-80.0, 0.5);
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (node==0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else if (node==1)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }
        else
        {
            NEVER_REACHED;
        }
    }

    ~MyCardiacCellFactory(void)
    {
        delete mpStimulus;
    }
};





class TestBidomainPde : public CxxTest::TestSuite
{
public:

    void TestBidomainPdeSolveCellSystems( void )
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(1);

        double big_time_step = 0.5;
        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainPde<1> monodomain_pde( &cell_factory );
        BidomainPde<1>     bidomain_pde( &cell_factory );

        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;

        // initial condition;
        Vec monodomain_vec = DistributedVector::CreateVec();
        DistributedVector monodomain_voltage(monodomain_vec);
        Vec bidomain_vec = DistributedVector::CreateVec(2);
        DistributedVector bidomain_ic(bidomain_vec);
        DistributedVector::Stripe bidomain_voltage(bidomain_ic,0);

        for (DistributedVector::Iterator index=DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            monodomain_voltage[index] = initial_voltage;
            bidomain_voltage[index] = initial_voltage;
        }

        monodomain_voltage.Restore();
        bidomain_ic.Restore();

        monodomain_pde.SolveCellSystems(monodomain_vec, 0, big_time_step);
        bidomain_pde.SolveCellSystems(bidomain_vec, 0, big_time_step);


        // Check that both the monodomain and bidomain PDE have the same ionic cache
        for (DistributedVector::Iterator index=DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            TS_ASSERT_EQUALS(monodomain_pde.rGetIionicCacheReplicated()[index.Global], bidomain_pde.rGetIionicCacheReplicated()[index.Global]);
        }

        // Check that the bidomain PDE has the right intracellular stimulus at node 0 and 1
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularStimulusCacheReplicated()[0], -80);
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularStimulusCacheReplicated()[1], 0);

        VecDestroy(monodomain_vec);
        VecDestroy(bidomain_vec);
    }
};

#endif /*TESTBIDOMAINPDE_HPP_*/
