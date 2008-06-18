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

#ifndef TESTMONODOMAINFASTSLOWHEART_HPP_
#define TESTMONODOMAINFASTSLOWHEART_HPP_

#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleStimulus.hpp"
#include "FastSlowLuoRudyIModel1991.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainFastSlowProblem.hpp"

class StimulateApexCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    SimpleStimulus *mpStimulus;
    bool mFastSlow;
    
public:
    StimulateApexCellFactory(bool fastSlow) : AbstractCardiacCellFactory<3>(0.0025)
    {
        mpStimulus = new SimpleStimulus(-1000.0*500, 0.5);
        mFastSlow = fastSlow;
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        AbstractStimulusFunction* p_stimulus = mpZeroStimulus;
        
        // Stimulate the apex
        if (mpMesh->GetNode(node)->rGetLocation()[0] > 0.94)
        {
            p_stimulus = mpStimulus;
        }
        
        if(mFastSlow)
        {
            // fast-slow cells
            return new FastSlowLuoRudyIModel1991(mpSolver, mTimeStep, p_stimulus); // state unset at the moment
        }
        else
        {
            // normal cells
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, p_stimulus);
        }
    }

    ~StimulateApexCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestMonodomainFastSlowHeart : public CxxTest::TestSuite
{
public:

    void TestHeartMonodomainFastSlowProblemAgainstNormal() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        double simulation_time = 0.0001;

        /*
         *  Read meshes from disk
         */
        TrianglesMeshReader<3,3> fine_mesh_reader("heart/test/data/heart");
        ConformingTetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructFromMeshReader(fine_mesh_reader);

        TrianglesMeshReader<3,3> coarse_mesh_reader("heart/test/data/HeartDecimation_1545nodes");
        MixedTetrahedralMesh<3,3> mixed_mesh;
        mixed_mesh.ConstructFromMeshReader(coarse_mesh_reader);

        mixed_mesh.SetFineMesh(&fine_mesh);


        //////////////////////////////////////////////////
        // solve a mixed mesh, fast/slow problem
        //////////////////////////////////////////////////
        ReplicatableVector voltage_fast_slow;
        {

    
            StimulateApexCellFactory cell_factory(true);
    
            MonodomainFastSlowProblem<3> monodomain_fast_slow_prob( &cell_factory, mixed_mesh, 1.0 );
    
            monodomain_fast_slow_prob.SetEndTime(simulation_time);   // ms
            monodomain_fast_slow_prob.SetOutputDirectory("HeartMonodomainFastSlow");
            monodomain_fast_slow_prob.SetOutputFilenamePrefix("res");
    
            monodomain_fast_slow_prob.Initialise();
            monodomain_fast_slow_prob.Solve();
    
            EventHandler::Headings();
            EventHandler::Report();
    
            voltage_fast_slow.ReplicatePetscVector( monodomain_fast_slow_prob.GetVoltage() );
            TS_ASSERT_EQUALS(voltage_fast_slow.size(), mixed_mesh.GetFineMesh()->GetNumNodes() );

            EventHandler::Reset();
        }

        //////////////////////////////////////////////////
        // solve using normal monodomain problem
        //////////////////////////////////////////////////
        ReplicatableVector voltage_normal;
        {
            StimulateApexCellFactory cell_factory_normal(false);
    
            MonodomainProblem<3> monodomain_prob( &cell_factory_normal);
            monodomain_prob.SetMesh(mixed_mesh.GetFineMesh());
    
            monodomain_prob.SetEndTime(simulation_time);   // ms
            monodomain_prob.SetOutputDirectory("HeartMonodomainNormalToCompareWithFastSlow");
            monodomain_prob.SetOutputFilenamePrefix("res");
    
            monodomain_prob.Initialise();
            monodomain_prob.Solve();
            EventHandler::Headings();
            EventHandler::Report();
            EventHandler::Reset();
            voltage_normal.ReplicatePetscVector( monodomain_prob.GetVoltage() );
        }

        //////////////////////////////////////////////////
        // compare 
        //////////////////////////////////////////////////
        TS_ASSERT_EQUALS(voltage_fast_slow.size(), voltage_normal.size() );

        bool some_voltage_greater_than_zero = false;
        for(unsigned i=0; i<voltage_fast_slow.size(); i++)
        {
            TS_ASSERT_DELTA(voltage_fast_slow[i], voltage_normal[i], 1.0);
            if(voltage_fast_slow[i] > 0.0)
            {
                some_voltage_greater_than_zero = true;
            }
        }
        TS_ASSERT(some_voltage_greater_than_zero);

        EventHandler::Enable();
    }
    
};

#endif /*TESTMONODOMAINFASTSLOWHEART_HPP_*/
