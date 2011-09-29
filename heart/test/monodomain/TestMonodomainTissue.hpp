/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef _TESTMONODOMAINTISSUE_HPP_
#define _TESTMONODOMAINTISSUE_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <vector>

#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudy1991.hpp"
#include "MonodomainTissue.hpp"
#include "OdeSolution.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ArchiveOpener.hpp"
#include "DiFrancescoNoble1985.hpp"


class MyCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:

    MyCardiacCellFactory()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-80.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (node==0)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    boost::shared_ptr<SimpleStimulus> GetStimulus()
    {
        return mpStimulus;
    }
};

/*
 *  With this cardiac cell factory we cover the case where one tells Mono/BidomainProblem
 *  to use halo exchange but the cardiac cell factory is not aware that it may be asked to
 *  create cells at halo nodes (i.e. using GetNode() instead of GetNodeOrHaloNode())
 */
class MyWrongCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:

    MyWrongCardiacCellFactory()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-80.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (GetMesh()->GetNode(node)->rGetLocation()[0] < 0.5)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    boost::shared_ptr<SimpleStimulus> GetStimulus()
    {
        return mpStimulus;
    }
};


class PurkinjeCellFactory : public AbstractPurkinjeCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PurkinjeCellFactory()
        : AbstractPurkinjeCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-6000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        ChastePoint<2> location = GetMesh()->GetNode(node)->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    AbstractCardiacCell* CreatePurkinjeCellForTissueNode(unsigned node)
    {
        return new CellDiFrancescoNoble1985FromCellML(mpSolver, mpZeroStimulus);
    }
};




class SimpleConductivityModifier : public AbstractConductivityModifier<1,1>
{
private:
    c_matrix<double,1,1> mTensor;

public:
	SimpleConductivityModifier()
	    : AbstractConductivityModifier<1,1>()
	{
	}

	c_matrix<double,1,1>& rGetModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,1,1>& rOriginalConductivity)
	{
	    mTensor(0,0) = (elementIndex+2.0)*rOriginalConductivity(0,0); //so conductivity on element 0 gets scaled by 2, and by 3 on element 1
		return mTensor;
	}
};


class TestMonodomainTissue : public CxxTest::TestSuite
{
public:
    void TestMonodomainTissueBasic() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        unsigned num_nodes=2;
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, ie 2 node mesh
        assert(mesh.GetNumNodes()==num_nodes);

        double start_time = 0;
        double big_time_step = 0.5;

        boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        // Stimulus function to use at node 0. Node 1 is not stimulated.
        boost::shared_ptr<SimpleStimulus> p_stimulus = cell_factory.GetStimulus();
        boost::shared_ptr<ZeroStimulus> p_zero_stim(new ZeroStimulus);

        MonodomainTissue<1> monodomain_tissue( &cell_factory );

        // check the purkinje cells vector is empty
        TS_ASSERT(!monodomain_tissue.HasPurkinje());
		TS_ASSERT_THROWS_ANYTHING(monodomain_tissue.rGetPurkinjeCellsDistributed().size());
        TS_ASSERT_THROWS_ANYTHING(monodomain_tissue.rGetPurkinjeIionicCacheReplicated().GetSize());

        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;

        // initial condition;
        Vec voltage = PetscTools::CreateAndSetVec(num_nodes, initial_voltage);

        // Solve 1 (PDE) timestep using MonodomainTissue
        monodomain_tissue.SolveCellSystems(voltage, start_time, start_time+big_time_step);

        // Check results by solving ODE systems directly
        // Check node 0
        double value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[0];
        CellLuoRudy1991FromCellML ode_system_stimulated(p_solver, p_stimulus);
        ode_system_stimulated.ComputeExceptVoltage(start_time, start_time + big_time_step);
        double value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_tissue, value_ode, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[0];
        TS_ASSERT_DELTA(value_tissue, value_ode, 0.000001);

        // Check node 1
        CellLuoRudy1991FromCellML ode_system_not_stim(p_solver, p_zero_stim);
        value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[1];
        ode_system_not_stim.ComputeExceptVoltage(start_time, start_time + big_time_step);
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_tissue, value_ode, 0.000001);

        // Reset the voltage vector from ODE systems
        DistributedVector dist_voltage = mesh.GetDistributedVectorFactory()->CreateDistributedVector(voltage);
        for (DistributedVector::Iterator index = dist_voltage.Begin();
             index != dist_voltage.End();
             ++index)
        {
            if (index.Global==0)
            {
                dist_voltage[index] = ode_system_stimulated.rGetStateVariables()[0];
            }
            if (index.Global==1)
            {
                dist_voltage[index] = ode_system_not_stim.rGetStateVariables()[0];
            }
        }
        dist_voltage.Restore();

        // Use MonodomainTissue to solve a second (PDE) time step
        monodomain_tissue.SolveCellSystems(voltage, start_time, start_time+big_time_step);
        value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[0];

        // Check node 0 by solving ODE system directly
        ode_system_stimulated.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_tissue, value_ode, 1e-10);

        // Check node 1 by solving ODE system directly
        ode_system_not_stim.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[1];
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_tissue, value_ode, 1e-10);

        TS_ASSERT_THROWS_THIS(monodomain_tissue.rGetExtracellularConductivityTensor(0),
                              "Monodomain tissues do not have extracellular conductivity tensors.");

        VecDestroy(voltage);
    }

    void TestMonodomainTissueGetCardiacCell() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, ie 2 node mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory );

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0))
        {
            AbstractCardiacCell* cell = monodomain_tissue.GetCardiacCell(0);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),-80,1e-10);
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1))
        {
            AbstractCardiacCell* cell = monodomain_tissue.GetCardiacCell(1);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),0,1e-10);
        }
    }

    void TestSolveCellSystemsInclUpdateVoltage() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, ie 2 node mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory );

        Vec voltage = PetscTools::CreateAndSetVec(2, -81.4354); // something that isn't resting potential
        monodomain_tissue.SolveCellSystems(voltage, 0, 1, false); // solve for 1ms without updating the voltage

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0))
        {
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -81.4354, 1e-3);
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1))
        {
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -81.4354, 1e-3);
        }

        Vec voltage2 = PetscTools::CreateAndSetVec(2, -75);
        monodomain_tissue.SolveCellSystems(voltage2, 1, 2, true); // solve another ms, using this new voltage, but now updating the voltage too
        ReplicatableVector voltage2_repl(voltage2); // should have changed following solve

        // check the new voltage in the cell is NEAR -75 (otherwise the passed in voltage wasn't used, but
        // NOT EXACTLY -75, ie that the voltage was solved for.
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0))
        {
            // check has been updated
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -75);
            // check near -75
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -75, 2.0); // within 2mV
            // check the passed in voltage was updated
            TS_ASSERT_DELTA(voltage2_repl[0], monodomain_tissue.GetCardiacCell(0)->GetVoltage(), 1e-10);
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1))
        {
        	TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -75);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -75, 2.0); // within 2mV
            TS_ASSERT_DELTA(voltage2_repl[1], monodomain_tissue.GetCardiacCell(1)->GetVoltage(), 1e-10);
        }

        VecDestroy(voltage);
        VecDestroy(voltage2);
    }

    void TestGetConductivityAndConductivityModifier() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, ie 2 node mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory );

        double orig_conductivity_0 = monodomain_tissue.rGetIntracellularConductivityTensor(0)(0,0);
        double orig_conductivity_1 = monodomain_tissue.rGetIntracellularConductivityTensor(1)(0,0);
        TS_ASSERT_DELTA(orig_conductivity_0, 1.75, 1e-9); // hard-coded using default
        TS_ASSERT_DELTA(orig_conductivity_1, 1.75, 1e-9); // hard-coded using default

        SimpleConductivityModifier modifier;
        monodomain_tissue.SetConductivityModifier(&modifier);

        TS_ASSERT_DELTA(monodomain_tissue.rGetIntracellularConductivityTensor(0)(0,0), 2*orig_conductivity_0, 1e-9);
        TS_ASSERT_DELTA(monodomain_tissue.rGetIntracellularConductivityTensor(1)(0,0), 3*orig_conductivity_1, 1e-9);
    }

    void TestNodeExchange() throw(Exception)
    {
        HeartConfig::Instance()->Reset();

        HeartConfig::Instance()->Reset();
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.1, 1.0); // [0,1] with h=0.1, ie 11 node mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory, true );


        if ( PetscTools::GetNumProcs() == 1 )
        {
            TS_ASSERT_EQUALS( mesh.GetNumHaloNodes(), 0u );
        }
        else
        {
            if ( PetscTools::AmMaster() || PetscTools::AmTopMost() )
            {
                TS_ASSERT_EQUALS( mesh.GetNumHaloNodes(), 1u );
            }
            else
            {
                TS_ASSERT_EQUALS( mesh.GetNumHaloNodes(), 2u );
            }
        }

        for (DistributedTetrahedralMesh<1,1>::HaloNodeIterator it=mesh.GetHaloNodeIteratorBegin();
                it != mesh.GetHaloNodeIteratorEnd();
                ++it)
        {
            AbstractCardiacCell* cell = monodomain_tissue.GetCardiacCellOrHaloCell( (*it)->GetIndex() );
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),0,1e-10);
        }

        if ( PetscTools::AmMaster() )
        {
            // Master owns node 0
            AbstractCardiacCell* cell = monodomain_tissue.GetCardiacCellOrHaloCell(0);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001), -80.0, 1e-10);
        }
        else
        {
            // Zero is not halo owned by any process (unless we have a lot of them).
            TS_ASSERT_THROWS_CONTAINS(monodomain_tissue.GetCardiacCellOrHaloCell(0),
                                      "Requested node/halo 0 does not belong to processor ");
        }
    }

    void TestNodeExchangeExceptions() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=0.1, ie 11 node mesh

        MyWrongCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        if ( ! PetscTools::IsSequential() )
        {
            TS_ASSERT_THROWS_CONTAINS( MonodomainTissue<1> monodomain_tissue( &cell_factory, true ),
                                       "does not belong to processor" );
        }
    }

    void TestSolveCellSystemsInclUpdateVoltageWithNodeExchange() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=0.1, ie 11 node mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory, true );

        Vec voltage = PetscTools::CreateAndSetVec(2, -81.4354); // something that isn't resting potential
        monodomain_tissue.SolveCellSystems(voltage, 0, 1, false); // solve for 1ms without updating the voltage

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0)) // Is node 0 locally owned?
        {
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -81.4354, 1e-3);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCellOrHaloCell(1)->GetVoltage(), -81.4354, 1e-3);
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1)) // Is node 1 locally owned?
        {
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCellOrHaloCell(0)->GetVoltage(), -81.4354, 1e-3);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -81.4354, 1e-3);
        }

        Vec voltage2 = PetscTools::CreateAndSetVec(2, -75);
        monodomain_tissue.SolveCellSystems(voltage2, 1, 2, true); // solve another ms, using this new voltage, but now updating the voltage too
        ReplicatableVector voltage2_repl(voltage2); // should have changed following solve

        // check the new voltage in the cell is NEAR -75 (otherwise the passed in voltage wasn't used, but
        // NOT EXACTLY -75, ie that the voltage was solved for.
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0)) // Is node 0 locally owned?
        {
            // check has been updated
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -75);
            // check near -75
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -75, 2.0); // within 2mV
            // check the passed in voltage was updated
            TS_ASSERT_DELTA(voltage2_repl[0], monodomain_tissue.GetCardiacCell(0)->GetVoltage(), 1e-10);
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1)) // Is node 1 locally owned?
        {
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -75);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -75, 2.0); // within 2mV
            TS_ASSERT_DELTA(voltage2_repl[1], monodomain_tissue.GetCardiacCell(1)->GetVoltage(), 1e-10);
        }

        // now check the new voltages have been communicated
        // check the new voltage in the cell is NEAR -75 (otherwise the passed in voltage wasn't used, but
        // NOT EXACTLY -75, ie that the voltage was solved for.
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0)) // Is node 0 locally owned?
        {
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCellOrHaloCell(1)->GetVoltage(), -75);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCellOrHaloCell(1)->GetVoltage(), -75, 2.0); // within 2mV
            TS_ASSERT_DELTA(voltage2_repl[1], monodomain_tissue.GetCardiacCellOrHaloCell(1)->GetVoltage(), 1e-10);
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1)) // Is node 1 locally owned?
        {
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCellOrHaloCell(0)->GetVoltage(), -75);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCellOrHaloCell(0)->GetVoltage(), -75, 2.0); // within 2mV
            TS_ASSERT_DELTA(voltage2_repl[0], monodomain_tissue.GetCardiacCellOrHaloCell(0)->GetVoltage(), 1e-10);
        }

        VecDestroy(voltage);
        VecDestroy(voltage2);
    }

    void TestSaveAndLoadCardiacTissue() throw (Exception)
    {
        HeartConfig::Instance()->Reset();
        // Archive settings
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "monodomain_tissue.arch";

        bool cache_replication_saved = false;
        double saved_printing_timestep = 2.0;
        double default_printing_timestep = HeartConfig::Instance()->GetPrintingTimeStep();

        // Info about the first cell on this process (if any)
        bool has_cell = false;
        unsigned cell_v_index = (unsigned)(-1);
        double cell_v = DBL_MAX;

        c_matrix<double, 1, 1> tensor_before_archiving;
        {
            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            MyCardiacCellFactory cell_factory;
            cell_factory.SetMesh(&mesh);

            MonodomainTissue<1> monodomain_tissue( &cell_factory );
            monodomain_tissue.SetCacheReplication(cache_replication_saved); // Not the default to check it is archived...

            tensor_before_archiving = monodomain_tissue.rGetIntracellularConductivityTensor(1);

            // Get some info about the first cell on this process (if any)
            const std::vector<AbstractCardiacCell*>& r_cells = monodomain_tissue.rGetCellsDistributed();
            has_cell = !r_cells.empty();
            if (has_cell)
            {
                cell_v_index = r_cells[0]->GetVoltageIndex();
                cell_v = r_cells[0]->GetVoltage();
            }

            // Some checks to make sure HeartConfig is being saved and loaded by this too.
            HeartConfig::Instance()->SetPrintingTimeStep(saved_printing_timestep);
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);

            // Save
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacTissue<1>* const p_archive_monodomain_tissue = &monodomain_tissue;
            (*p_arch) << p_archive_monodomain_tissue;

            HeartConfig::Reset();
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), default_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep);
        }

        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacTissue<1>* p_monodomain_tissue;
            (*p_arch) >> p_monodomain_tissue;

            // Test rGetIntracellularConductivityTensor
            const c_matrix<double, 1, 1>& tensor_after_archiving = p_monodomain_tissue->rGetIntracellularConductivityTensor(1);
            TS_ASSERT_DELTA(tensor_before_archiving(0,0), tensor_after_archiving(0,0), 1e-9);

            TS_ASSERT_EQUALS(cache_replication_saved, p_monodomain_tissue->GetDoCacheReplication());
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep); // Test we are testing something in case default changes

            // Test cardiac cells have also been archived
            const std::vector<AbstractCardiacCell*>& r_cells = p_monodomain_tissue->rGetCellsDistributed();
            TS_ASSERT_EQUALS(has_cell, !r_cells.empty());
            if (has_cell)
            {
                TS_ASSERT_EQUALS(cell_v_index, r_cells[0]->GetVoltageIndex());
                TS_ASSERT_EQUALS(cell_v, r_cells[0]->GetVoltage());
            }

            delete p_monodomain_tissue;
        }
    }

    void TestMonodomainTissueUsingPurkinjeCellFactory() throw(Exception)
	{
		HeartConfig::Instance()->Reset();

		TrianglesMeshReader<2,2> reader("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
		MixedDimensionMesh<2,2> mixed_mesh;
		mixed_mesh.ConstructFromMeshReader(reader);

		PurkinjeCellFactory cell_factory;
		cell_factory.SetMesh(&mixed_mesh);

		MonodomainTissue<2> tissue( &cell_factory );

		TS_ASSERT(tissue.HasPurkinje());
		TS_ASSERT_EQUALS(tissue.rGetPurkinjeCellsDistributed().size(), tissue.rGetCellsDistributed().size());
        TS_ASSERT_EQUALS(tissue.rGetPurkinjeIionicCacheReplicated().GetSize(),
                         tissue.rGetIionicCacheReplicated().GetSize());

        for (AbstractTetrahedralMesh<2,2>::NodeIterator current_node = mixed_mesh.GetNodeIteratorBegin();
             current_node != mixed_mesh.GetNodeIteratorEnd();
             ++current_node)
        {
            unsigned global_index = current_node->GetIndex();
            AbstractCardiacCell* p_purkinje_cell = tissue.GetPurkinjeCell(global_index);
            double y = current_node->rGetLocation()[1];

            // cable nodes are on y=0.05 (we don't test by index because indices may be permuted in parallel).
            if( fabs(y-0.05) < 1e-8 )
            {
                TS_ASSERT(dynamic_cast<CellDiFrancescoNoble1985FromCellML*>(p_purkinje_cell) != NULL);
            }
            else
            {
                TS_ASSERT(dynamic_cast<FakeBathCell*>(p_purkinje_cell) != NULL);
            }

            TS_ASSERT_EQUALS(tissue.rGetPurkinjeCellsDistributed()[global_index-mixed_mesh.GetDistributedVectorFactory()->GetLow()],
                             p_purkinje_cell);
        }

        // Test archiving too
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "monodomain_tissue_purkinje.arch";

        { // Save to archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Make sure at least one Purkinje cell has a non-initial-condition state variable to compare
            if (mixed_mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0u))
            {
                tissue.GetPurkinjeCell(0u)->SetVoltage(1234.5);
            }

            AbstractCardiacTissue<2>* const p_archive_tissue = &tissue;
            (*p_arch) << p_archive_tissue;
        }

        { // Load from archive and compare
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacTissue<2>* p_tissue;
            (*p_arch) >> p_tissue;

            TS_ASSERT(p_tissue->HasPurkinje());
            TS_ASSERT_EQUALS(p_tissue->rGetPurkinjeCellsDistributed().size(), tissue.rGetPurkinjeCellsDistributed().size());
            TS_ASSERT_EQUALS(p_tissue->rGetPurkinjeIionicCacheReplicated().GetSize(), tissue.rGetPurkinjeIionicCacheReplicated().GetSize());

            for (AbstractTetrahedralMesh<2,2>::NodeIterator current_node = p_tissue->mpMesh->GetNodeIteratorBegin();
                 current_node != p_tissue->mpMesh->GetNodeIteratorEnd();
                 ++current_node)
            {
                unsigned global_index = current_node->GetIndex();
                AbstractCardiacCell* p_purkinje_cell = p_tissue->GetPurkinjeCell(global_index);
                double y = current_node->rGetLocation()[1];

                // cable nodes are on y=0.05 (we don't test by index because indices may be permuted in parallel).
                if( fabs(y-0.05) < 1e-8 )
                {
                    TS_ASSERT(dynamic_cast<CellDiFrancescoNoble1985FromCellML*>(p_purkinje_cell) != NULL);
                }
                else
                {
                    TS_ASSERT(dynamic_cast<FakeBathCell*>(p_purkinje_cell) != NULL);
                }
                TS_ASSERT_EQUALS(p_purkinje_cell->GetVoltage(),
                                 tissue.GetPurkinjeCell(global_index)->GetVoltage());

                TS_ASSERT_EQUALS(p_tissue->rGetPurkinjeCellsDistributed()[global_index-p_tissue->mpMesh->GetDistributedVectorFactory()->GetLow()],
                                 p_purkinje_cell);
            }

            delete p_tissue;
        }
	}
};



#endif //_TESTMONODOMAINTISSUE_HPP_
