/*

Copyright (C) University of Oxford, 2005-2010

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


#ifndef _TESTMONODOMAINPDE_HPP_
#define _TESTMONODOMAINPDE_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <vector>

#include "MonodomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "ArchiveOpener.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"

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
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }
    }

    boost::shared_ptr<SimpleStimulus> GetStimulus()
    {
        return mpStimulus;
    }
};


class TestMonodomainPde : public CxxTest::TestSuite
{
public:
    void TestMonodomainPdeBasic( void )
    {
        unsigned num_nodes=2;
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(1);
        assert(mesh.GetNumNodes()==num_nodes);

        double start_time = 0;
        double big_time_step = 0.5;

        boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        // Stimulus function to use at node 0. Node 1 is not stimulated.
        boost::shared_ptr<SimpleStimulus> p_stimulus = cell_factory.GetStimulus();
        boost::shared_ptr<ZeroStimulus> p_zero_stim(new ZeroStimulus);

        MonodomainPde<1> monodomain_pde( &cell_factory );

        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;

        // initial condition;
        Vec voltage = PetscTools::CreateVec(num_nodes, initial_voltage);

        // Solve 1 (PDE) timestep using MonodomainPde
        monodomain_pde.SolveCellSystems(voltage, start_time, start_time+big_time_step);

        // Check results by solving ODE systems directly
        // Check node 0
        double value_pde = monodomain_pde.rGetIionicCacheReplicated()[0];
        LuoRudyIModel1991OdeSystem ode_system_stimulated(p_solver, p_stimulus);
        ode_system_stimulated.ComputeExceptVoltage(start_time, start_time + big_time_step);
        double value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value_pde = monodomain_pde.rGetIionicCacheReplicated()[0];
        TS_ASSERT_DELTA(value_pde, value_ode, 0.000001);

        // Check node 1
        LuoRudyIModel1991OdeSystem ode_system_not_stim(p_solver, p_zero_stim);
        value_pde = monodomain_pde.rGetIionicCacheReplicated()[1];
        ode_system_not_stim.ComputeExceptVoltage(start_time, start_time + big_time_step);
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 0.000001);

        // Reset the voltage vector from ODE systems
        DistributedVector dist_voltage = mesh.GetDistributedVectorFactory()->CreateDistributedVector(voltage);
        for (DistributedVector::Iterator index = dist_voltage.Begin();
             index != dist_voltage.End();
             ++index)
        {
            if (index.Global==0)
            {
                dist_voltage[index] = ode_system_stimulated.rGetStateVariables()[4];
            }
            if (index.Global==1)
            {
                dist_voltage[index] = ode_system_not_stim.rGetStateVariables()[4];
            }
        }
        dist_voltage.Restore();

        // Use MonodomainPde to solve a second (PDE) time step
        monodomain_pde.SolveCellSystems(voltage, start_time, start_time+big_time_step);
        value_pde = monodomain_pde.rGetIionicCacheReplicated()[0];

        // Check node 0 by solving ODE system directly
        ode_system_stimulated.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 1e-10);

        // Check node 1 by solving ODE system directly
        ode_system_not_stim.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_pde = monodomain_pde.rGetIionicCacheReplicated()[1];
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_pde, value_ode, 1e-10);

        VecDestroy(voltage);
    }

    void TestMonodomainPdeWithHeterogeneousConductivities() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        std::vector<ChasteCuboid<3> > heterogeneity_area;
        std::vector< c_vector<double,3> > intra_conductivities;
        std::vector< c_vector<double,3> > extra_conductivities;
        
        //first cuboid include element 0
        ChastePoint<3> cornerA(-1, -1, 0);
        ChastePoint<3> cornerB(0.1, 0.2, 0.2);
        ChasteCuboid<3> cuboid_1(cornerA, cornerB);
        heterogeneity_area.push_back(cuboid_1);
        
        //second cuboid include element 4
        ChastePoint<3> cornerC(0.11, 0.0, 0);
        ChastePoint<3> cornerD(0.2, 0.11, 0.2);
        ChasteCuboid<3> cuboid_2(cornerC, cornerD);
        
        heterogeneity_area.push_back(cuboid_2);
        
        //within the first area
        intra_conductivities.push_back( Create_c_vector(1.0, 2.0, 3.0) );   
        extra_conductivities.push_back( Create_c_vector(0.0, 0.0, 0.0) );

        //within the second area
        intra_conductivities.push_back( Create_c_vector(11.0, 22.0, 33.0) );   
        extra_conductivities.push_back( Create_c_vector(0.0, 0.0, 0.0) );
              
        HeartConfig::Instance()->SetConductivityHeterogeneities(heterogeneity_area, intra_conductivities, extra_conductivities); 
        
        //elsewhere
        double isotropic_conductivity=15.0;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(isotropic_conductivity, isotropic_conductivity, isotropic_conductivity));
        
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem,3> cell_factory;
        cell_factory.SetMesh(&mesh);
        
        //CreateIntracellularConductivityTensor called in the constructor
        MonodomainPde<3> monodomain_pde( &cell_factory );
        
        TS_ASSERT_EQUALS(monodomain_pde.rGetIntracellularConductivityTensor(0u)(0,0),1.0);//within first cuboid
        TS_ASSERT_EQUALS(monodomain_pde.rGetIntracellularConductivityTensor(4u)(0,0),11.0);//within second cuboid
        TS_ASSERT_EQUALS(monodomain_pde.rGetIntracellularConductivityTensor(4u)(1,1),22.0);//within second cuboid
        TS_ASSERT_EQUALS(monodomain_pde.rGetIntracellularConductivityTensor(8u)(0,0),15.0);//elsewhere, e.g. element 8
         
    }
    void TestMonodomainPdeGetCardiacCell( void )
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(1);

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainPde<1> monodomain_pde( &cell_factory );

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0))
        {
            AbstractCardiacCell* cell = monodomain_pde.GetCardiacCell(0);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),-80,1e-10);
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1))
        {
            AbstractCardiacCell* cell = monodomain_pde.GetCardiacCell(1);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),0,1e-10);
        }
    }

    void TestSaveAndLoadCardiacPDE() throw (Exception)
    {
        // Archive settings
        std::string archive_dir = "archive";
        std::string archive_file = "monodomain_pde.arch";

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

            MonodomainPde<1> monodomain_pde( &cell_factory );
            monodomain_pde.SetCacheReplication(cache_replication_saved); // Not the default to check it is archived...

            tensor_before_archiving = monodomain_pde.rGetIntracellularConductivityTensor(1);

            // Get some info about the first cell on this process (if any)
            const std::vector<AbstractCardiacCell*>& r_cells = monodomain_pde.rGetCellsDistributed();
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

            AbstractCardiacPde<1>* const p_archive_monodomain_pde = &monodomain_pde;
            (*p_arch) << p_archive_monodomain_pde;

            HeartConfig::Reset();
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), default_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep);
        }

        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacPde<1>* p_monodomain_pde;
            (*p_arch) >> p_monodomain_pde;

            // Test rGetIntracellularConductivityTensor
            const c_matrix<double, 1, 1>& tensor_after_archiving = p_monodomain_pde->rGetIntracellularConductivityTensor(1);
            TS_ASSERT_DELTA(tensor_before_archiving(0,0), tensor_after_archiving(0,0), 1e-9);

            TS_ASSERT_EQUALS(cache_replication_saved, p_monodomain_pde->GetDoCacheReplication());
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep); // Test we are testing something in case default changes

            // Test cardiac cells have also been archived
            const std::vector<AbstractCardiacCell*>& r_cells = p_monodomain_pde->rGetCellsDistributed();
            TS_ASSERT_EQUALS(has_cell, !r_cells.empty());
            if (has_cell)
            {
                TS_ASSERT_EQUALS(cell_v_index, r_cells[0]->GetVoltageIndex());
                TS_ASSERT_EQUALS(cell_v, r_cells[0]->GetVoltage());
            }

            delete p_monodomain_pde;
        }
    }
};



#endif //_TESTMONODOMAINPDE_HPP_
