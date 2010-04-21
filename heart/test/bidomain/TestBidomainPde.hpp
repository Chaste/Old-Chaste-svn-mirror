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


#ifndef TESTBIDOMAINPDE_HPP_
#define TESTBIDOMAINPDE_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <iostream>
#include <vector>

#include "ArchiveOpener.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
#include "BidomainPde.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "TetrahedralMesh.hpp"
#include <petsc.h>


// cell factory for creating 2 cells with both intra and extracellular stimuli
class MyCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<AbstractStimulusFunction> mpStimulus;
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

    ~MyCardiacCellFactory(void)
    {
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
        DistributedVectorFactory* p_factory = mesh.GetDistributedVectorFactory();
        Vec monodomain_vec = p_factory->CreateVec();
        DistributedVector monodomain_voltage = p_factory->CreateDistributedVector(monodomain_vec);
        Vec bidomain_vec = p_factory->CreateVec(2);
        DistributedVector bidomain_ic = p_factory->CreateDistributedVector(bidomain_vec);
        DistributedVector::Stripe bidomain_voltage(bidomain_ic,0);

        for (DistributedVector::Iterator index=monodomain_voltage.Begin();
             index != monodomain_voltage.End();
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
        for (unsigned node_index = mesh.GetDistributedVectorFactory()->GetLow();
             node_index < mesh.GetDistributedVectorFactory()->GetHigh();
             node_index++)
        {
            TS_ASSERT_EQUALS(monodomain_pde.rGetIionicCacheReplicated()[node_index], bidomain_pde.rGetIionicCacheReplicated()[node_index]);
        }

        // Check that the bidomain PDE has the right intracellular stimulus at node 0 and 1
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularStimulusCacheReplicated()[0], -80);
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularStimulusCacheReplicated()[1], 0);

        VecDestroy(monodomain_vec);
        VecDestroy(bidomain_vec);
    }
    
    void TestBidomainPdeWithHeterogeneousConductivities() throw (Exception)
    {
        HeartConfig::Instance()->Reset();
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
        extra_conductivities.push_back( Create_c_vector(51.0, 52.0, 53.0) );

        //within the second area
        intra_conductivities.push_back( Create_c_vector(11.0, 22.0, 33.0) );   
        extra_conductivities.push_back( Create_c_vector(151.0, 152.0, 153.0) );
              
        HeartConfig::Instance()->SetConductivityHeterogeneities(heterogeneity_area, intra_conductivities, extra_conductivities); 
        
        
        //elsewhere
        double isotropic_intra_conductivity=15.0;
        double isotropic_extra_conductivity=65.0;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(isotropic_intra_conductivity, isotropic_intra_conductivity, isotropic_intra_conductivity));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(isotropic_extra_conductivity, isotropic_extra_conductivity, isotropic_extra_conductivity));
        
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem,3> cell_factory_for_het;
        cell_factory_for_het.SetMesh(&mesh);
        
        //CreateIntracellularConductivityTensor called in the constructor
        BidomainPde<3> bidomain_pde( &cell_factory_for_het );
        
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularConductivityTensor(0u)(0,0),1.0);//within first cuboid
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularConductivityTensor(4u)(0,0),11.0);//within second cuboid
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularConductivityTensor(4u)(1,1),22.0);//within second cuboid
        TS_ASSERT_EQUALS(bidomain_pde.rGetIntracellularConductivityTensor(8u)(0,0),15.0);//elsewhere, e.g. element 8
        
        TS_ASSERT_EQUALS(bidomain_pde.rGetExtracellularConductivityTensor(0u)(0,0),51.0);//within first cuboid
        TS_ASSERT_EQUALS(bidomain_pde.rGetExtracellularConductivityTensor(4u)(0,0),151.0);//within second cuboid
        TS_ASSERT_EQUALS(bidomain_pde.rGetExtracellularConductivityTensor(4u)(1,1),152.0);//within second cuboid
        TS_ASSERT_EQUALS(bidomain_pde.rGetExtracellularConductivityTensor(8u)(0,0),65.0);//elsewhere, e.g. element 8
         
    }

    void TestSaveAndLoadCardiacPDE()
    {
        HeartConfig::Instance()->Reset();
        // Archive settings
        std::string archive_dir = "archive";
        std::string archive_file = "bidomain_pde.arch";

        bool cache_replication_saved = false;
        double saved_printing_timestep = 2.0;
        double default_printing_timestep = HeartConfig::Instance()->GetPrintingTimeStep();

        c_matrix<double, 1, 1> intra_tensor_before_archiving;
        c_matrix<double, 1, 1> extra_tensor_before_archiving;
        {
            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            MyCardiacCellFactory cell_factory;
            cell_factory.SetMesh(&mesh);

            BidomainPde<1> bidomain_pde( &cell_factory );
            bidomain_pde.SetCacheReplication(cache_replication_saved); // Not the default to check it is archived...

            // Some checks to make sure HeartConfig is being saved and loaded by this too.
            HeartConfig::Instance()->SetPrintingTimeStep(saved_printing_timestep);
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);

            intra_tensor_before_archiving = bidomain_pde.rGetIntracellularConductivityTensor(1);
            extra_tensor_before_archiving = bidomain_pde.rGetExtracellularConductivityTensor(1);

            // Save
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacPde<1>* const p_archive_bidomain_pde = &bidomain_pde;
            (*p_arch) << p_archive_bidomain_pde;

            HeartConfig::Reset();
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), default_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep);
        }

        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacPde<1>* p_bidomain_pde;
            (*p_arch) >> p_bidomain_pde;

            const c_matrix<double, 1, 1>& intra_tensor_after_archiving = p_bidomain_pde->rGetIntracellularConductivityTensor(1);
            TS_ASSERT_DELTA(intra_tensor_before_archiving(0,0), intra_tensor_after_archiving(0,0), 1e-9);
            const c_matrix<double, 1, 1>& extra_tensor_after_archiving = dynamic_cast<BidomainPde<1>*>(p_bidomain_pde)->rGetExtracellularConductivityTensor(1); //Naughty Gary using dynamic cast, but only for testing...
            TS_ASSERT_DELTA(extra_tensor_before_archiving(0,0), extra_tensor_after_archiving(0,0), 1e-9);

            TS_ASSERT_EQUALS(cache_replication_saved, p_bidomain_pde->GetDoCacheReplication());
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep); // Test we are testing something in case default changes

            delete p_bidomain_pde;
        }
    }
};

#endif /*TESTBIDOMAINPDE_HPP_*/
