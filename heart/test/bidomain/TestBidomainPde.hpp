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

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <iostream>
#include <vector>

//#include "ArchiveLocationInfo.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
#include "BidomainPde.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
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
    
    void TestSaveAndLoadCardiacPDE()
    {
        //Archive                           
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("bidomain_pde.arch");

        {
            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
    
            MyCardiacCellFactory cell_factory;
            cell_factory.SetMesh(&mesh);
    
            BidomainPde<1> bidomain_pde( &cell_factory );            
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            BidomainPde<1>* const p_archive_bidomain_pde = &bidomain_pde;
            output_arch << p_archive_bidomain_pde;  
    
            ofs.close();
        }

        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs); 
            
            BidomainPde<1> *p_bidomain_pde;
            input_arch >> p_bidomain_pde; 
            
            delete p_bidomain_pde;
            
            ifs.close();
        }
    }
};

#endif /*TESTBIDOMAINPDE_HPP_*/
