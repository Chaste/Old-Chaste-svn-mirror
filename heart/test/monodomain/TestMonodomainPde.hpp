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


#ifndef _TESTMONODOMAINPDE_HPP_
#define _TESTMONODOMAINPDE_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <vector>

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
        else if (node==1)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }
        else
        {
            NEVER_REACHED;
            return NULL;
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
    
    void TestSaveAndLoadCardiacPDE()
    {
        //Archive                           
        OutputFileHandler handler("Archive", false);
        std::string archive_filename;
        handler.SetArchiveDirectory();
        archive_filename = handler.GetOutputDirectoryFullPath() + "bidomain_pde.arch";       

        {
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructLinearMesh(1);
    
            MyCardiacCellFactory cell_factory;
            cell_factory.SetMesh(&mesh);
    
            MonodomainPde<1> monodomain_pde( &cell_factory );            
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            MonodomainPde<1>* const p_archive_monodomain_pde = &monodomain_pde;
            output_arch << p_archive_monodomain_pde;
            
            ofs.close();
        }  

        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs); 
            
            MonodomainPde<1> *p_monodomain_pde;
            input_arch >> p_monodomain_pde;
                
            // Test rGetIntracellularConductivityTensor
            /// \todo: #98 Can't test this yet. Intracellular tensors have not been generated. Check second constructor of AbstractCardiacPde.
    //        const c_matrix<double, 1, 1>& tensor_before_archiving = monodomain_pde.rGetIntracellularConductivityTensor(1);
    //        const c_matrix<double, 1, 1>& tensor_after_archiving = p_monodomain_pde->rGetIntracellularConductivityTensor(1);
    //        
    //        TS_ASSERT(tensor_before_archiving(0,0) == tensor_after_archiving(0,0));
            
            // Test GetCardiacCell
            // Test rGetIionicCacheReplicated
            // Test rGetIntracellularStimulusCacheReplicated
    
            delete p_monodomain_pde;
            
            ifs.close();
        }        
    }
    
};



#endif //_TESTMONODOMAINPDE_HPP_
