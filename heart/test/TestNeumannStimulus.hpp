#ifndef _TESTNEUMANNSTIMULUS_HPP_
#define _TESTNEUMANNSTIMULUS_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ReplicatableVector.hpp"
#include "InitialStimulus.hpp"
#include "StimulusBoundaryCondition.hpp"
#include "ZeroStimulusCellFactory.hpp"


class TestNeumannStimulus : public CxxTest::TestSuite
{
public:

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainConstantStimulus() throw(Exception)
    {
        ZeroStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        monodomain_problem.SetEndTime(2);   // ms
        monodomain_problem.SetOutputDirectory("MonoNeumannConst");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");

        // this parameters are a bit arbitrary, and chosen to get a good spread of voltages
        monodomain_problem.SetIntracellularConductivities(Create_c_vector(1.75));
        monodomain_problem.Initialise();
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1*1.75/0.0005);
        
        // create boundary conditions container
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1> *p_bc_stim = new ConstBoundaryCondition<1>(2*1.75/0.0005);
        ConstBoundaryCondition<1> *p_bc_no_flux = new ConstBoundaryCondition<1>(0);
                
        // get mesh
        ConformingTetrahedralMesh<1,1> &mesh = monodomain_problem.rGetMesh();
        // loop over boundary elements
        ConformingTetrahedralMesh<1, 1>::BoundaryElementIterator iter;
        iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            if (((*iter)->GetNodeLocation(0))[0]==0.0)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_stim);
            }
            else
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_no_flux);
            }
            
            //AddDirichletBoundaryCondition(*iter, p_boundary_condition, indexOfUnknown);
            iter++;
        }
            // if the element is on the left of the mesh, add a stimulus to the bcc
        
        // pass the bcc to the monodomain problem
        monodomain_problem.SetBoundaryConditionsContainer(&bcc);

        monodomain_problem.Solve();
        
        // check some voltages    
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        double atol=5e-3;
        
        TS_ASSERT_DELTA(voltage_replicated[1], 94.6426, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 49.7867, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 30.5954, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 21.6782, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -33.9983, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -52.2396, atol);
              
    }
    
    void TestMonodomainSquareWaveStimulus() throw(Exception)
    {
        ZeroStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        monodomain_problem.SetEndTime(2);   // ms
        monodomain_problem.SetOutputDirectory("MonoNeumannSquare");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");        
        
        // this parameters are a bit arbitrary, and chosen to get a good spread of voltages
        monodomain_problem.SetIntracellularConductivities(Create_c_vector(1.75));
        monodomain_problem.Initialise();
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1*1.75/0.0005);
        
        // create boundary conditions container
        BoundaryConditionsContainer<1,1,1> bcc;
        InitialStimulus stim(4*1.75/0.0005, 0.5);
        StimulusBoundaryCondition<1> *p_bc_stim = new StimulusBoundaryCondition<1>(&stim);
        ConstBoundaryCondition<1> *p_bc_no_flux = new ConstBoundaryCondition<1>(0);
                
        // get mesh
        ConformingTetrahedralMesh<1,1> &mesh = monodomain_problem.rGetMesh();
        // loop over boundary elements
        ConformingTetrahedralMesh<1, 1>::BoundaryElementIterator iter;
        iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            if (((*iter)->GetNodeLocation(0))[0]==0.0)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_stim);
            }
            else
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_no_flux);
            }
            
            //AddDirichletBoundaryCondition(*iter, p_boundary_condition, indexOfUnknown);
            iter++;
        }
            // if the element is on the left of the mesh, add a stimulus to the bcc
        
        // pass the bcc to the monodomain problem
        monodomain_problem.SetBoundaryConditionsContainer(&bcc);
        
        monodomain_problem.Solve();
        
        // check some voltages    
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        double atol=5e-3;
        
        TS_ASSERT_DELTA(voltage_replicated[1], 22.4940, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 22.6008, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 23.3054, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.4932, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], 14.5184, atol);
        TS_ASSERT_DELTA(voltage_replicated[10],3.7081, atol);
    }


    void TestBidomain1d() throw(Exception)
    {
        ZeroStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        
        bidomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        bidomain_problem.SetEndTime(2);   // ms
        bidomain_problem.SetOutputDirectory("BiNeuman1d");
        bidomain_problem.SetOutputFilenamePrefix("results");        
        
        bidomain_problem.Initialise();
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1*1.75/0.0005);

        //Double check to confirm the default values for conductivity
        c_matrix<double, 1,1> intra_tensor=bidomain_problem.GetBidomainPde()->rGetIntracellularConductivityTensor(0);
        TS_ASSERT_DELTA(intra_tensor(0,0), 1.75, 1e-10);
        c_matrix<double, 1,1> extra_tensor=bidomain_problem.GetBidomainPde()->rGetExtracellularConductivityTensor(0);
        TS_ASSERT_DELTA(extra_tensor(0,0), 6.2, 1e-10);
        
        // create boundary conditions container
        BoundaryConditionsContainer<1,1,2> bcc;
        InitialStimulus stim(4*1.75/0.0005, 0.5);
        StimulusBoundaryCondition<1> *p_bc_stim = new StimulusBoundaryCondition<1>(&stim);
        ConstBoundaryCondition<1> *p_bc_no_flux = new ConstBoundaryCondition<1>(0);
                
        // get mesh
        ConformingTetrahedralMesh<1,1>& r_mesh = bidomain_problem.rGetMesh();
        // loop over boundary elements
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter;
        iter = r_mesh.GetBoundaryElementIteratorBegin();
        while (iter != r_mesh.GetBoundaryElementIteratorEnd())
        {
            if (((*iter)->CalculateCentroid()[0])==0.0)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_stim);
            }
            else
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_no_flux);
            }
            
            // add zero-flux boundary condition for the extracellular potential
            bcc.AddNeumannBoundaryCondition(*iter, p_bc_no_flux, 1);
            iter++;
        }
        
        // pass the bcc to the monodomain problem
        bidomain_problem.SetBoundaryConditionsContainer(&bcc);
        
        bidomain_problem.Solve();
        
        // check some voltages    
        ReplicatableVector voltage_replicated(bidomain_problem.GetVoltage());
        double atol=5e-3;
        
        TS_ASSERT_DELTA(voltage_replicated[2*1], 23.7349, atol);
        TS_ASSERT_DELTA(voltage_replicated[2*3], 23.4607, atol);
        TS_ASSERT_DELTA(voltage_replicated[2*5], 24.0685, atol);
        TS_ASSERT_DELTA(voltage_replicated[2*7], 19.4519, atol);
        TS_ASSERT_DELTA(voltage_replicated[2*9], -46.0072, atol);
        TS_ASSERT_DELTA(voltage_replicated[2*10], -64.1003, atol);
              
    }
    void TestBidomain2d() throw(Exception) 
    {
        ZeroStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory;
        BidomainProblem<2> bidomain_problem( &cell_factory );
        
        bidomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_200_elements");
        bidomain_problem.SetEndTime(2);   // ms
        bidomain_problem.SetOutputDirectory("BiNeuman2d");
        bidomain_problem.SetOutputFilenamePrefix("results");        
        bidomain_problem.SetIntracellularConductivities(Create_c_vector(1.75,1.75));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(6.2,6.2));
        
        bidomain_problem.Initialise();
        bidomain_problem.GetBidomainPde()->SetSurfaceAreaToVolumeRatio(1*1.75/0.0005);

        //Double check to confirm the default values for conductivity
        c_matrix<double, 2,2> intra_tensor=bidomain_problem.GetBidomainPde()->rGetIntracellularConductivityTensor(0);
        TS_ASSERT_DELTA(intra_tensor(0,0), 1.75, 1e-10);
        TS_ASSERT_DELTA(intra_tensor(1,1), 1.75, 1e-10);
        c_matrix<double, 2,2> extra_tensor=bidomain_problem.GetBidomainPde()->rGetExtracellularConductivityTensor(0);
        TS_ASSERT_DELTA(extra_tensor(0,0), 6.2, 1e-10);
        TS_ASSERT_DELTA(extra_tensor(1,1), 6.2, 1e-10);
        
        // create boundary conditions container
        BoundaryConditionsContainer<2, 2, 2> bcc;
        InitialStimulus stim(10000.0, 0.5);
        StimulusBoundaryCondition<2> *p_bc_stim = new StimulusBoundaryCondition<2>(&stim);
        ConstBoundaryCondition<2> *p_bc_no_flux = new ConstBoundaryCondition<2>(0);
                
        // get mesh
        ConformingTetrahedralMesh<2,2>& r_mesh = bidomain_problem.rGetMesh();
        // loop over boundary elements
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter;
        iter = r_mesh.GetBoundaryElementIteratorBegin();
        while (iter != r_mesh.GetBoundaryElementIteratorEnd())
        {
            if (((*iter)->CalculateCentroid()[0])==0.0)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_stim);
            }
            else
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_no_flux);
            }
            
            // add zero-flux boundary condition for the extracellular potential
            bcc.AddNeumannBoundaryCondition(*iter, p_bc_no_flux, 1);
            iter++;
        }
        
        // pass the bcc to the monodomain problem
        bidomain_problem.SetBoundaryConditionsContainer(&bcc);
        
        bidomain_problem.Solve();
        
        // check some voltages    
        ReplicatableVector voltage_replicated(bidomain_problem.GetVoltage());
        double atol=2.0;

        for (unsigned node_index = 0; node_index < r_mesh.GetNumNodes(); node_index++)
        {
            double x = r_mesh.GetNode(node_index)->rGetLocation()[0];
            
            if (fabs(x)<1e-10)
            {
                TS_ASSERT_DELTA(voltage_replicated[2*node_index], 21.0, atol)
            }
            if (fabs(x-0.05)<1e-10)
            {
                TS_ASSERT_DELTA(voltage_replicated[2*node_index], 23.5, atol)
            }
            if (fabs(x-0.1)<1e-10)
            {
                TS_ASSERT_DELTA(voltage_replicated[2*node_index], -68.5, 2*atol)
            }            
        } 
        
//        TS_ASSERT_DELTA(voltage_replicated[2*1], 23.6028, atol);
//        TS_ASSERT_DELTA(voltage_replicated[2*3], 23.3720, atol);
//        TS_ASSERT_DELTA(voltage_replicated[2*5], 23.9703, atol);
//        TS_ASSERT_DELTA(voltage_replicated[2*7], 20.7020, atol);
//        TS_ASSERT_DELTA(voltage_replicated[2*9], -41.3815, atol);
//        TS_ASSERT_DELTA(voltage_replicated[2*10], -60.6728, atol);
              
    }
};

#endif //_TESTNEUMANNSTIMULUS_HPP_
