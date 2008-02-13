#ifndef _TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include <vector>

#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"


template <int SPACE_DIM>
class ZeroStimCellFactory : public AbstractCardiacCellFactory<SPACE_DIM>
{
public:
    ZeroStimCellFactory() : AbstractCardiacCellFactory<SPACE_DIM>(0.01)
    {}
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(this->mpSolver, this->mTimeStep, this->mpZeroStimulus);
    }
};




/*
 * A simple parabolic PDE used in this test.
 *
 * NOTE: The fischer pde is a parabolic linear pde, however we want to get the
 * MonodomainDg0Assembler, which expects a MonodomainPde passed it, to solve it.
 * Therefore, we get this pde to inherit from MonodomainPde, and overload all the
 * main functions. For this reason it has to take in a cell class, although this is
 * not used.
 */
template <int SPACE_DIM>
class FischerPde : public MonodomainPde<SPACE_DIM>
{
public:
    double ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& )
    {
        return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>&, double u)
    {
        return u*(1-u);
    }
    
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& , double u)
    {
        return u*(1-u);
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM> & , Element<SPACE_DIM,SPACE_DIM>* pElement)
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& )
    {
        return 1;
    }
    
    FischerPde(ZeroStimCellFactory<SPACE_DIM>* pCellFactory)
            : MonodomainPde<SPACE_DIM>(pCellFactory)
    {}
};




class TestMonodomainVersusSimpleDg0ParabolicAssembler : public CxxTest::TestSuite
{
public:
    void TestMonodomainDg0AssemblerWithFischer1DAgainstSimpleDg0Assembler()
    {
        double t_start = 0;
        double t_final = 1;
        double pde_timestep = 0.01;

        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("heart/test/data/heart_FHN_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        ZeroStimCellFactory<1> cell_factory;
        cell_factory.SetMesh(&mesh);
        FischerPde<1> pde(&cell_factory);
        
        // Boundary conditions: zero neumann on entire boundary (2 elements)
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_zero_condition = new ConstBoundaryCondition<1>(0.0);
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_zero_condition);
            iter++;
        }
        
        // Assembler
        MonodomainDg0Assembler<1,1> monodomain_assembler(&mesh,&pde);
        monodomain_assembler.SetBoundaryConditionsContainer(&bcc);
        SimpleDg0ParabolicAssembler<1,1, true> simple_assembler(&mesh,&pde,&bcc);
        
        // initial condition;
        Vec initial_condition_1, initial_condition_2;
        DistributedVector::SetProblemSize(mesh.GetNumNodes());
        initial_condition_1 = DistributedVector::CreateVec();
        VecDuplicate(initial_condition_1, &initial_condition_2);
        
        DistributedVector dist_ic(initial_condition_1);
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            double x=mesh.GetNode(index.Global)->GetPoint()[0];
            dist_ic[index] = exp(-(x*x)/100);
        }
        dist_ic.Restore();
        VecCopy(initial_condition_1, initial_condition_2); // Both assemblers use same initial cond'n
        
        monodomain_assembler.SetTimes(t_start, t_final, pde_timestep);
        simple_assembler.SetTimes(t_start, t_final, pde_timestep);
        
        monodomain_assembler.SetInitialCondition( initial_condition_1 );
        simple_assembler.SetInitialCondition( initial_condition_2 );
        
        Vec current_solution_1 = monodomain_assembler.Solve();
        Vec current_solution_2 = simple_assembler.Solve();

        DistributedVector dist_sol_1(current_solution_1);
        DistributedVector dist_sol_2(current_solution_2);
        
        // Compare the results      
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            TS_ASSERT_DELTA(dist_sol_1[index], dist_sol_2[index], 1e-3);
            if (index.Global==10) TS_ASSERT_DELTA(dist_sol_1[index], 5.8028e-07, 5e-8);
            if (index.Global==25) TS_ASSERT_DELTA(dist_sol_1[index], 0.00648079, 3e-5);
            if (index.Global==50) TS_ASSERT_DELTA(dist_sol_1[index], 0.992718, 1e-5);
            if (index.Global==75) TS_ASSERT_DELTA(dist_sol_1[index], 0.00648079, 3e-5);
        }
        
        VecDestroy(initial_condition_1);
        VecDestroy(initial_condition_2);
        VecDestroy(current_solution_1);
        VecDestroy(current_solution_2);
    }
    
    void TestMonodomainDg0AssemblerWithFischer2DAgainstSimpleDg0Assembler()
    {
        double t_start = 0;
        double t_final = 1;
        double pde_timestep = 0.01;
        
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        ZeroStimCellFactory<2> cell_factory;
        cell_factory.SetMesh(&mesh);
        FischerPde<2> pde(&cell_factory);
        
        // Boundary conditions: zero neumann on entire boundary (2 elements)
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_neumann_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);
            iter++;
        }
        
        // Assembler
        MonodomainDg0Assembler<2,2> monodomain_assembler(&mesh,&pde);
        monodomain_assembler.SetBoundaryConditionsContainer(&bcc);
        SimpleDg0ParabolicAssembler<2,2, true> simple_assembler(&mesh,&pde,&bcc);
        
        // initial condition;
        Vec initial_condition_1, initial_condition_2;
        DistributedVector::SetProblemSize(mesh.GetNumNodes());
        initial_condition_1 = DistributedVector::CreateVec();
        VecDuplicate(initial_condition_1, &initial_condition_2);
        
        DistributedVector dist_ic(initial_condition_1);
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            double x=mesh.GetNode(index.Global)->GetPoint()[0];
            dist_ic[index] = exp(-(x*x)/100);
        }
        dist_ic.Restore();
        VecCopy(initial_condition_1, initial_condition_2); // Both assemblers use same initial cond'n
        
        // Vars to hold current solutions at each iteration
        Vec current_solution_1, current_solution_2;
        
        double tCurrent = t_start;
        while ( tCurrent < t_final )
        {
            monodomain_assembler.SetTimes(tCurrent, tCurrent+pde_timestep, pde_timestep);
            simple_assembler.SetTimes(tCurrent, tCurrent+pde_timestep, pde_timestep);
            
            monodomain_assembler.SetInitialCondition( initial_condition_1 );
            simple_assembler.SetInitialCondition( initial_condition_2 );
            
            current_solution_1 = monodomain_assembler.Solve();
            current_solution_2 = simple_assembler.Solve();
            
            // Next iteration uses current solution as initial condition
            VecDestroy(initial_condition_1); // Free old initial condition
            VecDestroy(initial_condition_2);
            initial_condition_1 = current_solution_1;
            initial_condition_2 = current_solution_2;
            
            tCurrent += pde_timestep;
        }

        DistributedVector dist_sol_1(current_solution_1);
        DistributedVector dist_sol_2(current_solution_2);
        
        // Compare the results
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            TS_ASSERT_DELTA(dist_sol_1[index], dist_sol_2[index], 1e-3);
        }

        VecDestroy(current_solution_1);
        VecDestroy(current_solution_2);

    }
};

#endif //_TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_
