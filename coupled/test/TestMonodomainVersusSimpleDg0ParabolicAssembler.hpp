#ifndef _TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_
#define _TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include <vector>

#include "SimpleLinearSolver.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "VectorPortion.hpp"


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
    
    virtual unsigned GetNumberOfCells()
    {
        return 1;
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
    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
        return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double u)
    {
        return u*(1-u);
    }
    
    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& , double u)
    {
        return u*(1-u);
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> )
    {
        return identity_matrix<double>(SPACE_DIM);
    }
    
    double ComputeDuDtCoefficientFunction(Point<SPACE_DIM> )
    {
        return 1;
    }
    
    FischerPde(ZeroStimCellFactory<SPACE_DIM>* pCellFactory)
            : MonodomainPde<SPACE_DIM>(pCellFactory)
    {}
};




class TestMonodomainDg0Assembler : public CxxTest::TestSuite
{
private:

    /**
     * Refactor code to set up a PETSc vector holding the initial condition.
     */
    Vec CreateInitialConditionVec(int size)
    {
        Vec initial_condition;
        VecCreate(PETSC_COMM_WORLD, &initial_condition);
        VecSetSizes(initial_condition, PETSC_DECIDE, size);
        VecSetFromOptions(initial_condition);
        return initial_condition;
    }
    
public:
    void TestMonodomainDg0AssemblerWithFischer1DAgainstSimpleDg0Assembler()
    {
        double t_start = 0;
        double t_final = 1;
        double pde_timestep = 0.01;
        
        
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/heart_FHN_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        ZeroStimCellFactory<1> cell_factory;
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
        SimpleDg0ParabolicAssembler<1,1> simple_assembler(&mesh,&pde,&bcc);
        
        // initial condition;
        Vec initial_condition_1, initial_condition_2;
        initial_condition_1 = CreateInitialConditionVec(mesh.GetNumNodes());
        VecDuplicate(initial_condition_1, &initial_condition_2);

        VectorPortion portion(initial_condition_1);
        for (VectorPortion::Iterator index = portion.Begin();
             index != portion.End();
             ++index)
        {
            double x=mesh.GetNode(index.Global)->GetPoint()[0];
            *index = exp(-(x*x)/100);
        }
        portion.Restore();
        
        VecCopy(initial_condition_1, initial_condition_2); // Both assemblers use same initial cond'n
        
        monodomain_assembler.SetTimes(t_start, t_final, pde_timestep);
        simple_assembler.SetTimes(t_start, t_final, pde_timestep);
        
        monodomain_assembler.SetInitialCondition( initial_condition_1 );
        simple_assembler.SetInitialCondition( initial_condition_2 );
        
        Vec current_solution_1 = monodomain_assembler.Solve();
        Vec current_solution_2 = simple_assembler.Solve();
        
        
        // Compare the results
        VectorPortion solution_1_portion(current_solution_1);
        VectorPortion solution_2_portion(current_solution_2);
        VectorPortion::Iterator index_1 = solution_1_portion.Begin();        
        VectorPortion::Iterator index_2 = solution_2_portion.Begin();        
        while (index_1 != solution_1_portion.End()
               && index_2 != solution_1_portion.End())
        {
            TS_ASSERT_DELTA(*index_1, *index_2, 1e-3);
            switch (index_1.Global)
            {
                case 10:
                    TS_ASSERT_DELTA(*index_1, 5.8028e-07, 1e-9);
                    break;
                case 25:
                    TS_ASSERT_DELTA(*index_1, 0.00648079, 1e-5);
                    break;
            }
            ++index_1;
            ++index_2;
        }
        TS_ASSERT_EQUALS(index_1, solution_1_portion.End());
        TS_ASSERT_EQUALS(index_2, solution_2_portion.End());
        
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
        SimpleDg0ParabolicAssembler<2,2> simple_assembler(&mesh,&pde,&bcc);
        
        // initial condition;
        Vec initial_condition_1, initial_condition_2;
        initial_condition_1 = CreateInitialConditionVec(mesh.GetNumNodes());
        VecDuplicate(initial_condition_1, &initial_condition_2);


        VectorPortion portion(initial_condition_1);
        for (VectorPortion::Iterator index = portion.Begin();
             index != portion.End();
             ++index)
        {
            double x=mesh.GetNode(index.Global)->GetPoint()[0];
            *index = exp(-(x*x)/100);
        }
        portion.Restore();
        
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
        
        // Compare the results
        VectorPortion solution_1_portion(current_solution_1);
        VectorPortion solution_2_portion(current_solution_2);
        VectorPortion::Iterator index_1 = solution_1_portion.Begin();        
        VectorPortion::Iterator index_2 = solution_2_portion.Begin();        
        while (index_1 != solution_1_portion.End()
               && index_2 != solution_1_portion.End())
        {
            TS_ASSERT_DELTA(*index_1, *index_2, 1e-3);
            ++index_1;
            ++index_2;
        }
        TS_ASSERT_EQUALS(index_1, solution_1_portion.End());
        TS_ASSERT_EQUALS(index_2, solution_2_portion.End());

        VecDestroy(current_solution_1);
        VecDestroy(current_solution_2);
    }
};

#endif //_TESTMONODOMAINVERSUSSIMPLEDG0PARABOLICASSEMBLER_HPP_
