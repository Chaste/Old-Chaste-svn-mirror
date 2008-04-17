#ifndef _TESTTISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_
#define _TESTTISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_

#include "UblasCustomFunctions.hpp"

#include <cxxtest/TestSuite.h>

#include <petsc.h>
#include <cmath>
#include <pde/test/pdes/SimplePoissonEquation.hpp>

#include "SimpleLinearEllipticAssembler.hpp"
#include "TissueSimulationWithNutrientsAssembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestTissueSimulationWithNutrientsAssembler : public CxxTest::TestSuite
{
public:
    
    void Test2dHeatEquationOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Instantiate PDE object
        SimplePoissonEquation<2> pde;
        
        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(2), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(3), p_boundary_condition);
        
        // Assembler
        SimpleLinearEllipticAssembler<2,2> simple_assembler(&mesh,&pde,&bcc);
        TissueSimulationWithNutrientsAssembler<2> nutrients_assembler(&mesh,&pde,&bcc);
        
        Vec simple_result = simple_assembler.Solve();
        Vec nutrients_result = nutrients_assembler.Solve();
        
        ReplicatableVector simple_result_repl(simple_result);
        ReplicatableVector nutrients_result_repl(nutrients_result);
        
        for (unsigned i=0; i<simple_result_repl.size(); i++)
        {
            TS_ASSERT_EQUALS(simple_result_repl[i], nutrients_result_repl[i]);
        }

        VecDestroy(simple_result);
        VecDestroy(nutrients_result);
    }
    
};

#endif //_TESTTISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_
