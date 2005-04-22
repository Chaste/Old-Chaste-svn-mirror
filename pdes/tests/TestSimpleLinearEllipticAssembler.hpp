#ifndef _TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "LinearHeatEquationPde.hpp"
#include "SimpleLinearSolver.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include "Node.hpp"
#include "Element.hpp"

class TestSimpleLinearEllipticAssembler : public CxxTest::TestSuite
{
	public:
	
	void TestAssembleOnElement( void )
	{
		LinearHeatEquationPde<1> pde;
		std::vector<Node<1>*> nodes;
		nodes.push_back(new Node<1>(0, false, 1.0));
		nodes.push_back(new Node<1>(0, false, 3));
		Element<1,1> element(nodes);
		LinearBasisFunction<1> basis_function;
		MatrixDouble ael(2,2);
		VectorDouble bel(2);
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		assembler.AssembleOnElement(element, ael, bel, &pde, basis_function);
		
		TS_ASSERT_DELTA(ael(0,0),0.5, 1e-12);
		TS_ASSERT_DELTA(ael(0,1),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,1),0.5, 1e-12);
		
		TS_ASSERT_DELTA(bel(0),1, 1e-12);
		TS_ASSERT_DELTA(bel(1),1, 1e-12);
		
	}
	
	
	void TestWithHeatEquation()
	{
		PetscInitialize(0, NULL, 0, 0);
		
		// Create mesh (by hand!)
		const int num_elements = 10;
		ConformingTetrahedralMesh<1,1> mesh(num_elements);
		std::vector<Node<1>*> nodes;
		for (int i=0; i<num_elements+1; i++)
		{
			nodes.push_back(new Node<1>(i, false, 0.0 + 0.15*i));
			mesh.AddNode(*nodes[i]);
		}
		for (int i=0; i<num_elements; i++)
		{
			std::vector<Node<1>*> element_nodes;
			element_nodes.push_back(nodes[i]);
			element_nodes.push_back(nodes[i+1]);
			Element<1,1> element(element_nodes);
			mesh.AddElement(element);
		}
		
		// Instantiate PDE object
		LinearHeatEquationPde<1> pde;
		
		// Boundary conditions
		
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, /*bcs,*/ &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution should be u = 0.5*x*(3-x)
		for (int i=0; i < num_elements+1; i++)
		{
			double x = 0.0 + 0.15*i;
			double u = 0.5*x*(3-x);
			TS_ASSERT_DELTA(res[i], u, 0.001);
		}
		VecRestoreArray(result, &res);
//		PetscFinalize();
	}
};
                      
#endif //_TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_
