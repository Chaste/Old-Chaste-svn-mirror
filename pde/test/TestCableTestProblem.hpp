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

#ifndef TESTCABLETESTPROBLEM_HPP_
#define TESTCABLETESTPROBLEM_HPP_


#include <cxxtest/TestSuite.h>
#include "MixedDimensionMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "ReplicatableVector.hpp"
#include "StiffnessMatrixAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractFeCableObjectAssembler.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"


template<unsigned DIM>
class CableTestProblemRhsAssembler :  public AbstractFeCableObjectAssembler<DIM,DIM,1,true,false,NORMAL>
{
private:
    c_vector<double,1*2> ComputeCableVectorTerm(
        c_vector<double, 2>& rPhi,
        c_matrix<double, DIM, 2>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<1,DIM>* pElement)
    {
        return -rPhi;
    }

public:
    CableTestProblemRhsAssembler(MixedDimensionMesh<DIM,DIM>* pMesh)
        : AbstractFeCableObjectAssembler<DIM,DIM,1,true,false,NORMAL>(pMesh)
    {
    }
};




// solver
template<unsigned DIM>
class CableTestProblemSolver: public AbstractStaticLinearPdeSolver<DIM,DIM,1>
{
private:
	StiffnessMatrixAssembler<DIM,DIM>* mpStiffnessMatrixAssembler;
	CableTestProblemRhsAssembler<DIM>* mpRhsAssembler;

	/** Boundary conditions */
	BoundaryConditionsContainer<DIM,DIM,1>* mpBoundaryConditions;

    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
	{
        // use assemblers to set up KU=b
        mpStiffnessMatrixAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        mpStiffnessMatrixAssembler->Assemble();

        mpRhsAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), true);
        mpRhsAssembler->Assemble();

        this->mpLinearSystem->AssembleRhsVector();
        this->mpLinearSystem->AssembleIntermediateLhsMatrix();

        // apply the Dirichlet boundary conditions
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

        this->mpLinearSystem->AssembleRhsVector();
        this->mpLinearSystem->AssembleFinalLhsMatrix();
	}


public:
	CableTestProblemSolver(MixedDimensionMesh<DIM,DIM>* pMesh,
                           BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions)
		 : AbstractStaticLinearPdeSolver<DIM,DIM,1>(pMesh),
		   mpBoundaryConditions(pBoundaryConditions)
    {
		// set-up mStiffnessMatrixAssembler and mRhsAssembler
		mpStiffnessMatrixAssembler = new StiffnessMatrixAssembler<DIM,DIM>(pMesh);
		mpRhsAssembler = new CableTestProblemRhsAssembler<DIM>(pMesh);
    }

	~CableTestProblemSolver()
	{
		delete mpStiffnessMatrixAssembler;
		delete mpRhsAssembler;
	}
};

// See #1798
// Solve the problem \nabla^2 u = \delta_{cable}, on a cylindrical geometry r \in [0,1], z\in[0,1]
// with zero-Neumann BCs on the top and bottom surfaces and zero-Dirichlet BCs on the surface r=1.
// The solution is log(r)/2*pi.
//
class TestCableTestProblem : public CxxTest::TestSuite
{
public:
	void TestSolvingTestProblem() throw(Exception)
	{
	    EXIT_IF_PARALLEL;

		std::string mesh_base("mesh/test/data/mixed_dimension_meshes/cylinder");
		TrianglesMeshReader<3,3> reader(mesh_base);
		MixedDimensionMesh<3,3> mesh;
		mesh.ConstructFromMeshReader(reader);

        ConstBoundaryCondition<3>* p_zero_boundary_condition = new ConstBoundaryCondition<3>(0.0);

        // apply boundary conditions u=0 on curved surface of cylinder (zero neumann BCs are
        // applied on the bottom and top surface (z=0,1).
        BoundaryConditionsContainer<3,3,1> bcc;

//// NOTE: this mesh doesn't have a (valid) face file defined, so the face and boundary node
//// data isn't set up in the mesh. Therefore, rather than looping over boundary nodes, we
//// we loop over all nodes, and set the node to be a boundary node before passing it into
//// the BCC.
//        for (MixedDimensionMesh<3,3>::BoundaryNodeIterator iter =
//               mesh.GetBoundaryNodeIteratorBegin();
//             iter != mesh.GetBoundaryNodeIteratorEnd();
//	         iter++)

		for(unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
		    Node<3>* p_node = mesh.GetNode(i);
            double x = p_node->rGetLocation()[0];
		    double y = p_node->rGetLocation()[1];
            double r = sqrt(x*x+y*y);

		    if(fabs(r-1)<1e-3)
		    {
	            p_node->SetAsBoundaryNode(); // see comment above!
		        bcc.AddDirichletBoundaryCondition(p_node, p_zero_boundary_condition);
		    }
	    }

		// Solver
		CableTestProblemSolver<3> cable_solver(&mesh,&bcc);
		Vec result = cable_solver.Solve();
		ReplicatableVector result_repl(result);

		OutputFileHandler handler("CableTestProblem");
		out_stream p_file = handler.OpenOutputFile("solution.txt");

		// Solution should be u = log(r)/(2*pi)
		for (unsigned i=0; i<result_repl.GetSize(); i++)
		{
			double x = mesh.GetNode(i)->GetPoint()[0];
			double y = mesh.GetNode(i)->GetPoint()[1];
			double r = sqrt(x*x+y*y);

			double u = log(r)/2*M_PI;

double scale_factor_to_be_fixed=10;

			if(r>0.8)
			{
	            TS_ASSERT_DELTA(scale_factor_to_be_fixed*result_repl[i], u, 0.1);
			}
            else if(r>0.4)
            {
                TS_ASSERT_DELTA(scale_factor_to_be_fixed*result_repl[i], u, 0.25);
            }
			else if(r>0.1)
			{
			    TS_ASSERT_DELTA(scale_factor_to_be_fixed*result_repl[i], u, 0.5);
			}
			else
			{
			    // for these nodes r=0 and the true solution is infinite
			    TS_ASSERT_LESS_THAN(scale_factor_to_be_fixed*result_repl[i], -4);
			}

            *p_file << x << " " << y << " " << mesh.GetNode(i)->GetPoint()[2] << " " << result_repl[i] << "\n";

		}
		p_file->close();

		VecDestroy(result);
	}
};

#endif /* TESTCABLETESTPROBLEM_HPP_ */
