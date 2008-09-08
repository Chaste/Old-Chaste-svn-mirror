/*

Copyright (C) University of Oxford, 2008

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


#ifndef TESTNONLINEARELASTICITYASSEMBLER_HPP_
#define TESTNONLINEARELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "NonlinearElasticityAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ExponentialMaterialLaw2.hpp"
#include "MooneyRivlinMaterialLaw2.hpp"

//#define TS_ASSERT_DELTA_DEBUG(a,b,c,d) TS_ASSERT_DELTA((a),(b),(c));if(fabs((a)-(b))>(c)){d;};

class TestNonlinearElasticityAssembler : public CxxTest::TestSuite
{
public:
    void TestAssembleSystem() throw (Exception)
    {
        QuadraticMesh<2> mesh("mesh/test/data/square_128_elements_quadratic");
        ExponentialMaterialLaw2<2> law(2,3);
        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        
        NonlinearElasticityAssembler<2> assembler(&mesh, 
                                                  &law, 
                                                  zero_vector<double>(2),
                                                  1.0,
                                                  fixed_nodes,
                                                  "");
        assembler.AssembleSystem(true, true);
        
        ///////////////////////////////////////////////////////////////////    
        // test whether residual vector is currently zero (as
        // current solution should have been initialised to u=0, p=p0
        ///////////////////////////////////////////////////////////////////    
        ReplicatableVector rhs_vec(assembler.mpLinearSystem->rGetRhsVector());
        TS_ASSERT_EQUALS( (int)rhs_vec.size(), 2*289+81 );
        for(unsigned i=0; i<rhs_vec.size(); i++)
        {
            TS_ASSERT_DELTA(rhs_vec[i], 0.0, 1e-12);
        }
        
        ///////////////////////////////////////////////////////////////////    
        // compute numerical jacobian and compare with analytic jacobian
        // (about u=0, p=p0)
        ///////////////////////////////////////////////////////////////////    
        unsigned num_dofs = rhs_vec.size();
        double h = 1e-6;
        
        int lo, hi;
        MatGetOwnershipRange(assembler.mpLinearSystem->rGetLhsMatrix(), &lo, &hi);
        
        for(unsigned j=0; j<num_dofs; j++)
        {
            assembler.mCurrentSolution.clear(); 
            assembler.FormInitialGuess();
            assembler.mCurrentSolution[j] += h;

            assembler.AssembleSystem(true, false);
            
            ReplicatableVector perturbed_rhs( assembler.mpLinearSystem->rGetRhsVector() );
            
            for(unsigned i=0; i<num_dofs; i++)
            {
                if((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = assembler.mpLinearSystem->GetMatrixElement(i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec[i])/h;
                    if((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error                     
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-2);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-2);
                    }
                }
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);

        //////////////////////////////////////////////////////////
        // compare numerical and analytic jacobians again, this
        // time using a non-zero displacement, u=lambda x, v = mu y 
        // (lambda not equal to 1/nu), p = p0.
        //////////////////////////////////////////////////////////
        double lambda = 1.2;
        double mu = 1.0/1.3;

        assembler.mCurrentSolution.clear(); 
        assembler.FormInitialGuess();
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            assembler.mCurrentSolution[2*i]   = (lambda-1)*mesh.GetNode(i)->rGetLocation()[0];
            assembler.mCurrentSolution[2*i+1] = (mu-1)*mesh.GetNode(i)->rGetLocation()[1];
        }
        
        assembler.AssembleSystem(true, true);
        ReplicatableVector rhs_vec2(assembler.mpLinearSystem->rGetRhsVector());

        h=1e-8; // needs to be smaller for this one
                
        for(unsigned j=0; j<num_dofs; j++)
        {
            assembler.mCurrentSolution[j] += h;
            assembler.AssembleSystem(true, false);
            
            ReplicatableVector perturbed_rhs( assembler.mpLinearSystem->rGetRhsVector() );
            
            for(unsigned i=0; i<num_dofs; i++)
            {
                if((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = assembler.mpLinearSystem->GetMatrixElement(i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec2[i])/h;
                    if((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error                     
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-2);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-2);
                    }
                }
            }
            assembler.mCurrentSolution[j] -= h;
        }
    }
    
    // A test where the solution should be zero displacement
    // It mainly tests that the initial guess was set up correctly to
    // the final correct solution, ie u=0, p=zero_strain_pressure (!=0)
    void TestWithZeroDisplacement() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools
        
        QuadraticMesh<2> mesh("mesh/test/data/square_128_elements_quadratic");

        double c1 = 3.0;
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(c1);

        std::vector<unsigned> fixed_nodes;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }

        NonlinearElasticityAssembler<2> assembler(&mesh,
                                                  &mooney_rivlin_law,
                                                  zero_vector<double>(2),
                                                  1.0,
                                                  fixed_nodes,
                                                  "");

        assembler.Solve();

        TS_ASSERT_EQUALS(assembler.GetNumNewtonIterations(), 0u);

        // get deformed position
        std::vector<c_vector<double,2> >& r_deformed_position
            = assembler.rGetDeformedPosition();

        for (unsigned i=0; i<r_deformed_position.size(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], r_deformed_position[i](0), 1e-8);
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[1], r_deformed_position[i](1), 1e-8);
        }

        // check the final pressure
        std::vector<double>& r_pressures = assembler.rGetPressures();
        TS_ASSERT_EQUALS(r_pressures.size(), mesh.GetNumVertices());
        for(unsigned i=0; i<r_pressures.size(); i++)
        {
            TS_ASSERT_DELTA(r_pressures[i], 2*c1, 1e-6);
        }
    }
    
    void TestSettingUpHeterogeneousProblem() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools
        
        // two element quad mesh on the square
        QuadraticMesh<2> mesh(1.0, 1.0, 1, 1);

        MooneyRivlinMaterialLaw2<2> law_1(1.0);
        MooneyRivlinMaterialLaw2<2> law_2(5.0);
        std::vector<AbstractIncompressibleMaterialLaw2<2>*> laws;
        laws.push_back(&law_1);
        laws.push_back(&law_2);
        
        std::vector<unsigned> fixed_nodes;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }

        NonlinearElasticityAssembler<2> assembler(&mesh, 
                                                  laws, 
                                                  zero_vector<double>(2),
                                                  1.0,
                                                  fixed_nodes,
                                                  "");
                                                  
        TS_ASSERT_EQUALS(assembler.mMaterialLaws.size(), 2u);
        TS_ASSERT_DELTA(assembler.mMaterialLaws[0]->GetZeroStrainPressure(), 2.0, 1e-6);
        TS_ASSERT_DELTA(assembler.mMaterialLaws[1]->GetZeroStrainPressure(), 10.0, 1e-6);
        
        unsigned num_nodes = 9; 
        // pressure for node 0 (in elem 0)
        TS_ASSERT_DELTA(assembler.mCurrentSolution[2*num_nodes + 0], 2.0, 1e-6); 
        // pressure for node 3 (in elem 1)
        TS_ASSERT_DELTA(assembler.mCurrentSolution[2*num_nodes + 3], 10.0, 1e-6); 
    }

    void TestSolve() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        QuadraticMesh<2> mesh("mesh/test/data/square_128_elements_quadratic");

        MooneyRivlinMaterialLaw2<2> law(0.02);
        c_vector<double,2> body_force;
        body_force(0) = 0.06;
        body_force(1) = 0.0;
        
        std::vector<unsigned> fixed_nodes;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }
        
        NonlinearElasticityAssembler<2> assembler(&mesh, 
                                                  &law, 
                                                  body_force,
                                                  1.0,
                                                  fixed_nodes,
                                                  "simple_nonlin_elas");
                                                  
        assembler.Solve();
        
        std::vector<c_vector<double,2> >& r_solution = assembler.rGetDeformedPosition();
        
        double xend = 1.17199;
        double yend = 0.01001;
        
        ///////////////////////////////////////////////////////////
        // compare the solution at the corners with the values 
        // obtained using the dealii finite elasticity assembler
        //        
        // Results have been visually checked to see they agree 
        // (they do, virtually or completely overlapping.
        ///////////////////////////////////////////////////////////

        // bottom lhs corner should still be at (0,0)
        assert( fabs(mesh.GetNode(0)->rGetLocation()[0] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(0)->rGetLocation()[1] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[0](0), 0.0, 1e-6 );
        TS_ASSERT_DELTA( r_solution[0](1), 0.0, 1e-6 );
        
        // top lhs corner should still be at (0,1)
        assert( fabs(mesh.GetNode(3)->rGetLocation()[0] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(3)->rGetLocation()[1] - 1) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[3](0), 0.0, 1e-6 );
        TS_ASSERT_DELTA( r_solution[3](1), 1.0, 1e-6 );

        // DEALII value for bottom rhs corner is (1.17199,0.01001)
        assert( fabs(mesh.GetNode(1)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(1)->rGetLocation()[1] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[1](0), xend, 1e-3 );
        TS_ASSERT_DELTA( r_solution[1](1), yend, 1e-3 );

        // DEALII value for top rhs corner is (1.17199,0.98999)
        assert( fabs(mesh.GetNode(2)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(2)->rGetLocation()[1] - 1) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[2](0), xend,   1e-3 );
        TS_ASSERT_DELTA( r_solution[2](1), 1-yend, 1e-3 );
    }
};

#endif /*TESTNONLINEARELASTICITYASSEMBLER_HPP_*/
