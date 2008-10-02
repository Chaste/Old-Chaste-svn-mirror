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


#ifndef TESTIMPLICITCARDIACMECHANICSASSEMBLER2_HPP_
#define TESTIMPLICITCARDIACMECHANICSASSEMBLER2_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "ImplicitCardiacMechanicsAssembler2.hpp"
#include "MooneyRivlinMaterialLaw2.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestImplicitCardiacMechanicsAssembler2 : public CxxTest::TestSuite
{
public:
    void TestCompareJacobians() throw(Exception)
    {
        QuadraticMesh<2> mesh(1.0, 1.0, 1, 1);
        MooneyRivlinMaterialLaw2<2> law(0.02);
        
        std::vector<unsigned> fixed_nodes;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }

        ImplicitCardiacMechanicsAssembler2<2> assembler(&mesh,"",fixed_nodes,&law);

        std::vector<double> calcium_conc(assembler.GetTotalNumQuadPoints(), 0.0);

        for(unsigned i=0; i<calcium_conc.size(); i++)
        {
            calcium_conc[i] = 0.05;
        }

        assembler.SetIntracellularCalciumConcentrations(calcium_conc);

        // NOTE: calling CompareJacobians below bypasses calling Solve(t0,t1,dt), hence the
        // time info will not be set. We therefore must explicitly set them here.
        assembler.mCurrentTime = 0.0;
        assembler.mNextTime = 0.01;
        assembler.mOdeTimestep = 0.01;

   
        ///////////////////////////////////////////////////////////////////    
        // compute numerical jacobian and compare with analytic jacobian
        // (about u=0, p=p0)
        ///////////////////////////////////////////////////////////////////    
        assembler.AssembleSystem(true, true);
        ReplicatableVector rhs_vec(assembler.mpLinearSystem->rGetRhsVector());
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
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-4);
                    }
                }
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }


//
//    void deleteTestCompareWithExplicitSolver()
//    {
//        // solve an implicit deformation, with some [Ca] forcing term. Then
//        // get the active tension, pass that into an explicit solver, and
//        // check the result is the same.
//        Triangulation<2> mesh;
//        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
//        mesh.refine_global(3);
//
//        Point<2> zero;
//        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
//
//        // specify this material law so the test continues to pass when the default
//        // material law is changed.
//        MooneyRivlinMaterialLaw<2> material_law(0.02);
//
//        ImplicitCardiacMechanicsAssembler<2> implicit_assembler(&mesh,
//                                                                "ImplicitCardiacMech",
//                                                                &material_law);
//
//        std::vector<double> calcium_conc(implicit_assembler.GetTotalNumQuadPoints(), 1);
//        implicit_assembler.SetForcingQuantity(calcium_conc);
//
//        implicit_assembler.Solve(0,0.01,0.01);
//
//        // get the active tensions
//        std::vector<double> active_tension(implicit_assembler.GetTotalNumQuadPoints());
//        for(unsigned i=0; i<active_tension.size(); i++)
//        {
//            active_tension[i] = implicit_assembler.mCellMechSystems[i].GetActiveTension();
//        }
//
//        CardiacMechanicsAssembler<2> explicit_assembler(&mesh,"",&material_law);
//        explicit_assembler.SetForcingQuantity(active_tension);
//
//        // overwrite the current solution with what should be the true solution
//        // from the implicit solve
//        explicit_assembler.mCurrentSolution = implicit_assembler.mCurrentSolution;
//        // need to call this otherwise the assembler thinks it is the first time
//        // so it guesses the ZeroDeformation solution for the pressure
//        explicit_assembler.mADeformedHasBeenSolved=true;
//
//        explicit_assembler.Solve(0,0.01,0.01);
//
//        // check there were no newton iterations needed, ie the solution of the implicit
//        // method with these active tensions is the same as the solution of the explicit
//        TS_ASSERT_EQUALS(explicit_assembler.GetNumNewtonIterations(),0u);
//
//        //// NOTES:
//        // if we don't provide the implicit solution to the explicit assembler, and let the
//        // explicit assembler find it's own solution, the explicit solution can be compared
//        // to the implicit solution with the below code. With
//        // FiniteElasticityAssembler::NEWTON_ABS_TOL = 1e-13 and AbstractDealiiAssembler->gmres
//        // tol = 1e-6 the results are not that close visually. 2e-15 and 1e-8 is much
//        // better.
//        //
//        // std::vector<Vector<double> >& r_impl_solution = implicit_assembler.rGetDeformedPosition();
//        // std::vector<Vector<double> >& r_expl_solution = explicit_assembler.rGetDeformedPosition();
//        // for(unsigned i=0; i<r_impl_solution[0].size(); i++)
//        // {
//        //     for(unsigned j=0; j<2; j++)
//        //     {
//        //         double tol = fabs(r_impl_solution[j](i)/100);
//        //         TS_ASSERT_DELTA(r_impl_solution[j](i),r_expl_solution[j](i),tol);
//        //     }
//        // }
//    }


};

#endif /*TESTIMPLICITCARDIACMECHANICSASSEMBLER2_HPP_*/
