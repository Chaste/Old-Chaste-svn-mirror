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

#ifndef TESTQUADHEXFINITEELASTICITYEXPERIMENTS_HPP_
#define TESTQUADHEXFINITEELASTICITYEXPERIMENTS_HPP_


#include <cxxtest/TestSuite.h>
#include "FiniteElasticityAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "FiniteElasticityTools.hpp"
#include "UblasCustomFunctions.hpp"

double MATERIAL_PARAM = 0.05;
double ALPHA = 0.2;

double SHEARS_ALPHA = 0.1;
double SHEARS_C1 = 0.01;

Vector<double> TractionForMixedShears(const Point<2>& X)
{
    assert(X[0]==1 || X[1]==0 || X[1]==1);
    
    Vector<double> traction(2);
    if(X[1]==1)
    {
        traction(0) = -2*SHEARS_ALPHA * SHEARS_C1;
        traction(1) =  0;
    }
    else if (X[1]==0)
    {
        traction(0) = 2*SHEARS_ALPHA*SHEARS_C1;
        traction(1) = 0;
    }
    else if (X[0]==1)
    {
        traction(0) = 0;
        traction(1) = -2*SHEARS_ALPHA*SHEARS_C1; // also equal to -2*c2*beta;
    }
    else
    {
        NEVER_REACHED;
    }
    return traction;
}

Vector<double> MyBodyForce(const Point<2>& X)
{
    assert(X[0]>=0 && X[0]<=1 && X[1]>=0 && X[1]<=1);

    Vector<double> body_force(2);
    double lam = 1+ALPHA*X[0];
    body_force(0) = -2*MATERIAL_PARAM * ALPHA;
    body_force(1) = -2*MATERIAL_PARAM * 2*ALPHA*ALPHA*X[1]/(lam*lam*lam);
    return body_force;
}

Vector<double> MyTraction(const Point<2>& X)
{
    Vector<double> traction(2);
    
    double lam = 1+ALPHA*X[0];
    if(X[0]==1)
    {
        traction(0) =  2*MATERIAL_PARAM * (lam - 1.0/lam);
        traction(1) = -2*MATERIAL_PARAM * X[1]*ALPHA/(lam*lam);
    }
    else if(X[1]==0)
    {
        traction(0) =  2*MATERIAL_PARAM * X[1]*ALPHA/(lam*lam);
        traction(1) = -2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else if(X[1]==1)
    {
        traction(0) = -2*MATERIAL_PARAM * X[1]*ALPHA/(lam*lam);
        traction(1) =  2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else
    {
        NEVER_REACHED;
    }
    return traction;
}


class TestQuadHexFiniteElasticityExperiments : public CxxTest::TestSuite
{
public:
    void TestWithMixedShears() throw(Exception)
    {
        double c1 = SHEARS_C1;
        double c2 = 0.02;
        double alpha = SHEARS_ALPHA;
        double beta = (c1/c2)*alpha;

        Vector<double> body_force(2);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(4);

        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);


        MooneyRivlinMaterialLaw<2> law1(c1);
        MooneyRivlinMaterialLaw<2> law2(c2);

        std::vector<AbstractIncompressibleMaterialLaw<2>*> laws;
        laws.push_back(&law1);
        laws.push_back(&law2);

        Triangulation<2>::cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            if (element_iter->center()[0] <= 0.5)
            {
                element_iter->set_material_id(5);
            }
            else
            {
                element_iter->set_material_id(6);
            }
            element_iter++;
        }

        std::vector<unsigned> material_ids;
        material_ids.push_back(5);
        material_ids.push_back(6);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       NULL,
                                                       body_force,
                                                       1.0,
                                                       "dealii_finite_elas/mixed_shears");

        finite_elasticity.SetMaterialLawsForHeterogeneousProblem(laws, material_ids);

        finite_elasticity.SetFunctionalTractionBoundaryCondition(TractionForMixedShears);
//        finite_elasticity.SetFunctionalBodyForce(MyBodyForce);


        // solve
        finite_elasticity.StaticSolve();
                
                
        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        
        for(unsigned i=0; i < r_deformed_position[0].size(); i++)
        {
            double X = r_undeformed_position[0](i);
            double Y = r_undeformed_position[1](i);
    
            double exact_y = X<0.5 ? Y - alpha*X : Y - beta*X + (beta-alpha)/2;
             
            double tol = 1e-2;
            if(fabs(X)<1e-6)
            {
                tol = 1e-9;
            }
            
            TS_ASSERT_DELTA( r_deformed_position[0](i), X, tol );
            TS_ASSERT_DELTA( r_deformed_position[1](i), exact_y, tol );
        }
        
        // don't check the final pressure
        //  if quad-hex is like quad-tet they won't be anywhere near...
    }


    void TestWithFunctionalData() throw(Exception)
    {
        Vector<double> body_force(2);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(4);

        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        MooneyRivlinMaterialLaw<2> law(MATERIAL_PARAM);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &law,
                                                       body_force,
                                                       1.0,
                                                       "dealii_finite_elas/functional");

        finite_elasticity.SetFunctionalTractionBoundaryCondition(MyTraction);
        finite_elasticity.SetFunctionalBodyForce(MyBodyForce);

        // solve
        finite_elasticity.StaticSolve();
                
                
        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        
        for(unsigned i=0; i < r_deformed_position[0].size(); i++)
        {
            double X = r_undeformed_position[0](i);
            double Y = r_undeformed_position[1](i);
    
            double exact_x = X + 0.5*ALPHA*X*X;
            double exact_y = Y/(1+ALPHA*X);
             
            TS_ASSERT_DELTA( r_deformed_position[0](i), exact_x, 1e-3 );
            TS_ASSERT_DELTA( r_deformed_position[1](i), exact_y, 1e-3 );
        }

        // check the final pressure
        Vector<double>& full_solution = finite_elasticity.rGetCurrentSolution();
        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();
        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
        while (!vertex_iter.ReachedEnd())
        {
            double pressure = full_solution(vertex_iter.GetDof(2));
            TS_ASSERT_DELTA(pressure/(2*MATERIAL_PARAM), 1.0, 1.3e-3);
            vertex_iter.Next();
        }
    }



    void dontTestConvergence() throw(Exception)
    {
        unsigned num_sims = 6;
        
        LogFile::Instance()->Set(1, "quad_hex_convergence");
        LogFile::Instance()->SetPrecision(9);
        
        unsigned num_elem_in_each_dir = 1;
        std::vector<unsigned> num_elems;
        std::vector<c_vector<double,2> > results;
        
        for(unsigned i=0; i<num_sims; i++)
        {
            Vector<double> body_force(2);
            body_force(0) = 0.06;
            MooneyRivlinMaterialLaw<2> mooney_rivlin_law(0.02);
    
            Triangulation<2> mesh;
            GridGenerator::hyper_cube(mesh, 0.0, 1.0);

            if(i>0)
            {
                mesh.refine_global(i);
            }

            FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
            
            std::stringstream dir;
            dir << "quad_hex_convergence/" << num_elem_in_each_dir;
    
            FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                           &mooney_rivlin_law,
                                                           body_force,
                                                           1.0,
                                                           dir.str());
            
            double tol = 1e-4*finite_elasticity.CalculateResidualNorm();
            if(tol > 1e-8)
            {
                tol = 1e-8;
            }
            else if(tol < 1e-12)
            {
                tol = 1e-12;
            }

            finite_elasticity.StaticSolve(true, tol);

            std::vector<Vector<double> >& r_solution = finite_elasticity.rGetDeformedPosition();
            std::vector<Vector<double> >& r_undef_position = finite_elasticity.rGetUndeformedPosition();

            unsigned top_corner_index = 2;
            TS_ASSERT_DELTA(r_undef_position[0](top_corner_index), 1.0, 1e-6);
            TS_ASSERT_DELTA(r_undef_position[1](top_corner_index), 1.0, 1e-6);

            LOG_AND_COUT(1,num_elem_in_each_dir << " " << r_solution[0](top_corner_index)<< " " << r_solution[1](top_corner_index));
            num_elems.push_back(num_elem_in_each_dir);
           
            c_vector<double,2> result;
            result(0) = r_solution[0](top_corner_index);
            result(1) = r_solution[1](top_corner_index);
            results.push_back(result);
            
            num_elem_in_each_dir*=2;
        }

        std::cout << "\n\n";
        for(unsigned i=0; i<results.size(); i++)
        {
            std::cout << std::setprecision(9) << num_elems[i] << " " << results[i](0) << " " << results[i](1) << std::endl;
        }
    }
};
#endif /*TESTQUADHEXFINITEELASTICITYEXPERIMENTS_HPP_*/
