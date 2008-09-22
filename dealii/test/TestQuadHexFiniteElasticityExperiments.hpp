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


class ModelProblem2
{
public:
    static const double ALPHA = 0.1;
    static const double C1 = 0.02;
    static const double C2 = 0.1;
};


Vector<double> TractionForMixedShears(const Point<2>& X)
{
    assert(X[0]==1 || X[1]==0 || X[1]==1);
    
    Vector<double> traction(2);
    if(X[1]==1)
    {
        traction(0) = -2*ModelProblem2::ALPHA * ModelProblem2::C1;
        traction(1) =  0;
    }
    else if (X[1]==0)
    {
        traction(0) = 2*ModelProblem2::ALPHA*ModelProblem2::C1;
        traction(1) = 0;
    }
    else if (X[0]==1)
    {
        traction(0) = 0;
        traction(1) = -2*ModelProblem2::ALPHA*ModelProblem2::C1; // also equal to -2*c2*beta;
    }
    else
    {
        NEVER_REACHED;
    }
    return traction;
}


// class for non-zero dirichlet boundary conditions
class DirichletFunctionForModelProblem1 : public Function<2>
{
private:
    double mLambda;
public:
    DirichletFunctionForModelProblem1(double lambda)
     : Function<2>(3),
       mLambda(lambda)
    {
    }
    
    void vector_value(const Point<2>& p, Vector<double> &values) const
    {
        assert(values.size()==3);
        values(0) = 0.0;
        values(1) = (mLambda-1)*p(1);
        values(2) = 0.0;
    }
};
    

class DealiiModelProblem3
{
    static const double d = 0.05;
    static const double k = 10;
public:
    static const double c1 = 0.1;
    static double Mu(double X)
    {
        double t = tanh(k*(X-0.5));
        return (d/2)*(t+1);
    }
    static double MuDash(double X)
    {
        double t = tanh(k*(X-0.5));
        return (d*k/2)*(1-t*t);
    }
    static double MuDashDash(double X)
    {
        double t = tanh(k*(X-0.5));
        return (d*k*k/2)*(-2*t + 2*t*t*t);
    }
    static double MuDashDashDash(double X)
    {
        double t = tanh(k*(X-0.5));
        return (d*k*k*k/2)*(-2 + 8*t*t - 6*t*t*t*t);
    }

    static Vector<double> GetBodyForce(const Point<2>& X)
    {
        assert(X[0]>=0 && X[0]<=1 && X[1]>=0 && X[1]<=1);
    
        Vector<double> body_force(2);
    
        double lam = 1-MuDash(X[0]);
        double ddmu = MuDashDash(X[0]);
        double dddmu = MuDashDashDash(X[0]);
    
        body_force(0) =  2*c1*ddmu;
        body_force(1) =  2*c1*(-X[1]*dddmu/(lam*lam) - 2*ddmu*ddmu*X[1]/(lam*lam*lam));
        return body_force;
    }

    static Vector<double> GetTraction(const Point<2>& X)
    {
        Vector<double> traction(2);
    
        double lam = 1-MuDash(X[0]);
        if(X[0]==1)
        {
            traction(0) =  2*c1*(lam - 1.0/lam);
            traction(1) =  2*c1*X[1]*MuDashDash(1)/(lam*lam);
        }
        else if(X[1]==0)
        {
            traction(1) =  2*c1*(lam - 1.0/lam);
        }
        else if(X[1]==1)
        {
            traction(0) =  2*c1*MuDashDash(X[0])/(lam*lam); //*Y where Y=1
            traction(1) =  2*c1*(-lam + 1.0/lam);
        }
        else
        {
            NEVER_REACHED;
        }
        return traction;
    }
};

// class for non-zero dirichlet boundary conditions
class ModelProblem3DirichletValue : public Function<2>
{
public:
    ModelProblem3DirichletValue()
     : Function<2>(3)
    {
    }
    
    // Note: here we provide displacement..
    void vector_value(const Point<2>& p, Vector<double> &values) const
    {
        assert(values.size()==3);
        double exact_u = -DealiiModelProblem3::Mu(p[0]);
        double exact_v = p[1] - p[1]/(1-DealiiModelProblem3::MuDash(p[0]));

        values(0) = exact_u;
        values(1) = exact_v;
        values(2) = 0.0;
    }
};


class TestQuadHexFiniteElasticityExperiments : public CxxTest::TestSuite
{
private:
    unsigned FindWatchedPointIndex(Point<2> watchedPoint, std::vector<Vector<double> >& rUndeformedPosition)
    {
        // find the index corresponding to given point
        for(unsigned i=0; i<rUndeformedPosition[0].size(); i++)
        {
            double X = rUndeformedPosition[0](i);
            double Y = rUndeformedPosition[1](i);
            
            if(    (fabs(X-watchedPoint[0]) < 1e-9)
                && (fabs(Y-watchedPoint[1]) < 1e-9) )
            {
                return i;
            }
        }

        LOG_AND_COUT(1,"Cannot find node matching watched point");
        TS_FAIL("");
        return 0;
    }


    Point<2> RunModelProblem1(unsigned numRefinements, 
                              Point<2> watchedPoint,
                              bool testing = false)
    {
        const double LAMBDA = 0.75;
        double c1 = 0.02;
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw<2> law(c1);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);


        if(numRefinements>0)
        {
            mesh.refine_global(numRefinements);
        }
        ////////////////////////////////////////////////////////
        // define dirichlet (X=0) and Neumann (X=1) boundaries
        ////////////////////////////////////////////////////////
        unsigned component = 0;
        double value = 1.0;
        Triangulation<2>::cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    double component_val = element_iter->face(face_index)->center()(component);
                    if (fabs(component_val)<1e-4)
                    {
                        // X=0, label as dirichlet boundary
                        element_iter->face(face_index)->set_boundary_indicator(DIRICHLET_BOUNDARY);
                    }
                    else if (fabs(component_val - value)<1e-4)
                    {
                        // X=1, label as neumann boundary
                        element_iter->face(face_index)->set_boundary_indicator(NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }


  
        std::stringstream dir;
        if(!testing)
        {
            dir << "DealiiModelProblem1/" << numRefinements;
        }
        else
        {
            dir << "TestVerifySolveDealiiModelProblem1";
        }
 
        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &law,
                                                       body_force,
                                                       1.0,
                                                       dir.str());

        // apply traction
        Vector<double> traction(2);
        traction(0) = 2*c1*(pow(LAMBDA,-1) - LAMBDA*LAMBDA*LAMBDA);
        traction(1) = 0;
        finite_elasticity.SetConstantSurfaceTraction(traction);

        //////////////////////////////////////////////////
        // define non-zero dirichlet boundary conditions
        //////////////////////////////////////////////////
        std::map<unsigned,double> boundary_values;

        std::vector<bool> component_mask(2+1); // dim+1
        component_mask[0] = true;
        component_mask[1] = true;
        component_mask[2] = false;

        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 DIRICHLET_BOUNDARY,
                                                 DirichletFunctionForModelProblem1(LAMBDA),
                                                 boundary_values,
                                                 component_mask);

        assert(!boundary_values.empty());
        finite_elasticity.SetBoundaryValues(boundary_values);


        // solve
        finite_elasticity.StaticSolve();
                
        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        
        if(testing)
        {
            for(unsigned i=0; i < r_deformed_position[0].size(); i++)
            {
                double X = r_undeformed_position[0](i);
                double Y = r_undeformed_position[1](i);
                double exact_x = (1.0/LAMBDA)*X;
                double exact_y = LAMBDA*Y;
                
                double tol = 1e-4;
                if(fabs(X)<1e-6)
                {
                    tol = 1e-9;
                }
                
                TS_ASSERT_DELTA( r_deformed_position[0](i), exact_x, tol );
                TS_ASSERT_DELTA( r_deformed_position[1](i), exact_y, tol );
            }
            
            // check the final pressure
            Vector<double>& full_solution = finite_elasticity.rGetCurrentSolution();
            DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
    
            while (!vertex_iter.ReachedEnd())
            {
                // get the pressure at this node
                double pressure = full_solution(vertex_iter.GetDof(2));
                TS_ASSERT_DELTA(pressure, 2*c1*LAMBDA*LAMBDA, 1e-4);
                vertex_iter.Next();
            }
        }

        unsigned watched_point_index = FindWatchedPointIndex(watchedPoint, r_undeformed_position);
        
        Point<2> error;
        double X = r_undeformed_position[0](watched_point_index);
        double Y = r_undeformed_position[1](watched_point_index);
        double exact_x = (1.0/LAMBDA)*X;
        double exact_y = LAMBDA*Y;

        error[0] = r_deformed_position[0](watched_point_index) - exact_x;
        error[1] = r_deformed_position[1](watched_point_index) - exact_y;

        return error;
    }
    
    Point<2> RunModelProblem2(unsigned numRefinements, 
                              Point<2> watchedPoint,
                              bool testing = false)
    {
        double beta = (ModelProblem2::C1/ModelProblem2::C2)*ModelProblem2::ALPHA;

        Vector<double> body_force(2);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        if(numRefinements>0)
        {
            mesh.refine_global(numRefinements);
        }

        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        MooneyRivlinMaterialLaw<2> law1(ModelProblem2::C1);
        MooneyRivlinMaterialLaw<2> law2(ModelProblem2::C2);

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

        std::stringstream dir;
        if(!testing)
        {
            dir << "DealiiModelProblem2/" << numRefinements;
        }
        else
        {
            dir << "TestVerifySolveDealiiModelProblem2";
        }
 
        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       NULL,
                                                       body_force,
                                                       1.0,
                                                       dir.str());

        finite_elasticity.SetMaterialLawsForHeterogeneousProblem(laws, material_ids);
        finite_elasticity.SetFunctionalTractionBoundaryCondition(TractionForMixedShears);

        // solve
        finite_elasticity.StaticSolve();
                
        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        
        if(testing)
        {
            for(unsigned i=0; i < r_deformed_position[0].size(); i++)
            {
                double X = r_undeformed_position[0](i);
                double Y = r_undeformed_position[1](i);
        
                double exact_y = X<0.5 ? Y - ModelProblem2::ALPHA*X : Y - beta*X + (beta-ModelProblem2::ALPHA)/2;
                 
                double tol = 3e-2;
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

        unsigned watched_point_index = FindWatchedPointIndex(watchedPoint, r_undeformed_position);
        
        Point<2> error;
        double X = r_undeformed_position[0](watched_point_index);
        double Y = r_undeformed_position[1](watched_point_index);
        double exact_x = X;
        double exact_y = X<0.5 ? Y - ModelProblem2::ALPHA*X : Y - beta*X + (beta-ModelProblem2::ALPHA)/2;

        error[0] = r_deformed_position[0](watched_point_index) - exact_x;
        error[1] = r_deformed_position[1](watched_point_index) - exact_y;

        return error;
    }
    
    

    Point<2> RunModelProblem3(unsigned numRefinements, 
                              Point<2> watchedPoint,
                              bool testing = false)
    {
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw<2> law(DealiiModelProblem3::c1);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
   
        if(numRefinements>0)
        {
            mesh.refine_global(numRefinements);
        }

        ////////////////////////////////////////////////////////
        // define dirichlet (X=0) and Neumann (X!=0) boundaries
       ////////////////////////////////////////////////////////
        Triangulation<2>::cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    double component_val = element_iter->face(face_index)->center()(0);
                    if (fabs(component_val)<1e-6)
                    {
                        // X=0, label as dirichlet boundary
                        element_iter->face(face_index)->set_boundary_indicator(DIRICHLET_BOUNDARY);
                    }
                    else
                    {
                        // X!=0, label as neumann boundary
                        element_iter->face(face_index)->set_boundary_indicator(NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }

        std::stringstream dir;
        if(!testing)
        {
            dir << "DealiiModelProblem3/" << numRefinements;
        }
        else
        {
            dir << "TestVerifySolveDealiiModelProblem3";
        }
        
        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &law,
                                                       body_force,
                                                       1.0,
                                                       dir.str());

        // define boundary values
        std::map<unsigned,double> boundary_values;
        std::vector<bool> component_mask(2+1); // dim+1
        component_mask[0] = true;
        component_mask[1] = true;
        component_mask[2] = false;

        ModelProblem3DirichletValue dirich_func;

        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 DIRICHLET_BOUNDARY,
                                                 dirich_func,
                                                 boundary_values,
                                                 component_mask);

        assert(!boundary_values.empty());        

        // set body force, traction and boundary values..
        finite_elasticity.SetFunctionalTractionBoundaryCondition(DealiiModelProblem3::GetTraction);
        finite_elasticity.SetFunctionalBodyForce(DealiiModelProblem3::GetBodyForce);
        finite_elasticity.SetBoundaryValues(boundary_values);

        // solve
        finite_elasticity.StaticSolve();

        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        
        if(testing)
        {
            for(unsigned i=0; i < r_deformed_position[0].size(); i++)
            {
                double X = r_undeformed_position[0](i);
                double Y = r_undeformed_position[1](i);
        
                double exact_x = X - DealiiModelProblem3::Mu(X);
                double exact_y = Y/(1-DealiiModelProblem3::MuDash(X));
                 
                TS_ASSERT_DELTA( r_deformed_position[0](i), exact_x, 1e-3 );
                TS_ASSERT_DELTA( r_deformed_position[1](i), exact_y, 2e-3 );
            }
    
            // check the final pressure
            Vector<double>& full_solution = finite_elasticity.rGetCurrentSolution();
            DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
            while (!vertex_iter.ReachedEnd())
            {
                double pressure = full_solution(vertex_iter.GetDof(2));
                TS_ASSERT_DELTA(pressure/(2*DealiiModelProblem3::c1), 1.0, 0.12);//tol = 0.04 works fine except at one node -  the middle?
                vertex_iter.Next();
            }
        }

        unsigned watched_point_index = FindWatchedPointIndex(watchedPoint, r_undeformed_position);
        
        Point<2> error;
        double X = r_undeformed_position[0](watched_point_index);
        double Y = r_undeformed_position[1](watched_point_index);
        double exact_x = X - DealiiModelProblem3::Mu(X);
        double exact_y = Y/(1-DealiiModelProblem3::MuDash(X));

        error[0] = r_deformed_position[0](watched_point_index) - exact_x;
        error[1] = r_deformed_position[1](watched_point_index) - exact_y;

        return error;

    }

public:
    void TestVerifyRunModelProblem1() throw(Exception)
    {
        Point<2> X;
        X[0] = 1.0;
        X[1] = 1.0;

        Point<2> error = RunModelProblem1(3,X,true);
        TS_ASSERT_DELTA(error[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(error[1], 0.0, 1e-3);
    }

    void TestVerifyRunModelProblem2() throw(Exception)
    {
        Point<2> X;
        X[0] = 0.0;
        X[1] = 1.0;

        Point<2> error = RunModelProblem2(3,X,true);
        TS_ASSERT_DELTA(error[0], 0.0, 1e-2);
        TS_ASSERT_DELTA(error[1], 0.0, 1e-2);
    }

    void TestVerifyRunModelProblem3() throw(Exception)
    {
        Point<2> X;
        X[0] = 0.5;
        X[1] = 1.0;

        Point<2> error = RunModelProblem3(4,X,true);
        TS_ASSERT_DELTA(error[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(error[1], 0.0, 1e-3);
    }




//    void dontTestConvergence() throw(Exception)
//    {
//        unsigned num_sims = 6;
//        
//        LogFile::Instance()->Set(1, "quad_hex_convergence");
//        LogFile::Instance()->SetPrecision(9);
//        
//        unsigned num_elem_in_each_dir = 1;
//        std::vector<unsigned> num_elems;
//        std::vector<c_vector<double,2> > results;
//        
//        for(unsigned i=0; i<num_sims; i++)
//        {
//            Vector<double> body_force(2);
//            body_force(0) = 0.06;
//            MooneyRivlinMaterialLaw<2> mooney_rivlin_law(0.02);
//    
//            Triangulation<2> mesh;
//            GridGenerator::hyper_cube(mesh, 0.0, 1.0);
//
//            if(i>0)
//            {
//                mesh.refine_global(i);
//            }
//
//            FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
//            
//            std::stringstream dir;
//            dir << "quad_hex_convergence/" << num_elem_in_each_dir;
//    
//            FiniteElasticityAssembler<2> finite_elasticity(&mesh,
//                                                           &mooney_rivlin_law,
//                                                           body_force,
//                                                           1.0,
//                                                           dir.str());
//            
//            double tol = 1e-4*finite_elasticity.CalculateResidualNorm();
//            if(tol > 1e-8)
//            {
//                tol = 1e-8;
//            }
//            else if(tol < 1e-12)
//            {
//                tol = 1e-12;
//            }
//
//            finite_elasticity.StaticSolve(true, tol);
//
//            std::vector<Vector<double> >& r_solution = finite_elasticity.rGetDeformedPosition();
//            std::vector<Vector<double> >& r_undef_position = finite_elasticity.rGetUndeformedPosition();
//
//            unsigned top_corner_index = 2;
//            TS_ASSERT_DELTA(r_undef_position[0](top_corner_index), 1.0, 1e-6);
//            TS_ASSERT_DELTA(r_undef_position[1](top_corner_index), 1.0, 1e-6);
//
//            LOG_AND_COUT(1,num_elem_in_each_dir << " " << r_solution[0](top_corner_index)<< " " << r_solution[1](top_corner_index));
//            num_elems.push_back(num_elem_in_each_dir);
//           
//            c_vector<double,2> result;
//            result(0) = r_solution[0](top_corner_index);
//            result(1) = r_solution[1](top_corner_index);
//            results.push_back(result);
//            
//            num_elem_in_each_dir*=2;
//        }
//
//        std::cout << "\n\n";
//        for(unsigned i=0; i<results.size(); i++)
//        {
//            std::cout << std::setprecision(9) << num_elems[i] << " " << results[i](0) << " " << results[i](1) << std::endl;
//        }
//    }
};
#endif /*TESTQUADHEXFINITEELASTICITYEXPERIMENTS_HPP_*/
