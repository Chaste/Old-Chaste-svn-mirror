/*

Copyright (C) University of Oxford, 2005-2009

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
#include "MooneyRivlinMaterialLaw2.hpp"
#include "PolynomialMaterialLaw3d2.hpp"
#include "ExponentialMaterialLaw2.hpp"
#include "FiniteElasticityTools.hpp"
#include "UblasCustomFunctions.hpp"




Vector<double> ModelProblem4GetTraction(const Point<2>& X)
{
    assert(X[0]==0 || X[0]==1 || X[1]==1);
    
    Vector<double> traction(2);
    traction(0) = 0;
    traction(1) = 0;

    if(X[1]==1)
    {
        if( (X[0]>=0.25) && (X[0]<=0.75) )
        {
            traction(1) = -0.06;
        }
    }

    return traction;
}

class ModelProblem2
{
public:
    static const double ALPHA = 0.1;
    static const double C1 = 0.1;
    static const double C2 = 0.02;
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




class DealiiModelProblem3
{
    static const double a = 0.2;
    static const double d = 0.05;
    static const double k = 10;
public:
    static const double c1 = 0.1;
    static double Mu(double X)
    {
//        return -a*0.5*X*X;
        double t = tanh(k*(X-0.5));
        return (d/2)*(t+1);
    }
    static double MuDash(double X)
    {
//        return -a*X;
        double t = tanh(k*(X-0.5));
        return (d*k/2)*(1-t*t);
    }
    static double MuDashDash(double X)
    {
//        return -a;
        double t = tanh(k*(X-0.5));
        return (d*k*k/2)*(-2*t + 2*t*t*t);
    }
    static double MuDashDashDash(double X)
    {
//        return 0.0;
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



class FiniteElasticityWithWatchedInfo : public FiniteElasticityAssembler<2>
{
private:
    Point<2> mWatchedPoint;
    unsigned mNodeIndexOfWatchedPoint;
    
    void PostNewtonStep(unsigned counter, double normResidual)
    {
        assert(mNodeIndexOfWatchedPoint!=UINT_MAX);

        //// MP1-3
        Point<2> error;
        double X = mWatchedPoint[0];
        double Y = mWatchedPoint[1];

        //// MP1 or 3
        //double exact_x = X - DealiiModelProblem3::Mu(X);
        //double exact_y = Y/(1-DealiiModelProblem3::MuDash(X));

        //// MP2
        double exact_x = X;
        double exact_y = Y - 0.5*X + (0.5-0.1)/2;

        //// MP1-3
        std::vector<Vector<double> >& r_deformed_position = this->rGetDeformedPosition();
        error(0) = r_deformed_position[0](mNodeIndexOfWatchedPoint) - exact_x;
        error(1) = r_deformed_position[1](mNodeIndexOfWatchedPoint) - exact_y;
        double norm2_err = sqrt(error(0)*error(0)+error(1)*error(1));
        LOG_AND_COUT(1,"\tNewton step " << counter << ": " << normResidual << " " << norm2_err);

        //// MP4
        //std::vector<Vector<double> >& r_deformed_position = this->rGetDeformedPosition();
        //LOG_AND_COUT(1,"\tNewton step " << counter << ": " << normResidual << " " << r_deformed_position[0](mNodeIndexOfWatchedPoint) << " " << r_deformed_position[1](mNodeIndexOfWatchedPoint));
        //CalculateForceOnTopSurface(counter);
    }
        

public:
    FiniteElasticityWithWatchedInfo(Triangulation<2>* pMesh,
                                    AbstractIncompressibleMaterialLaw2<2>*  pMaterialLaw,
                                    Vector<double> bodyForce,
                                    double density,
                                    std::string outputDirectory)
        : FiniteElasticityAssembler<2>(pMesh,pMaterialLaw,bodyForce,density,outputDirectory)
    {
        mNodeIndexOfWatchedPoint=UINT_MAX;
    }
    
    
    void AddWatchedLocation(Point<2> location, unsigned nodeIndex)
    {
        mWatchedPoint = location;
        mNodeIndexOfWatchedPoint = nodeIndex;
    }
    

    void CalculateForceOnTopSurface(unsigned counter)
    {
        // only write output if the flag mWriteOutput has been set
        if (!this->mWriteOutput)
        {
            return;
        }
    
        // create stresses file
        std::stringstream ss;
        ss << this->mOutputDirectoryFullPath << "/top_forces_" << counter << ".txt";
        std::string stress_filename = ss.str();
        std::ofstream stress_output(stress_filename.c_str());
    
        static QSimpson<2>   quadrature_formula;
        const unsigned n_q_points = quadrature_formula.n_quadrature_points;
    
        FEValues<2> fe_values(mFeSystem, quadrature_formula,
                                UpdateFlags(update_values    |
                                            update_gradients |
                                            update_q_points  |     // needed for interpolating u and u' on the quad point
                                            update_JxW_values));
    
        std::vector< Vector<double> >                  local_solution_values(n_q_points);
        std::vector< std::vector< Tensor<1,2> > >    local_solution_gradients(n_q_points);
    
        Tensor<2,2> identity;
        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            local_solution_values[q_point].reinit(2+1);
            local_solution_gradients[q_point].resize(2+1);
        }
    
        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                identity[i][j] = i==j ? 1.0 : 0.0;
            }
        }
    
        DoFHandler<2>::active_cell_iterator  element_iter = this->mDofHandler.begin_active();
    
        while (element_iter!=this->mDofHandler.end())
        {
            fe_values.reinit(element_iter); // compute fe values for this element
            fe_values.get_function_values(this->mCurrentSolution, local_solution_values);
            fe_values.get_function_grads(this->mCurrentSolution, local_solution_gradients);
    
            AbstractIncompressibleMaterialLaw2<2>* p_material_law = GetMaterialLawForElement(element_iter);
    
            for (unsigned q_point=0; q_point<n_q_points; q_point++)
            {
                const Point<2>& X = fe_values.quadrature_point(q_point);
                if(X[1]==1)
                {
                    unsigned N = (unsigned)(round(sqrt(this->mpMesh->n_active_cells())));
    
                    bool bdy = false;
                    double h = 1.0/N;
                    for(unsigned k=0; k<=N; k++)
                    {
                        if( fabs(X[0] - k*h)<1e-6)
                        {
                            bdy = true;
                        }
                    }
                    
                    if(!bdy)
                    {                    
                        const std::vector< Tensor<1,2> >& grad_u_p = local_solution_gradients[q_point];
        
                        Vector<double> x;
                        x.reinit(2);
                        for(unsigned i=0; i<2; i++)
                        {
                            x(i) = local_solution_values[q_point](i);
                        }
        
                        double p = local_solution_values[q_point](PRESSURE_COMPONENT_INDEX);
        
                        static Tensor<2,2> F;
                        static Tensor<2,2> S;
        
                        for (unsigned i=0; i<2; i++)
                        {
                            for (unsigned j=0; j<2; j++)
                            {
                                F[i][j] = identity[i][j] + grad_u_p[i][j];
                            }
                        }
        
                        p_material_law->Compute1stPiolaKirchoffStress(F,p,S);
                        
                        stress_output << X[0] << " " << X[1] << " " << S[0][1] << " " << S[1][1] << "\n";
                    }
                }
            }
    
            element_iter++;
        }
        stress_output.close();
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

        MooneyRivlinMaterialLaw2<2> law1(ModelProblem2::C1);
        MooneyRivlinMaterialLaw2<2> law2(ModelProblem2::C2);

        std::vector<AbstractIncompressibleMaterialLaw2<2>*> laws;
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
 
        FiniteElasticityWithWatchedInfo finite_elasticity(&mesh,
                                                          NULL,
                                                          body_force,
                                                          1.0,
                                                          dir.str());

        finite_elasticity.SetMaterialLawsForHeterogeneousProblem(laws, material_ids);
        finite_elasticity.SetFunctionalTractionBoundaryCondition(TractionForMixedShears);
        
        finite_elasticity.UseDirectSolver();

        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        unsigned watched_point_index = FindWatchedPointIndex(watchedPoint, r_undeformed_position);
        finite_elasticity.AddWatchedLocation(watchedPoint,watched_point_index);


        // solve
        try
        {
            finite_elasticity.StaticSolve(true, 1e-18);
        }
        catch(Exception& e)
        {
            LOG_AND_COUT(1,e.GetMessage());
        }
                
        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        
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
        MooneyRivlinMaterialLaw2<2> law(DealiiModelProblem3::c1);

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
        
        FiniteElasticityWithWatchedInfo finite_elasticity(&mesh,
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

        finite_elasticity.UseDirectSolver();

        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        unsigned watched_point_index = FindWatchedPointIndex(watchedPoint, r_undeformed_position);
        finite_elasticity.AddWatchedLocation(watchedPoint,watched_point_index);


        // solve
        try
        {
            finite_elasticity.StaticSolve(true, 1e-19);
        }
        catch(Exception& e)
        {
            LOG_AND_COUT(1,e.GetMessage());
        }

        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        
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

        
        Point<2> error;
        double X = r_undeformed_position[0](watched_point_index);
        double Y = r_undeformed_position[1](watched_point_index);
        double exact_x = X - DealiiModelProblem3::Mu(X);
        double exact_y = Y/(1-DealiiModelProblem3::MuDash(X));

        error[0] = r_deformed_position[0](watched_point_index) - exact_x;
        error[1] = r_deformed_position[1](watched_point_index) - exact_y;

        return error;

    }







    Point<2> RunModelProblem4(unsigned numRefinements, 
                              Point<2> watchedPoint,
                              bool testing = false)
    {
        Vector<double> body_force(2);
        body_force(0) = 0;
        body_force(1) = 0;
        MooneyRivlinMaterialLaw2<2> law(0.02);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
   
        if(numRefinements>0)
        {
            mesh.refine_global(numRefinements);
        }

        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 1, 0.0);
        
        std::stringstream dir;
        if(!testing)
        {
            dir << "DealiiModelProblem4/" << numRefinements;
        }
        else
        {
            dir << "TestVerifySolveDealiiModelProblem4";
        }
        
        FiniteElasticityWithWatchedInfo finite_elasticity(&mesh,
                                                          &law,
                                                          body_force,
                                                          1.0,
                                                          dir.str());


        // set traction
        finite_elasticity.SetFunctionalTractionBoundaryCondition(ModelProblem4GetTraction);

        finite_elasticity.UseDirectSolver();

        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        unsigned watched_point_index = FindWatchedPointIndex(watchedPoint, r_undeformed_position);
        finite_elasticity.AddWatchedLocation(watchedPoint,watched_point_index);

        // solve
        try
        {
            finite_elasticity.StaticSolve(true, 1e-19);
        }
        catch(Exception& e)
        {
            LOG_AND_COUT(1,e.GetMessage());
        }

        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        
        if(testing)
        {
            assert(0);
        }

        
        Point<2> ret;
        ret[0] = r_deformed_position[0](watched_point_index);
        ret[1] = r_deformed_position[1](watched_point_index);

        return ret;
    }


public:

    void DONT_TestVerifyRun__NEW__ModelProblem1() throw(Exception)
    {
        assert(DealiiModelProblem3::Mu(1)==-0.1);

        Point<2> X;
        X[0] = 1.0;
        X[1] = 1.0;

        Point<2> error = RunModelProblem3(4,X,true);
        TS_ASSERT_DELTA(error[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(error[1], 0.0, 1e-3);
    }


    void dont__TestVerifyRunModelProblem2() throw(Exception)
    {
        Point<2> X;
        X[0] = 0.0;
        X[1] = 1.0;

        Point<2> error = RunModelProblem2(3,X,true);
        TS_ASSERT_DELTA(error[0], 0.0, 1e-2);
        TS_ASSERT_DELTA(error[1], 0.0, 1e-2);
    }

    void DONT_TestVerifyRunModelProblem3() throw(Exception)
    {
        // check ModelProblem3 is in correct mode
        assert(fabs(DealiiModelProblem3::Mu(1)-0.05)<1e-2);

        Point<2> X;
        X[0] = 0.5;
        X[1] = 1.0;

        Point<2> error = RunModelProblem3(4,X,true);
        TS_ASSERT_DELTA(error[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(error[1], 0.0, 1e-3);
    }

    void DONT_TestRun_New_ModelProblem1() throw(Exception)
    {
        assert(DealiiModelProblem3::Mu(1)==-0.1);
        
        LogFile::Instance()->Set(1,"MechComp","dealii_prob1.txt");
        LogFile::Instance()->SetPrecision(10);

        Point<2> X;
        X[0] = 1.0; // 1.0 not 0.5
        X[1] = 1.0;

        for(unsigned i=1; i<=7; i++)
        {
            unsigned N = (unsigned)pow(2,i);
            unsigned num_nodes = (unsigned)((2*N-1)*(2*N-1));
            unsigned num_vertices = 2*(N+1)*(N+1);
            unsigned num_dofs = 2*num_nodes+num_vertices; 

            Point<2> error = RunModelProblem3(i,X);
    
            double norm_err = sqrt(error[0]*error[0]+error[1]*error[1]);
        
            LOG_AND_COUT(1,N<<" "<<num_nodes<<" "<<num_dofs<<" "<<error[0]<<" "<<error[1]<<" "<<norm_err);
        }
        LogFile::Close();
    }



    void TestRunModelProblem2() throw(Exception)
    {
        /// ***to run this properly, make sure PostNewtonStep() is doing the right thing***

        LogFile::Instance()->Set(1,"MechComp","dealii_prob2.txt");
        LogFile::Instance()->SetPrecision(10);

        Point<2> X;
        X[0] = 1.0;
        X[1] = 1.0;

        for(unsigned i=1; i<=4; i++)
        {
            unsigned N = (unsigned)pow(2,i);
            unsigned num_nodes = (unsigned)((2*N-1)*(2*N-1));
            unsigned num_vertices = 2*(N+1)*(N+1);
            unsigned num_dofs = 2*num_nodes+num_vertices; 

            Point<2> error = RunModelProblem2(i,X);
            double norm_err = sqrt(error[0]*error[0]+error[1]*error[1]);
        
            LOG_AND_COUT(1,N<<" "<<num_nodes<<" "<<num_dofs<<" "<<error[0]<<" "<<error[1]<<" "<<norm_err);
        }
        LogFile::Close();
    }

    void DONT_TestRunModelProblem3() throw(Exception)
    {
        /// ***to run this properly, make sure PostNewtonStep() is doing the right thing***

        // check ModelProblem3 is in correct mode
        assert(fabs(DealiiModelProblem3::Mu(1)-0.05)<1e-2);
        
        LogFile::Instance()->Set(1,"MechComp","dealii_prob3.txt");
        LogFile::Instance()->SetPrecision(10);

        Point<2> X;
        X[0] = 0.5;   // 0.5
        X[1] = 1.0;

        for(unsigned i=1; i<=4; i++)
        {
            unsigned N = (unsigned)pow(2,i);
            unsigned num_nodes = (unsigned)((2*N-1)*(2*N-1));
            unsigned num_vertices = 2*(N+1)*(N+1);
            unsigned num_dofs = 2*num_nodes+num_vertices; 

            Point<2> error = RunModelProblem3(i,X);
    
            double norm_err = sqrt(error[0]*error[0]+error[1]*error[1]);
        
            LOG_AND_COUT(1,N<<" "<<num_nodes<<" "<<num_dofs<<" "<<error[0]<<" "<<error[1]<<" "<<norm_err);
        }
        LogFile::Close();
    }
    
    void DONT_TestRunModelProblem4() throw(Exception)
    {
        /// ***to run this properly, make sure PostNewtonStep() is doing the right thing***
        
        LogFile::Instance()->Set(1,"MechComp","dealii_prob4.txt");
        LogFile::Instance()->SetPrecision(10);

        Point<2> X;
        X[0] = 0.5;   // 0.5
        X[1] = 1.0;

        for(unsigned i=2; i<=4; i++)
        {
            unsigned N = (unsigned)pow(2,i);
            unsigned num_nodes = (unsigned)((2*N-1)*(2*N-1));
            unsigned num_vertices = 2*(N+1)*(N+1);
            unsigned num_dofs = 2*num_nodes+num_vertices; 

            Point<2> def = RunModelProblem4(i,X);
    
            LOG_AND_COUT(1,N<<" "<<num_nodes<<" "<<num_dofs<<" "<<def[0]<<" "<<def[1]);
        }
        LogFile::Close();
    }

};
#endif /*TESTQUADHEXFINITEELASTICITYEXPERIMENTS_HPP_*/
