#ifndef TESTFINITEELASTICITYASSEMBLER_HPP_
#define TESTFINITEELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "FiniteElasticityAssembler.cpp"

#include "TriangulationVertexIterator.hpp"

#include "DofVertexIterator.hpp"

#define DIMENSION 2

class TestFiniteElasticityAssembler : public CxxTest::TestSuite
{
public:
    void test2dProblemOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(0) = 6.0;
    
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);


        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/simple2d");
        finite_elasticity.Solve();


        Vector<double>& solution = finite_elasticity.GetSolutionVector();
        DoFHandler<2>& dof_handler = finite_elasticity.GetDofHandler();


        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
        
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<2> old_posn = vertex_iter.GetVertex();
            
            Point<2> new_posn;
            new_posn(0) = old_posn(0)+solution(vertex_iter.GetDof(0));
            new_posn(1) = old_posn(1)+solution(vertex_iter.GetDof(1));
            
            // todo: TEST THESE!!

            std::cout << vertex_index << " " << old_posn(0) << " " << old_posn(1)
                                      << " " << new_posn(0) << " " << new_posn(1) << "\n";
            vertex_iter.Next();
        }
    }
    
    
    void test3dProblemOnSquare() throw(Exception)
    {
        Vector<double> body_force(3);
        body_force(0) = 6.0;
    
        MooneyRivlinMaterialLaw<3> mooney_rivlin_law(2.0,2.0);


        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(2);

        FiniteElasticityAssembler<3> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/simple3d");
        finite_elasticity.Solve();


        Vector<double>& solution = finite_elasticity.GetSolutionVector();
        DoFHandler<3>& dof_handler = finite_elasticity.GetDofHandler();


        DofVertexIterator<3> vertex_iter(&mesh, &dof_handler);
        
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<3> old_posn = vertex_iter.GetVertex();
            
            Point<3> new_posn;
            new_posn(0) = old_posn(0)+solution(vertex_iter.GetDof(0));
            new_posn(1) = old_posn(1)+solution(vertex_iter.GetDof(1));
            new_posn(2) = old_posn(2)+solution(vertex_iter.GetDof(2));
            
            // todo: TEST THESE!!

            std::cout << vertex_index << " " << old_posn(0) << " " << old_posn(1) << " " << old_posn(2) 
                                      << " " << new_posn(0) << " " << new_posn(1) << " " << new_posn(2)
                                      << "\n";
            vertex_iter.Next();
        }
    }
};
#endif /*TESTFINITEELASTICITYASSEMBLER_HPP_*/
