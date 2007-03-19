#ifndef TESTFINITEELASTICITYASSEMBLER_HPP_
#define TESTFINITEELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "FiniteElasticityAssembler.cpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

#include "MooneyRivlinMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "ExponentialMaterialLaw.hpp"

#include "FiniteElasticityTools.hpp"


// todos: proper test of answers


class TestFiniteElasticityAssembler : public CxxTest::TestSuite
{
public :
    void testExceptions() throw(Exception)
    {
        Vector<double> body_force(2);
        Vector<double> bad_body_force(3);
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);

        // should throw because the mesh has no surface elements set as the fixed boundary                
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssembler<2> bad_fe1(&mesh,&mooney_rivlin_law,body_force,1.0,""));

        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);

        // should throw because of the negative density
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssembler<2> bad_fe2(&mesh,&mooney_rivlin_law,body_force,-1.0,""));
        // should throw because the body force is the wrong size
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssembler<2> bad_fe3(&mesh,&mooney_rivlin_law,bad_body_force,1.0,""));

        // should be ok now
        TS_ASSERT_THROWS_NOTHING(FiniteElasticityAssembler<2> fe3(&mesh,&mooney_rivlin_law,body_force,1.0,""));

        std::vector<unsigned> material_ids;
        material_ids.push_back(0);

        std::vector<AbstractIncompressibleMaterialLaw<2>*> material_laws;
        material_laws.push_back(&mooney_rivlin_law);
        material_laws.push_back(&mooney_rivlin_law);
        
        FiniteElasticityAssembler<2> fe(&mesh,&mooney_rivlin_law,body_force,1.0,"");
        
        // should thrown because material_laws and material_ids are not the same size
        TS_ASSERT_THROWS_ANYTHING(fe.SetMaterialLawsForHeterogeneousProblem(material_laws,material_ids));

        // check for exception is mesh contains elements whose material id is not
        // equal to either of those passed in in material_ids
        material_ids.clear();
        material_ids.push_back(5);
        material_ids.push_back(6);
        TS_ASSERT_THROWS_ANYTHING(fe.SetMaterialLawsForHeterogeneousProblem(material_laws,material_ids));
    }
    
    
    void TestCompareJacobians() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(0) = 6.0;
    
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");                                    
        TS_ASSERT_THROWS_NOTHING( finite_elasticity.CompareJacobians() );
    }


    // just tests whether the method rGetUndeformedPosition() returns a data structure 
    // that consistent with the mesh..
    void TestGetUndeformedPosition() throw(Exception)
    {
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");


        // get undeformed position
        std::vector<Vector<double> >& undeformed_position = finite_elasticity.rGetUndeformedPosition(); 
        TS_ASSERT_EQUALS(undeformed_position.size(), 2);
        TS_ASSERT_EQUALS(undeformed_position[0].size(), mesh.n_vertices());
        TS_ASSERT_EQUALS(undeformed_position[1].size(), mesh.n_vertices());
        
        TriangulationVertexIterator<2> vertex_iter(&mesh);
        
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();

            TS_ASSERT_DELTA(undeformed_position[0](vertex_index), posn(0), 1e-12);
            TS_ASSERT_DELTA(undeformed_position[1](vertex_index), posn(1), 1e-12);
            
            vertex_iter.Next();
        } 
    }
    
    // A test where the solution should be zero displacement
    // It mainly tests that the initial guess was set up correctly to
    // the final correct solution, ie u=0, p=zero_strain_pressure (!=0)
    void TestWithZeroDisplacement() throw(Exception)   
    {
        Vector<double> body_force(2); //zero
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(3.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        
        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");
                                                         
        finite_elasticity.Solve();

        TS_ASSERT_EQUALS(finite_elasticity.GetNumNewtonIterations(), 0);

        // get undeformed position
        std::vector<Vector<double> >& undeformed_position 
            = finite_elasticity.rGetUndeformedPosition();

        // get deformed position
        std::vector<Vector<double> >& deformed_position 
            = finite_elasticity.rGetDeformedPosition();
        for(unsigned i=0; i<deformed_position[0].size(); i++)
        {
            TS_ASSERT_DELTA(undeformed_position[0](i), deformed_position[0](i), 1e-8);
            TS_ASSERT_DELTA(undeformed_position[1](i), deformed_position[1](i), 1e-8);
        }
    }


    void Test2dProblemOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(0) = 6.0;
    
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/simple2d");
                                                         
        finite_elasticity.Solve();

        // get deformed position
        std::vector<Vector<double> >& deformed_position = finite_elasticity.rGetDeformedPosition(); 
        TS_ASSERT_EQUALS(deformed_position.size(), 2);
        TS_ASSERT_EQUALS(deformed_position[0].size(), mesh.n_vertices());
        TS_ASSERT_EQUALS(deformed_position[1].size(), mesh.n_vertices());
        


        // also get the solution vector directly and check the deformed position
        // object was set up correctly...
        Vector<double>& solution = finite_elasticity.rGetSolutionVector();
        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();

        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
        
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<2> old_posn = vertex_iter.GetVertex();
            
            Point<2> new_posn;
            new_posn(0) = old_posn(0)+solution(vertex_iter.GetDof(0));
            new_posn(1) = old_posn(1)+solution(vertex_iter.GetDof(1));
            
            TS_ASSERT_DELTA(deformed_position[0](vertex_index), new_posn(0), 1e-12);
            TS_ASSERT_DELTA(deformed_position[1](vertex_index), new_posn(1), 1e-12);
            
            // todo: TEST THESE!!

            std::cout << vertex_index << " " << old_posn(0) << " " << old_posn(1)
                                      << " " << new_posn(0) << " " << new_posn(1) << "\n";
          
          
            // some hardcoded tests
            if(vertex_index==62)
            {
                TS_ASSERT_DELTA(new_posn(0),0.31559,1e-4);
                TS_ASSERT_DELTA(new_posn(1),0.77431,1e-4);
            }
            if(vertex_index==38)
            {
                TS_ASSERT_DELTA(new_posn(0),0.47755,1e-4);
                TS_ASSERT_DELTA(new_posn(1),0.87779,1e-4);
            }
            if(vertex_index==15)
            {
                TS_ASSERT_DELTA(new_posn(0),0.32141,1e-4);
                TS_ASSERT_DELTA(new_posn(1),0.86898,1e-4);
            }
            if(vertex_index==80)
            {
                TS_ASSERT_DELTA(new_posn(0),0.14713,1e-4);
                TS_ASSERT_DELTA(new_posn(1),0.79705,1e-4);
            }
            if(vertex_index==32)
            {
                TS_ASSERT_DELTA(new_posn(0),0.00000,1e-4);
                TS_ASSERT_DELTA(new_posn(1),0.87500,1e-4);
            }
                                      

            //// UPDATE THE NODE POSITIONS
            // GetVertex returns a reference to a Point<DIM>, so this changes the mesh
            // directly. Do this so the new volume can be computed
            vertex_iter.GetVertex()[0] = new_posn(0);         
            vertex_iter.GetVertex()[1] = new_posn(1);         
                                      
            vertex_iter.Next();
        }

        // compute the deformed volume 
        // NOTE: this aren't very accurate volumes, since we have lost the
        // positions of the extra nodes (those used with quadratic basis functions)
        // and the measure() function below must use linear interpolation. Hence
        // the high tolerances
        double deformed_volume = 0.0;
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        while(element_iter!=mesh.end())
        {
            double element_volume = element_iter->measure();
            TS_ASSERT_DELTA(element_volume, 1.0/mesh.n_active_cells(), 1e-2); 
            
            deformed_volume += element_volume;
            element_iter++;
        }
        
        TS_ASSERT_DELTA(deformed_volume, 1.0, 1e-2);
    }
};
#endif /*TESTFINITEELASTICITYASSEMBLER_HPP_*/
