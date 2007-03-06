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


// todos: proper test of answers, fix compare numerical jacobian 


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
                                                       "finite_elas/simple2d");

        // to be fixed                                             
        //finite_elasticity.CompareJacobians();
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
    


    // Run same simulation on two meshes (one more refined than the other)
    // and test they agree on shared gridpoints
    void Test2dProblemOnSquareForConvergence() throw(Exception)
    {
        ////////////////////////////////////////////////
        // run 1: on a mesh which is refined 3 times..
        ////////////////////////////////////////////////
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

        Vector<double>& solution = finite_elasticity.GetSolutionVector();
        DoFHandler<2>& dof_handler = finite_elasticity.GetDofHandler();
        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
        
        unsigned num_vertices = mesh.n_vertices(); 
        std::vector<double> new_x(num_vertices);
        std::vector<double> new_y(num_vertices);

        for(unsigned i=0; i<num_vertices; i++)
        {
            new_x[i] = 0.0;
            new_y[i] = 0.0;
        }
        
        // store the new position
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<2> old_posn = vertex_iter.GetVertex();
            
            new_x[vertex_index] = old_posn(0)+solution(vertex_iter.GetDof(0));                          
            new_y[vertex_index] = old_posn(1)+solution(vertex_iter.GetDof(1));

            vertex_iter.Next();
        }


        //////////////////////////////////////////////////////////////
        // run 2: same problem, on a mesh which is refined 4 times..
        //////////////////////////////////////////////////////////////
        Triangulation<2> mesh_refined;
        GridGenerator::hyper_cube(mesh_refined, 0.0, 1.0); 
        mesh_refined.refine_global(4);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh_refined, 0);

        FiniteElasticityAssembler<2> finite_elasticity_ref(&mesh_refined,
                                                           &mooney_rivlin_law,
                                                           body_force,
                                                           1.0,
                                                           "finite_elas/simple2d");
        finite_elasticity_ref.Solve();



        //////////////////////////////////////////////////////////////
        // compare the solution with that on the previous mesh
        //////////////////////////////////////////////////////////////
        Vector<double>& solution_ref = finite_elasticity_ref.GetSolutionVector();
        DoFHandler<2>& dof_handler_ref = finite_elasticity_ref.GetDofHandler();

        DofVertexIterator<2> vertex_iter_ref(&mesh_refined, &dof_handler_ref);
        while(!vertex_iter_ref.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter_ref.GetVertexGlobalIndex();
            
            // if vertex_index < num_vertices (in the first mesh), this
            // node is in both meshes, so compare the results
            if(vertex_index < num_vertices)
            {
                Point<2> old_posn = vertex_iter_ref.GetVertex();
            
                Point<2> new_posn;
                new_posn(0) = old_posn(0)+solution_ref(vertex_iter_ref.GetDof(0));
                new_posn(1) = old_posn(1)+solution_ref(vertex_iter_ref.GetDof(1));
                
                TS_ASSERT_DELTA( new_posn(0), new_x[vertex_index], 1e-2 );
                TS_ASSERT_DELTA( new_posn(1), new_y[vertex_index], 1e-2 );
            }
            
            vertex_iter_ref.Next();
        }
        
        // check nothing has changed
        TS_ASSERT_DELTA( new_x[6], 1.2158, 1e-3); 
        TS_ASSERT_DELTA( new_y[6], 0.5, 1e-3); 
    }


    
    void Test3dProblemOnCube() throw(Exception)
    {
        Vector<double> body_force(3);
        body_force(1) = 10.0;
    
        MooneyRivlinMaterialLaw<3> mooney_rivlin_law(2.0,2.0);

        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 0.1); 
        mesh.refine_global(2);
        FiniteElasticityTools<3>::SetFixedBoundary(mesh,0);

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



    void Test3dProblemOnCubeFixedDisplacement() throw(Exception)
    {
        Vector<double> body_force(3); // zero vector
        body_force(2)=5;
    
        MooneyRivlinMaterialLaw<3> mooney_rivlin_law(2.0,2.0);

        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(1);

        Triangulation<3>::cell_iterator element_iter = mesh.begin_active();
        while(element_iter!=mesh.end())
        {
            for(unsigned face_index=0; face_index<GeometryInfo<3>::faces_per_cell; face_index++)
            {
                if(element_iter->face(face_index)->at_boundary()) 
                {
                    double z = element_iter->face(face_index)->center()(2);
                    
                    if(fabs(z)<1e-4)
                    {
                        // z=0, label as fixed boundary
                        element_iter->face(face_index)->set_boundary_indicator(FIXED_BOUNDARY);
                    }
                    else if(fabs(z-1)<1e-4)
                    {
                        // z=1, label as dirichlet boundary
                        element_iter->face(face_index)->set_boundary_indicator(DIRICHLET_BOUNDARY);
                    }
                    else 
                    {
                        // z!=0 or 1, label as neumann boundary
                        element_iter->face(face_index)->set_boundary_indicator(NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }

        FiniteElasticityAssembler<3> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/fixeddisp3d");

        DoFHandler<3>& dof_handler = finite_elasticity.GetDofHandler();


        std::map<unsigned,double> boundary_values;
        
        std::vector<bool> component_mask(3+1); // dim+1
        component_mask[0] = true;
        component_mask[1] = true;
        component_mask[2] = true;
        component_mask[3] = false;

        VectorTools::interpolate_boundary_values(dof_handler,
                                                 FIXED_BOUNDARY,
                                                 ZeroFunction<3>(3+1),  // note the "+1" here! - number of components
                                                 boundary_values,
                                                 component_mask);

        assert(!boundary_values.empty());


        VectorTools::interpolate_boundary_values(dof_handler,
                                                 DIRICHLET_BOUNDARY,
                                                 ComponentSelectFunction<3>(2,-0.1,4),
                                                 boundary_values,
                                                 component_mask);



        finite_elasticity.SetBoundaryValues(boundary_values);

        
        finite_elasticity.Solve();


        Vector<double>& solution = finite_elasticity.GetSolutionVector();


        DofVertexIterator<3> dof_vertex_iter(&mesh, &dof_handler);

        
        while(!dof_vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = dof_vertex_iter.GetVertexGlobalIndex();
            Point<3> old_posn = dof_vertex_iter.GetVertex();
            
            Point<3> new_posn;
            new_posn(0) = old_posn(0)+solution(dof_vertex_iter.GetDof(0));
            new_posn(1) = old_posn(1)+solution(dof_vertex_iter.GetDof(1));
            new_posn(2) = old_posn(2)+solution(dof_vertex_iter.GetDof(2));
            
            // todo: TEST THESE!!

            std::cout << vertex_index << " " << old_posn(0) << " " << old_posn(1) << " " << old_posn(2) 
                                      << " " << new_posn(0) << " " << new_posn(1) << " " << new_posn(2)
                                      << "\n";

            //// UPDATE THE NODE POSITIONS
            // GetVertex returns a reference to a Point<DIM>, so this changes the mesh
            // directly. Do this so the new volume can be computed
            dof_vertex_iter.GetVertex()[0] = new_posn(0);         
            dof_vertex_iter.GetVertex()[1] = new_posn(1);         
            dof_vertex_iter.GetVertex()[2] = new_posn(2);         
                                      
            dof_vertex_iter.Next();
        }

        // compute the deformed volume
        // NOTE: this aren't very accurate volumes, since we have lost the
        // positions of the extra nodes (those used with quadratic basis functions)
        // and the measure() function below must use linear interpolation. Hence
        // the high tolerances
        double deformed_volume = 0.0;
        element_iter = mesh.begin_active();
        while(element_iter!=mesh.end())
        {
            double element_volume = element_iter->measure();
            TS_ASSERT_DELTA(element_volume, 1.0/mesh.n_active_cells(), 1e-2); 
            
            deformed_volume += element_volume;
            element_iter++;
        }
        
        TS_ASSERT_DELTA(deformed_volume, 1.0, 0.05);
    }
    
    
// fails - something odd happens with this choice of params - needs fixing..    
    void xTestOnHeterogeneousProblem()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        

        Triangulation<2>::cell_iterator element_iter = mesh.begin_active();
        while(element_iter!=mesh.end())
        {
            if(element_iter->center()[0] < 0.5) 
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


        Vector<double> body_force(2);
        body_force(1) = 3.0;
    
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law_stiff(10.0);
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law_weak(2.0);

        std::vector<AbstractIncompressibleMaterialLaw<2>*> material_laws;
        material_laws.push_back(&mooney_rivlin_law_stiff);
        material_laws.push_back(&mooney_rivlin_law_weak);


        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law_weak,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/heterogeneous2d");
 
        

        finite_elasticity.SetMaterialLawsForHeterogeneousProblem(material_laws, material_ids);

                                                         
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
                                      

            //// UPDATE THE NODE POSITIONS
            // GetVertex returns a reference to a Point<DIM>, so this changes the mesh
            // directly. Do this so the new volume can be computed
            vertex_iter.GetVertex()[0] = new_posn(0);         
            vertex_iter.GetVertex()[1] = new_posn(1);         
                                      
            vertex_iter.Next();
        }

//        // compute the deformed volume 
//        // NOTE: this aren't very accurate volumes, since we have lost the
//        // positions of the extra nodes (those used with quadratic basis functions)
//        // and the measure() function below must use linear interpolation. Hence
//        // the high tolerances
//        double deformed_volume = 0.0;
//        //Triangulation<2>::active_cell_iterator 
//        element_iter = mesh.begin_active();
//        while(element_iter!=mesh.end())
//        {
//            double element_volume = element_iter->measure();
//            TS_ASSERT_DELTA(element_volume, 1.0/mesh.n_active_cells(), 1e-2); 
//            
//            deformed_volume += element_volume;
//            element_iter++;
//        }
//        
//        TS_ASSERT_DELTA(deformed_volume, 1.0, 1e-2);
    }
};
#endif /*TESTFINITEELASTICITYASSEMBLER_HPP_*/
