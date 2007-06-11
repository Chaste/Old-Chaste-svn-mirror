#ifndef TESTFINITEELASTICITYASSEMBLERWITHGROWTH_HPP_
#define TESTFINITEELASTICITYASSEMBLERWITHGROWTH_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include "FiniteElasticityAssemblerWithGrowth.cpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

#include "MooneyRivlinMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "ExponentialMaterialLaw.hpp"

#include "FiniteElasticityTools.hpp"
#include "ConcentrationBasedTumourSourceModel.hpp"
#include "ConstantTumourSourceModel.hpp"

#include "grid/tria_boundary_lib.h"



// todos: proper test of answers, compare numerical jacobian
// sensible test once s set up. change constructor

class TumourGrowingDyingSourceModel : public AbstractGrowingTumourSourceModel<2>
{
private :
    Point<2> mCentre;
public :    
    TumourGrowingDyingSourceModel(double centre_x=0, double centre_y=0)
    {
        mCentre[0] = centre_x;
        mCentre[1] = centre_y;
    }
    
    void Run(double tStart, double tEnd, FiniteElasticityAssembler<2>* pFiniteElasticityAssembler)
    {
        std::map<unsigned,EvaluationPointInfo<2> >::iterator iter
            = this->mEvaluationPoints.begin();
        while (iter!=this->mEvaluationPoints.end())
        {
            //unsigned mesh_index = iter->first;
            Point<2>& position = iter->second.OldPosition;
            Point<2> diff = position-mCentre;
    
            double distance_to_centre = std::sqrt(diff.square());
            double source_value = 2*(distance_to_centre - 0.7);
            iter->second.SourceValue = source_value;
            iter++;
        }
    }
};


class TestFiniteElasticityAssemblerWithGrowth : public CxxTest::TestSuite
{
private :
    void MakeRectangularMeshWithTwoCircles(Triangulation<2>& mesh)
    {
        Point<2> zero;
        Point<2> opposite_corner;
        opposite_corner[0] = 1.3;
        opposite_corner[1] = 1;
        
        unsigned num_elem_x = 40;
        unsigned num_elem_y = 20;
        
        std::vector<unsigned> repetitions;
        repetitions.push_back(num_elem_x);
        repetitions.push_back(num_elem_y);
        
        GridGenerator::subdivided_hyper_rectangle(mesh, repetitions, zero, opposite_corner);
        
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        double radius = 0.21;
        
        Point<2> centre1;
        centre1[0]=0.5;
        centre1[1]=0.5;
        
        Point<2> centre2;
        centre2[0]=0.8;
        centre2[1]=0.5;
        
        // set all elements as non growing initially, then set circular region as growing
        FiniteElasticityTools<2>::SetAllElementsAsNonGrowingRegion(mesh);
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, centre1, radius);
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, centre2, radius);
    }
    
    
public :
    void TestExceptions() throw(Exception)
    {
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        ConstantTumourSourceModel<2> source_model(1.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        // elements haven't been set as GROWING_REGION or NON_GROWING_REGION
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssemblerWithGrowth<2> bad_fe_with_growth(&mesh,&mooney_rivlin_law,body_force,1.0,"",&source_model));
        
        FiniteElasticityTools<2>::SetAllElementsAsNonGrowingRegion(mesh);
        
        // no elements set as GROWING_REGION
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssemblerWithGrowth<2> bad_fe_with_growth2(&mesh,&mooney_rivlin_law,body_force,1.0,"",&source_model));
        
        // set the first element as growing
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        element_iter->set_material_id(GROWING_REGION);
        
        // should construct ok now
        FiniteElasticityAssemblerWithGrowth<2> fe_with_growth(&mesh,&mooney_rivlin_law,body_force,1.0,"",&source_model);
        
        // set times not been called
        TS_ASSERT_THROWS_ANYTHING(fe_with_growth.Run());
        
        // start time > end time
        TS_ASSERT_THROWS_ANYTHING(fe_with_growth.SetTimes(1.0, 0.0, 0.01));
        
        // dt negative
        TS_ASSERT_THROWS_ANYTHING(fe_with_growth.SetTimes(0.0, 1.0, -0.01));
        
        // none of the above should throw now
        TS_ASSERT_THROWS_NOTHING(fe_with_growth.SetTimes(0.0, 1.0, 0.01));
    }
    
    
    void NO__TestGrowingDyingTumourModel()
    {
        return;
        
        Vector<double> body_force(2); // zero
        double density = 1.0;
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(0.02);

        Triangulation<2> mesh;
        Point<2> zero;

//        GridGenerator::hyper_cube(mesh, 0.0, 2.0);
//        mesh.refine_global(3);
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<2> boundary(zero);
        mesh.set_boundary(0, boundary);
        mesh.refine_global(4);

        TriangulationVertexIterator<2> iter(&mesh);
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, iter.GetVertex());
        
        // set all elements as growing region (using a circle with a big radius)
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, zero, 100);
        
        // a source model which means death if the point is in the centre of the square
        // defined in this file
        TumourGrowingDyingSourceModel source_model;
        
        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/tumour_only",
                                                                      &source_model); 
                
        finiteelas_with_growth.SetTimes(0.0, 5, 0.05);
        finiteelas_with_growth.Run();
    }
    
    
    
    
    
    
    void NO__Test2dProblemOnSquare() throw(Exception)
    {
        Vector<double> body_force(2); // zero
        double density = 1.0;
        
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(4);
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        Point<2> centre;
        centre[0]=0.5;
        centre[1]=0.5;
        
        // set all elements as non growing initially, then set circular region as growing
        FiniteElasticityTools<2>::SetAllElementsAsNonGrowingRegion(mesh);
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, centre, 0.2);
        
        
        //ConstantTumourSourceModel<2> source_model(1.0);
        ConcentrationBasedTumourSourceModel<2> source_model(mesh);
        
        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/simple2d",
                                                                      &source_model);
                
                
        // loop over all the elements, and if it is in the growing region, check
        // each node has an ode system associated with it...
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            if (element_iter->material_id()==GROWING_REGION)
            {
                for (unsigned i=0; i<GeometryInfo<2>::vertices_per_cell; i++)
                {
                    unsigned vertex_index = element_iter->vertex_index(i);
                    TS_ASSERT_EQUALS(finiteelas_with_growth.IsGrowingNode(vertex_index), true);
                }
            }
            element_iter++;
        }
        
        
        finiteelas_with_growth.SetTimes(0.0, 10.0, 0.1);
        
        finiteelas_with_growth.Run();  
    }
    
    
    
    
    void NO_Test2dPolypFormation() throw(Exception)
    {
        Vector<double> body_force(2); // zero
        double density = 1.0;

        double length = 50;
        double height = 2;

        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(0.02);

        Triangulation<2> mesh;
        Point<2> zero;
        Point<2> opposite_corner;
        opposite_corner[0] = length;
        opposite_corner[1] = height;
        
        unsigned num_elem_x = 50;
        unsigned num_elem_y = 1;
        
        std::vector<unsigned> repetitions;
        repetitions.push_back(num_elem_x);
        repetitions.push_back(num_elem_y);
        
        GridGenerator::subdivided_hyper_rectangle(mesh, repetitions, zero, opposite_corner);


//
//        double alpha = 1;
//        double pi = 3.14159265;
//		TriangulationVertexIterator<2> vertex_iter(&mesh);
//		while(!vertex_iter.End())
//		{
//			Point<2>& position = vertex_iter.GetVertex();
//			position[1] += alpha*sin(3*pi*position[0]/length); 
//			vertex_iter.Next();
//		}
    


        // set all elements as growing region (using a circle with a big radius)
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, zero, 10*length);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, length, false);

        ConcentrationBasedTumourSourceModel<2> source_model(mesh);
        
        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/polyp2d",
                                                                      &source_model);
                
        finiteelas_with_growth.SetTimes(0.0, 10, 0.1);
        finiteelas_with_growth.Run();
    }
    
    
    
    
    
    void NO__Test2dProblemOnRectangleTwoCircles() throw(Exception)
    {
        Vector<double> body_force(2); // zero
        double density = 1.0;

        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(0.002);
        
        Triangulation<2> mesh;
        MakeRectangularMeshWithTwoCircles(mesh);
        
        //ConstantTumourSourceModel<2> source_model(1.0);
        ConcentrationBasedTumourSourceModel<2> source_model(mesh);
        
        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth(&mesh,
                                                                      NULL,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/rectanle2d_2circles",
                                                                      &source_model);
                
                
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law_stiff(0.02);
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law_weak(0.002);
        
        std::vector<unsigned> material_ids;
        std::vector<AbstractIncompressibleMaterialLaw<2>*> material_laws;
        
        material_ids.push_back(GROWING_REGION);
        material_laws.push_back(&mooney_rivlin_law_weak);
        
        material_ids.push_back(NON_GROWING_REGION);
        material_laws.push_back(&mooney_rivlin_law_stiff);
        
        finiteelas_with_growth.SetMaterialLawsForHeterogeneousProblem(material_laws,
                                                                      material_ids);
                
        finiteelas_with_growth.SetTimes(0.0, 10.0, 0.1);
        finiteelas_with_growth.Run();
    }
};
#endif /*TESTFINITEELASTICITYASSEMBLERWITHGROWTH_HPP_*/
