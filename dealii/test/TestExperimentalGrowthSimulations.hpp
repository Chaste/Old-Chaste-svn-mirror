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


#ifndef TESTEXPERIMENTGROWTHSIMULATIONS_HPP_
#define TESTEXPERIMENTGROWTHSIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include "FiniteElasticityAssemblerWithGrowth.cpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

#include "MooneyRivlinMaterialLaw2.hpp"
#include "PolynomialMaterialLaw3d2.hpp"
#include "ExponentialMaterialLaw2.hpp"

#include "FiniteElasticityTools.hpp"
#include "ConcentrationBasedTumourSourceModel.hpp"
#include "ConstantTumourSourceModel.hpp"
#include "grid/tria_boundary_lib.h"
#include "ConstantTumourSourceModel.hpp"


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

class TestExperimentalGrowthSimulations : public CxxTest::TestSuite
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
    // silly hack to avoid any 'no tests defined' errors if all tests are 'NO_Test'ed out
    // (this test suite is in the ExtraSimulations test pack so won't be run)
    void TestOnePlusOneEqualsTwo()
    {
        TS_ASSERT_EQUALS(1+1,2);
    }



    void NO__TestGrowingDyingTumourModel()
    {
        Vector<double> body_force(2); // zero
        double density = 1.0;
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(0.02);

        Triangulation<2> mesh;
        Point<2> zero;

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

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(2.0);

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
                                                                      "finite_elas_growth/square2d",
                                                                      &source_model);


        finiteelas_with_growth.SetTimes(0.0, 10.0, 0.1);

        finiteelas_with_growth.Run();
    }


    void Test2dPolypFormationWithElasticUnderside() throw(Exception)
    {
        Vector<double> body_force(2); // zero
        double density = 1.0;

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(0.02);

        Triangulation<2> mesh;
        Point<2> zero;
        Point<2> opposite_corner;
        opposite_corner[0] = 1.0;
        opposite_corner[1] = 1.0;

        unsigned num_elem_x = 20;
        unsigned num_elem_y = 20;

        std::vector<unsigned> repetitions;
        repetitions.push_back(num_elem_x);
        repetitions.push_back(num_elem_y);

        GridGenerator::subdivided_hyper_rectangle(mesh, repetitions, zero, opposite_corner);

        // set all elements as growing region (using a circle with a big radius)
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 1.0, false);

        FiniteElasticityTools<2>::SetAllElementsAsNonGrowingRegion(mesh);

        bool found_an_element = false;

        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            for (unsigned i=0; i<GeometryInfo<2>::vertices_per_cell; i++)
            {
                double x = element_iter->vertex(i)[0];
                double y = element_iter->vertex(i)[1];
                if ( (fabs(y-1.0) < 1e-4) && (x>0.4) && (x<0.6) )
                {
                    found_an_element = true;
                    element_iter->set_material_id(GROWING_REGION);
                }
            }
            element_iter++;
        }

        assert(found_an_element);

        double source_value = 1;
        ConstantTumourSourceModel<2> source_model(source_value);

        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/polyp2d_elastic_under",
                                                                      &source_model);

        finiteelas_with_growth.SetTimes(0.0, 3, 0.1);
        finiteelas_with_growth.Run();
    }


    void NO_Test2dPolypFormation() throw(Exception)
    {
        Vector<double> body_force(2); // zero
        double density = 1.0;

        double length = 50;
        double height = 2;

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(0.02);

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
//      TriangulationVertexIterator<2> vertex_iter(&mesh);
//      while(!vertex_iter.ReachedEnd())
//      {
//          Point<2>& position = vertex_iter.GetVertex();
//          position[1] += alpha*sin(3*pi*position[0]/length);
//          vertex_iter.Next();
//      }



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

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(0.002);

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


        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law_stiff(0.02);
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law_weak(0.002);

        std::vector<unsigned> material_ids;
        std::vector<AbstractIncompressibleMaterialLaw2<2>*> material_laws;

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

#endif /*TESTEXPERIMENTGROWTHSIMULATIONS_HPP_*/
