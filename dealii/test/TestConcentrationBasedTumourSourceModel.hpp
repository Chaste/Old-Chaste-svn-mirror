/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTCONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_
#define TESTCONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_
#include <cxxtest/TestSuite.h>

#include "AbstractGrowingTumourSourceModel.hpp"
#include "SimpleTumourSourceModel.hpp"
#include "ConstantTumourSourceModel.hpp"
#include "ConcentrationBasedTumourSourceModel.hpp"


#include "FiniteElasticityAssembler.hpp"
#include "FiniteElasticityTools.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "TriangulationVertexIterator.hpp"

class TestConcentrationBasedTumourSourceModel : public CxxTest::TestSuite
{
public:

//    void testConcentrationBasedTumourSourceModel()
//    {
//        Vector<double> body_force(2);
//        body_force(0) = 1.0;
//    
//        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
//
//        Triangulation<2> mesh;
//        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
//        mesh.refine_global(1);
//        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
//        
//        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
//                                                       &mooney_rivlin_law,
//                                                       body_force,
//                                                       1.0,
//                                                       "");
//
//        ConcentrationBasedTumourSourceModel<2> source_model(mesh);
//        
//        // add the first few nodes
//        TriangulationVertexIterator<2> iter(&mesh); 
//
//        unsigned indices[2];
//        
//        indices[0] = iter.GetVertexGlobalIndex();
//        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";  
//        source_model.AddEvaluationPoint(indices[0], iter.GetVertex());
//        iter.Next();
//
//        indices[1] = iter.GetVertexGlobalIndex();
//        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";  
//        source_model.AddEvaluationPoint(indices[1], iter.GetVertex());
//
//        source_model.Run(0,1,&finite_elasticity);
//    }


    void TestUpdateMesh()
    {
        Vector<double> body_force(2);
        body_force(0) = 1.0;
        
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");
                                                       
        ConcentrationBasedTumourSourceModel<2> source_model(mesh);
        
        TriangulationVertexIterator<2> fe_mesh_iter(finite_elasticity.GetMesh());
        TriangulationVertexIterator<2> source_mesh_iter(&source_model.mDeformedMesh);
        
        TS_ASSERT_EQUALS(finite_elasticity.GetMesh()->n_vertices(), source_model.mDeformedMesh.n_vertices());
        while (!fe_mesh_iter.ReachedEnd())
        {
            TS_ASSERT_EQUALS(fe_mesh_iter.GetVertexGlobalIndex(), source_mesh_iter.GetVertexGlobalIndex());
            
            Point<2> fe_point = fe_mesh_iter.GetVertex();
            Point<2> source_point = source_mesh_iter.GetVertex();
            Point<2> diff = fe_point-source_point;
            
            TS_ASSERT_DELTA(std::sqrt(diff.square()), 0.0, 1e-10);
            
            fe_mesh_iter.Next();
            source_mesh_iter.Next();
        }
        
        // solve the static fe problem...
        finite_elasticity.StaticSolve();
        
        // update again..
        source_model.Run(0,1,&finite_elasticity);
        
        // check (fe_mesh+deformed_posn) = source_mesh
        fe_mesh_iter.Reset();
        source_mesh_iter.Reset();
        TS_ASSERT_EQUALS(finite_elasticity.GetMesh()->n_vertices(), source_model.mDeformedMesh.n_vertices());
        std::vector<Vector<double> > deformed_position = finite_elasticity.rGetDeformedPosition();
        while (!fe_mesh_iter.ReachedEnd())
        {
            TS_ASSERT_EQUALS(fe_mesh_iter.GetVertexGlobalIndex(), source_mesh_iter.GetVertexGlobalIndex());
            
            Point<2> fe_deformed_posn;
            fe_deformed_posn[0] = deformed_position[0](fe_mesh_iter.GetVertexGlobalIndex());
            fe_deformed_posn[1] = deformed_position[1](fe_mesh_iter.GetVertexGlobalIndex());

            Point<2> source_point = source_mesh_iter.GetVertex();
            Point<2> diff = fe_deformed_posn-source_point;
            
            TS_ASSERT_DELTA(std::sqrt(diff.square()), 0.0, 1e-10);
            
            fe_mesh_iter.Next();
            source_mesh_iter.Next();
        }
    }
};
#endif /*TESTCONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_*/
