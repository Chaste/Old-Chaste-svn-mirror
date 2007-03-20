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
    void TestUpdateMesh()
    {
        Vector<double> body_force(2);
        body_force(0) = 1.0;
        
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        
        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");
                                                       
        ConcentrationBasedTumourSourceModel<2> source_model(mesh);
        
        // add the first few nodes
        TriangulationVertexIterator<2> iter(&mesh);
        
        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";
        source_model.AddEvaluationPoint(0, iter.GetVertex(), iter.GetVertexGlobalIndex());
        iter.Next();
        
        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";
        source_model.AddEvaluationPoint(1, iter.GetVertex(), iter.GetVertexGlobalIndex());
        
        source_model.Run(0,1,&finite_elasticity);
        
        TriangulationVertexIterator<2> fe_mesh_iter(finite_elasticity.GetMesh());
        TriangulationVertexIterator<2> source_mesh_iter(&source_model.mMesh);
        
        TS_ASSERT_EQUALS(finite_elasticity.GetMesh()->n_vertices(), source_model.mMesh.n_vertices());
        while (!fe_mesh_iter.ReachedEnd())
        {
            Point<2> fe_point = fe_mesh_iter.GetVertex();
            Point<2> source_point = source_mesh_iter.GetVertex();
            Point<2> diff = fe_point-source_point;
            
            TS_ASSERT_DELTA(std::sqrt(diff.square()), 0.0, 1e-10);
            
            fe_mesh_iter.Next();
            source_mesh_iter.Next();
        }
        
        finite_elasticity.Solve();
        
        source_model.Run(0,1,&finite_elasticity);
// check visually and write test...
    }
    
};
#endif /*TESTCONCENTRATIONBASEDTUMOURSOURCEMODEL_HPP_*/
