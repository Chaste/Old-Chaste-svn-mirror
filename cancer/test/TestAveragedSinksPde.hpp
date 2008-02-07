#ifndef TESTAVERAGEDSINKSPDE_HPP_
#define TESTAVERAGEDSINKSPDE_HPP_

#include <cxxtest/TestSuite.h>

#include "AveragedSinksPde.hpp"
#include "CellsGenerator.hpp"

 
class TestAveragedSinksPde : public CxxTest::TestSuite
{
public:
    void TestAveragedSinksPde1()
    {
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);
        CancerParameters::Instance()->Reset();
        
        // Set up tissue        
        HoneycombMeshGenerator generator(5, 5, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, *p_mesh);
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        
        // create a coarse mesh - element 1 contains all the cells, 
        // element 0 contains none
        ConformingTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(1,1);
        coarse_mesh.Scale(100,100);

        // Set up PDE
        AveragedSinksPde<2> pde(tissue, -1.0);
        pde.SetupSourceTerms(coarse_mesh);

        // test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(0)); 
        double value_at_elem_1 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(1)); 

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(value_at_elem_1, -(tissue.GetNumRealCells()/coarse_mesh.GetElement(1)->GetVolume()), 1e-6);

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTAVERAGEDSINKSPDE_HPP_*/
