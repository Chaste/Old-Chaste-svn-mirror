#ifndef TESTNUTRIENTPDES_HPP_
#define TESTNUTRIENTPDES_HPP_

#include <cxxtest/TestSuite.h>

#include "MeshBasedTissue.hpp"
#include "SimpleNutrientPde.hpp"
#include "CellwiseNutrientSinkPde.hpp"
#include "AveragedSinksPde.hpp"
#include "CellsGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"
 
class TestNutrientPdes : public AbstractCancerTestSuite
{
public:

    void TestSimpleNutrientPde()
    {
        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
                
        // Set up PDE
        SimpleNutrientPde<2> pde(1.0);
        
        // Test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = pde.ComputeConstantInUSourceTerm(unused_point); 
        double value_at_elem_1 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,p_mesh->GetElement(0)); 

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(value_at_elem_1, -1.0, 1e-6);
    }
    
    void TestCellwiseNutrientSinkPde()
    {        
        // Set up tissue        
        HoneycombMeshGenerator generator(5, 5, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, *p_mesh);
        
        // Make one cell necrotic
        cells[0].SetCellType(NECROTIC);
        
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        
        // Set up PDE
        CellwiseNutrientSinkPde<2> pde(tissue, 1.0);
        
        // Test compute source terms
        ChastePoint<2> unused_point;
        double constant_in_u_source_term = pde.ComputeConstantInUSourceTerm(unused_point);
        
        TS_ASSERT_DELTA(constant_in_u_source_term, 0.0, 1e-6);
        
        Node<2>* p_node_0 = tissue.GetNodeCorrespondingToCell(tissue.rGetCellAtNodeIndex(0));
        Node<2>* p_node_1 = tissue.GetNodeCorrespondingToCell(tissue.rGetCellAtNodeIndex(1));
        
        double source_term_at_node_0 = pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node_0);
        double source_term_at_node_1 = pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_node_1);
                
        TS_ASSERT_DELTA(source_term_at_node_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(source_term_at_node_1, -1.0, 1e-6);

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void TestAveragedSinksPde()
    {        
        // Set up tissue        
        HoneycombMeshGenerator generator(5, 5, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, *p_mesh);
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        
        // Create a coarse mesh - element 1 contains all the cells, 
        // element 0 contains none
        ConformingTetrahedralMesh<2,2> coarse_mesh;
        coarse_mesh.ConstructRectangularMesh(1,1);
        coarse_mesh.Scale(100,100);

        // Set up PDE
        AveragedSinksPde<2> pde(tissue, -1.0);
        pde.SetupSourceTerms(coarse_mesh);

        // Test Compute source term
        ChastePoint<2> unused_point;
        double value_at_elem_0 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(0)); 
        double value_at_elem_1 = pde.ComputeLinearInUCoeffInSourceTerm(unused_point,coarse_mesh.GetElement(1)); 

        TS_ASSERT_DELTA(value_at_elem_0, 0.0, 1e-6);
        TS_ASSERT_DELTA(value_at_elem_1, -(tissue.GetNumRealCells()/coarse_mesh.GetElement(1)->GetVolume()), 1e-6);

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTNUTRIENTPDES_HPP_*/
