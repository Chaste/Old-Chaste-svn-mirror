#ifndef TESTCELLSGENERATOR_HPP_
#define TESTCELLSGENERATOR_HPP_
#include <cxxtest/TestSuite.h>
#include "CellsGenerator.hpp"
#include "TissueCell.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "SimulationTime.hpp"

/**
 * Note that all these tests call setUp() and tearDown() before running,
 * so if you copy them into a new test suite be sure to copy these methods
 * too.
 */
class TestCellsGenerator : public CxxTest::TestSuite
{   
private:

    void setUp()
    {
        // Initialise singleton classes
        SimulationTime::Instance()->SetStartTime(0.0);
    }
    void tearDown()
    {
        // Clear up singleton classes
        SimulationTime::Destroy();
    }

public:

    void TestCellsGeneratorBasic() throw(Exception)
    {        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        std::vector<TissueCell> cells;        
        
        CellsGenerator<2> generator;
        generator.GenerateBasic(cells, mesh);
        
        TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());
        
        for (unsigned i=0 ; i<cells.size() ; i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetNodeIndex(), i);
            TS_ASSERT_DELTA(cells[i].GetBirthTime(), -(double)(i), 1e-9);   
        }
    }
    
    
    void TestSimpleCellsGeneratorForCryptRandom() throw(Exception)
    {        
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;
        
        CellsGenerator<2> generator;
        
        std::vector<TissueCell> cells;  
              
        double y0 = 0.2;
        double y1 = 1.0;
        double y2 = 2.0;
        double y3 = 3.0;
        
        generator.GenerateForCrypt(cells, *p_mesh, FIXED, true, y0, y1, y2 ,y3 );
        
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());
        
        for (unsigned i=0 ; i<cells.size() ; i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetNodeIndex(), i);
            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = cells[i].GetCellCycleModel()->GetGeneration();
            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
            }
        }
    }
    
    
    void TestSimpleCellsGeneratorForCryptNonRandom() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;
        
        CellsGenerator<2> generator;
        
        std::vector<TissueCell> cells;        
        generator.GenerateForCrypt(cells, *p_mesh, STOCHASTIC, false);
        
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());
        double y0 = 0.3;
        double y1 = 2.0;
        double y2 = 3.0;
        double y3 = 4.0;
        for (unsigned i=0 ; i<cells.size() ; i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetNodeIndex(), i);
            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = cells[i].GetCellCycleModel()->GetGeneration();
            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
            }
            
            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);   
        }
    }
    
    
    void TestOdeCellsGeneratorForCryptRandom() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;
        
        CellsGenerator<2> generator;
        
        std::vector<TissueCell> cells;        
        generator.GenerateForCrypt(cells, *p_mesh, TYSONNOVAK, true);
        
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());
        
        for (unsigned i=0 ; i<cells.size() ; i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetNodeIndex(), i);
        }
    }
    
    
    void TestOdeCellsGeneratorForCryptNonRandom() throw(Exception)
    {
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();;
        
        CellsGenerator<2> generator;
        
        std::vector<TissueCell> cells;        
        generator.GenerateForCrypt(cells, *p_mesh, WNT, false);
        
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());
        
        for (unsigned i=0 ; i<cells.size() ; i++)
        {
            TS_ASSERT_EQUALS(cells[i].GetNodeIndex(), i);
            
            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);   
        }
        
        // Coverage of different cell types, all work like WNT.
        generator.GenerateForCrypt(cells, *p_mesh, SIMPLE_WNT, false);
        generator.GenerateForCrypt(cells, *p_mesh, INGE_WNT_SWAT_HYPOTHESIS_ONE, false);
        generator.GenerateForCrypt(cells, *p_mesh, INGE_WNT_SWAT_HYPOTHESIS_TWO, false);
        generator.GenerateForCrypt(cells, *p_mesh, STOCHASTIC_WNT, false);
    }
    
};

#endif /*TESTCELLSGENERATOR_HPP_*/
