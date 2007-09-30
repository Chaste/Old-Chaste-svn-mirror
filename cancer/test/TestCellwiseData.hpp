#ifndef TESTCELLWISEDATA_HPP_
#define TESTCELLWISEDATA_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <vector>
#include "Crypt.cpp"
#include "CellwiseData.cpp"
#include "CellsGenerator.hpp"


class TestCellwiseData : public CxxTest::TestSuite
{
public:
    void TestCellwiseDataSimple()
    {
        // set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<MeinekeCryptCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, mesh);

        // create a crypt
        Crypt<2> crypt(mesh,cells);

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());
        
        // 1 variable tests
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        
        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());
        
        TS_ASSERT_THROWS_ANYTHING(p_data->SetCrypt(crypt)); 
        
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());

        p_data->SetCrypt(crypt);     

        TS_ASSERT(CellwiseData<2>::Instance()->IsSetUp());
        
        p_data->SetValue(1.23, mesh.GetNode(0));
        Crypt<2>::Iterator iter = crypt.Begin();
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter)), 1.23, 1e-12);

        p_data->SetValue(2.23, mesh.GetNode(1));
        ++iter;
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter)), 2.23, 1e-12);
        
        p_data->Destroy();

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());

        // 2 variable test

        p_data = CellwiseData<2>::Instance();
        
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 2);
        p_data->SetCrypt(crypt);     
        TS_ASSERT_THROWS_ANYTHING(p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1));

        TS_ASSERT(CellwiseData<2>::Instance()->IsSetUp());
        
        p_data->SetValue(3.23, mesh.GetNode(0), 1);
        Crypt<2>::Iterator iter2 = crypt.Begin();
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 1), 3.23, 1e-12);

        p_data->SetValue(4.23, mesh.GetNode(1), 1);
        ++iter2;
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 1), 4.23, 1e-12);

        //  other values should have been initialised to zero        
        ++iter2;
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 0), 0.0, 1e-12);
    }
};        

#endif /*TESTCELLWISEDATA_HPP_*/
