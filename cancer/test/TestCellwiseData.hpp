#ifndef TESTCELLWISEDATA_HPP_
#define TESTCELLWISEDATA_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include <cmath>
#include <vector>
#include "Tissue.cpp"
#include "CellwiseData.cpp"
#include "CellsGenerator.hpp"


class TestCellwiseData : public CxxTest::TestSuite
{
public:

    void TestCellwiseDataSimple() throw(Exception)
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
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, mesh);

        // create a tissue
        Tissue<2> tissue(mesh,cells);

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());
        
        // 1 variable tests
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        
        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());
                
        TS_ASSERT_THROWS_ANYTHING(p_data->SetTissue(tissue)); 
        
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());

        p_data->SetTissue(tissue);     

        TS_ASSERT(CellwiseData<2>::Instance()->IsSetUp());
        
        p_data->SetValue(1.23, mesh.GetNode(0));
        Tissue<2>::Iterator iter = tissue.Begin();
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter)), 1.23, 1e-12);

        p_data->SetValue(2.23, mesh.GetNode(1));
        ++iter;
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter)), 2.23, 1e-12);
        
        // test ReallocateMemory method
        TissueCell new_cell(STEM, HEALTHY, new FixedCellCycleModel());
        new_cell.SetBirthTime(-1);
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 0.2;
        new_cell_location[1] = 0.3;
        tissue.AddCell(new_cell,new_cell_location); 
                
        TS_ASSERT_THROWS_NOTHING(p_data->ReallocateMemory());
        TS_ASSERT_EQUALS(p_data->mData.size(), tissue.rGetMesh().GetNumNodes());
                
        p_data->Destroy();

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());
        
        // 2 variable test

        p_data = CellwiseData<2>::Instance();
        
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 2);
        p_data->SetTissue(tissue);     
        TS_ASSERT_THROWS_ANYTHING(p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1));

        TS_ASSERT(CellwiseData<2>::Instance()->IsSetUp());
        
        p_data->SetValue(3.23, mesh.GetNode(0), 1);
        Tissue<2>::Iterator iter2 = tissue.Begin();
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 1), 3.23, 1e-12);

        p_data->SetValue(4.23, mesh.GetNode(1), 1);
        ++iter2;
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 1), 4.23, 1e-12);

        //  other values should have been initialised to zero        
        ++iter2;
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 0), 0.0, 1e-12);
        
        SimulationTime::Destroy();
        CellwiseData<2>::Destroy();
    }
    
    
    void TestArchiveCellwiseData()
    {
        // Set up the simulation time object so the cells can be created
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells, one for each node. Get each a birth time of -node_index, so the age = node_index
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, mesh);

        // Create a tissue
        Tissue<2> tissue(mesh,cells);
        
        // Work out where to put the archive
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "cellwise_data.arch";

        {
            // Set up the data store
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();
            p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);
            p_data->SetTissue(tissue);  
            
            // Put some data in
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                p_data->SetValue((double) i, mesh.GetNode(i), 0);
            }
            
            // Write to the archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const CellwiseData<2>&>(*CellwiseData<2>::Instance());
            
            CellwiseData<2>::Destroy();
        }
        
        {
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // Restore from the archive
            Tissue<2>::meshPathname = "mesh/test/data/square_4_elements";
            input_arch >> *p_data;
            
            // Check the data
            TS_ASSERT(CellwiseData<2>::Instance()->IsSetUp());
            TS_ASSERT(p_data->IsSetUp());
            TS_ASSERT(!p_data->mUseConstantDataForTesting);
            
            for (Tissue<2>::Iterator iter = tissue.Begin();
                 iter != tissue.End();
                 ++iter)
            {
                TS_ASSERT_DELTA(p_data->GetValue(&(*iter), 0),(double) iter->GetNodeIndex(), 1e-12);
            }
            
            CellwiseData<2>::Destroy();
        }
        
        SimulationTime::Destroy();
    }
    
};        

#endif /*TESTCELLWISEDATA_HPP_*/
