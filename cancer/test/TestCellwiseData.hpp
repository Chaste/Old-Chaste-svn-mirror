/*

Copyright (C) University of Oxford, 2008

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
#ifndef TESTCELLWISEDATA_HPP_
#define TESTCELLWISEDATA_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "CellwiseData.hpp"
#include "FixedCellCycleModelCellsGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "MeshArchiveInfo.hpp"

/**
 * This class contains tests for methods on the class CellwiseData.
 */
class TestCellwiseData : public AbstractCancerTestSuite
{
public:
	

    void TestCellwiseDataSimple() throw(Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh);

        // Create a tissue
        MeshBasedTissue<2> tissue(mesh,cells);

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());

        // One variable tests
        
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());
        TS_ASSERT_THROWS_ANYTHING(p_data->SetTissue(tissue));

        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);

        TS_ASSERT(!CellwiseData<2>::Instance()->IsSetUp());

        p_data->SetTissue(tissue);
        
        TS_ASSERT(CellwiseData<2>::Instance()->IsSetUp());

        p_data->SetValue(1.23, mesh.GetNode(0));
        AbstractTissue<2>::Iterator iter = tissue.Begin();
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter)), 1.23, 1e-12);

        p_data->SetValue(2.23, mesh.GetNode(1));
        ++iter;
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter)), 2.23, 1e-12);

        // Test ReallocateMemory method
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

        // Two variable tests

        p_data = CellwiseData<2>::Instance();

        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 2);
        p_data->SetTissue(tissue);
        
        TS_ASSERT_THROWS_ANYTHING(p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1));
        TS_ASSERT(CellwiseData<2>::Instance()->IsSetUp());

        p_data->SetValue(3.23, mesh.GetNode(0), 1);
        AbstractTissue<2>::Iterator iter2 = tissue.Begin();
        
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 1), 3.23, 1e-12);

        p_data->SetValue(4.23, mesh.GetNode(1), 1);
        ++iter2;
        
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 1), 4.23, 1e-12);

        // Other values should have been initialised to zero
        ++iter2;
        TS_ASSERT_DELTA( p_data->GetValue(&(*iter2), 0), 0.0, 1e-12);

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestArchiveCellwiseData()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Get each a birth time of -node_index, so the age = node_index
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh);

        // Create a tissue
        MeshBasedTissue<2> tissue(mesh,cells);

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

            TS_ASSERT(p_data->IsSetUp());

            // Create an ouput archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            // Write to the archive
            output_arch << static_cast<const CellwiseData<2>&>(*CellwiseData<2>::Instance());

            CellwiseData<2>::Destroy();
        }

        {
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            MeshArchiveInfo::meshPathname = "mesh/test/data/square_4_elements";
            input_arch >> *p_data;

            // Check the data
            TS_ASSERT(CellwiseData<2>::Instance()->IsSetUp());
            TS_ASSERT(p_data->IsSetUp());
            TS_ASSERT(!p_data->mUseConstantDataForTesting);

            for (AbstractTissue<2>::Iterator iter = tissue.Begin();
                 iter != tissue.End();
                 ++iter)
            {
                TS_ASSERT_DELTA(p_data->GetValue(&(*iter), 0), (double) tissue.GetNodeCorrespondingToCell(*iter)->GetIndex(), 1e-12);
            }

            // Tidy up
            delete p_data->mpTissue;
            CellwiseData<2>::Destroy();
        }
    }

};

#endif /*TESTCELLWISEDATA_HPP_*/
