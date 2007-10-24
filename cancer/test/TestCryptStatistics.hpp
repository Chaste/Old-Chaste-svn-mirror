#ifndef TESTCRYPTSTATISTICS_HPP_
#define TESTCRYPTSTATISTICS_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "Tissue.cpp"
#include "CryptStatistics.hpp"
#include "CellsGenerator.hpp"



class TestCryptStatistics : public CxxTest::TestSuite
{
public:
    void TestGetSection() throw (Exception)
    {        
        CancerParameters::Instance()->Reset();

        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateForCrypt(cells, *p_mesh, FIXED, true);// true = mature cells

        Tissue<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        CryptStatistics crypt_statistics(crypt);
        
        std::vector< TissueCell* > test_section=crypt_statistics.GetCryptSection(0.5,1.5,sqrt(3));
        
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section.size(), 6u);

        unsigned expected_indices[6]={0,1,3,4,7,8};

        for(unsigned i=0; i<test_section.size(); i++)
        {
            TissueCell* cell = test_section[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices[i]);
        }

        
        // Test that we get a valid section when the x-values are the same
        std::vector< TissueCell* > test_section_vertical=crypt_statistics.GetCryptSection(0.5,0.5,sqrt(3));
        
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section_vertical.size(), 5u);

        unsigned expected_indices_vertical[6]={0,1,3,6,7};

        for(unsigned i=0; i<test_section_vertical.size(); i++)
        {
            TissueCell* cell = test_section_vertical[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices_vertical[i]);
        }
        
        std::vector< TissueCell* > test_section_periodic=crypt_statistics.GetCryptSectionPeriodic(0.5,2.5,sqrt(3));
        
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic.size(), 6u);
        
        unsigned expected_indices_periodic[6]={0,1,3,5,6,8};
        
        for(unsigned i=0; i<test_section_periodic.size(); i++)
        {
            TissueCell* cell = test_section_periodic[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices_periodic[i]);
        }
        
             
        std::vector< TissueCell* > test_section_periodic_2=crypt_statistics.GetCryptSectionPeriodic(2.5,0.5,sqrt(3));
        
        //Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic_2.size(), 6u);
        
        unsigned expected_indices_periodic_2[6]={0,2,3,5,6,7};
        
        for(unsigned i=0; i<test_section_periodic_2.size(); i++)
        {
            TissueCell* cell = test_section_periodic_2[i];
            TS_ASSERT_EQUALS( crypt.GetNodeCorrespondingToCell(*cell)->GetIndex(), 
                              expected_indices_periodic_2[i]);
        }
        
        SimulationTime::Destroy();
    }
   
};

#endif /*TESTCRYPTSTATISTICS_HPP_*/
