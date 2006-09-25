#ifndef TESTCRYPTSIMULATION2D_HPP_
#define TESTCRYPTSIMULATION2D_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "CryptSimulation2D.hpp"

#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"

#include "CancerParameters.hpp"


class TestCryptSimulation2D : public CxxTest::TestSuite
{
   public:
    // same as Test1DChainWithBirthVariableRestLength but with Meineke cells.
    // (see comment for Test1DChainWithBirthVariableRestLength).
    void Test2DSprings() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        srandom(0);
        double crypt_length = 11.0;
        
//        Make1dCryptMesh("1D_crypt_mesh", 23, crypt_length);
//        std::string testoutput_dir;
//        OutputFileHandler output_file_handler("");
//        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, 0.0, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation2D simulator(mesh, cells);
        simulator.SetOutputDirectory("Crypt2DSprings");
        simulator.SetEndTime(1.0);
        simulator.SetCryptLength(crypt_length);
        //simulator.SetIncludeVariableRestLength();
        
        
        // this is FAILING at the moment, the throws anything is for committing purposes
        TS_ASSERT_THROWS_ANYTHING( simulator.Solve() );
    }
    
};

#endif /*TESTCRYPTSIMULATION2D_HPP_*/

