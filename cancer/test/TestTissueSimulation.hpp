#ifndef TESTTISSUESIMULATION_HPP_
#define TESTTISSUESIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <stdio.h>
#include <time.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <vector>
#include <iostream>

#include "TissueSimulation.cpp"
#include "FixedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ColumnDataReader.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/*
 *  
 *  
 *  Note: Most tests of TissueSimulation are in TestCryptSimulation2d
 * 
 * 
 */

class TestTissueSimulation : public AbstractCancerTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCancerTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCancerTestSuite::tearDown();
    }
    
public:

    void TestOutputStatistics() throw(Exception)
    {
        EXIT_IF_PARALLEL; //defined in PetscTools
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        CancerParameters::Instance()->SetHepaOneParameters();
        
        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
                    
        // Set up cells
        std::vector<TissueCell> cells;        
        
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new TysonNovakCellCycleModel());
            double birth_time = -1.0*p_gen->ranf();
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        // Set up tissue        
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        tissue.SetWriteTissueAreas(true); // record the spheroid radius and necrotic radius
              
        Meineke2001SpringSystem<2>* p_spring_system = new Meineke2001SpringSystem<2>(tissue);
        p_spring_system->UseCutoffPoint(1.5);
                  
        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, p_spring_system);
        simulator.SetOutputDirectory("TissueSimulationWritingProteins");
        simulator.SetEndTime(0.5);
        simulator.SetOutputCellVariables(true);
        simulator.SetOutputCellCyclePhases(true);

        // Run tissue simulation 
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        OutputFileHandler handler("TissueSimulationWritingProteins",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellvariables.dat";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/TissueSimulationWritingProteins/cellvariables.dat").c_str()), 0);
    }
};
#endif /*TESTTISSUESIMULATION_HPP_*/
