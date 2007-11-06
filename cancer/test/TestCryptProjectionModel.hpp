#ifndef TESTCRYPTPROJECTIONMODEL_HPP_
#define TESTCRYPTPROJECTIONMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "ColumnDataReader.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "CryptProjectionCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"

// This is Alex F experimental work and will be moved to my project folder 
// as soon as it can be!
    
class TestCryptProjectionModel : public CxxTest::TestSuite
{
public:

	/*
	 *  Test a 2D tissue simulation with a radial sloughing cell killer 
	 *  and the crypt projection cell cycle model, which depends on a radial 
	 *  Wnt gradient
	 */
	void TestMonolayerWithKillerAndCellCycleModel() throw (Exception)
	{        
        // instantiate singleton objects
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // set up mesh
        int num_cells_depth = 10;
        int num_cells_width = 10;     
        unsigned thickness_of_ghost_layer = 2;   
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);
        
        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(width_extremes[1] - width_extremes[0]);
		double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(height_extremes[1] - height_extremes[0]);
		
        p_mesh->Translate(-width_of_mesh/2,-height_of_mesh/2);                   
        
        // set up cells
        std::vector<TissueCell> cells;
        
        unsigned counter = 0;
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
        	if ( norm_2(p_mesh->GetNode(i)->rGetLocation()) <= sqrt(2) )
        	{
        		if (counter==1)
        		{
        			TissueCell cell(STEM, LABELLED, 0, new CryptProjectionCellCycleModel());
	            	double birth_time = -p_gen->ranf()*(p_params->GetStemCellG1Duration()
	                                             +p_params->GetSG2MDuration());
	                cell.SetNodeIndex(i);
	            	cell.SetBirthTime(birth_time);
	            	cell.SetSymmetricDivision();
	            	cells.push_back(cell);
	            	counter++;
        		}
        		else
        		{
        			TissueCell cell(STEM, HEALTHY, 0, new CryptProjectionCellCycleModel());
	            	double birth_time = -p_gen->ranf()*(p_params->GetStemCellG1Duration()
	                                             +p_params->GetSG2MDuration());
	                cell.SetNodeIndex(i);
	            	cell.SetBirthTime(birth_time);
	            	cell.SetSymmetricDivision();
	            	cells.push_back(cell);
        		}            	
        	}
        	else
        	{
        		TissueCell cell(TRANSIT, HEALTHY, 0, new CryptProjectionCellCycleModel());
            	double birth_time = -p_gen->ranf()*(p_params->GetTransitCellG1Duration()
                                             +p_params->GetSG2MDuration());
                cell.SetNodeIndex(i);
            	cell.SetBirthTime(birth_time);
            	cells.push_back(cell);            	
        	}
        }       
                                  
        // make a tissue and pass in the Wnt gradient
        Tissue<2> crypt(*p_mesh, cells);         
        crypt.SetGhostNodes(ghost_node_indices);          
        
        WntGradient::Instance()->SetType(RADIAL); 
        WntGradient::Instance()->SetTissue(crypt);

        // make a tissue simulation
        TissueSimulation<2> crypt_projection_simulator(crypt);
        
        // create cell killer and pass it in to the tissue simulation
        // find centre of crypt
        c_vector<double,2> centre(2);
        centre[0] = 0.0;
        centre[1] = 0.0;         
                
        RadialSloughingCellKiller* p_killer = new RadialSloughingCellKiller(&crypt, centre, 8.0);
        crypt_projection_simulator.AddCellKiller(p_killer);
        
        // set up other stuff
        crypt_projection_simulator.SetOutputDirectory("TestMonolayerWithKillerAndCellCycleModel");
        crypt_projection_simulator.SetEndTime(1.0);
        crypt_projection_simulator.SetMaxCells(1000);
        crypt_projection_simulator.SetMaxElements(2000);

        // TODO: should we use a cut-off?
        
        double start_time = std::clock();
        
        // run tissue simulation         
        TS_ASSERT_THROWS_NOTHING(crypt_projection_simulator.Solve());
                
        double end_time = std::clock();
        
        // print out time taken to run tissue simulation    
        double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "Time to perform simulation was = " << elapsed_time << "\n";
        
        // tidy up
        delete p_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntGradient::Destroy();
	}

};


#endif /*TESTCRYPTPROJECTIONMODEL_HPP_*/
