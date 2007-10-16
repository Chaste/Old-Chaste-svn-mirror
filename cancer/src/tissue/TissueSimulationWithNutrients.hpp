#ifndef TISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TISSUESIMULATIONWITHNUTRIENTS_HPP_

#include "TissueSimulation.cpp"
#include "EllipticFlaggedMeshAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"
#include "SimpleDataWriter.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "CellwiseData.cpp"

template<unsigned DIM>
class TissueSimulationWithNutrients : public TissueSimulation<DIM>
{
private :
    AbstractLinearEllipticPde<DIM>* mpPde;

    void PostSolve()
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = this->mrTissue.rGetMesh();

        CellwiseData<DIM>::Instance()->ReallocateMemory();
        
        std::set<unsigned> ghost_node_indices = this->mrTissue.GetGhostNodeIndices();
        r_mesh.FlagElementsNotContainingNodes(ghost_node_indices);
        r_mesh.SetupSmasrmMap();
        
        // Set up boundary conditions
        FlaggedMeshBoundaryConditionsContainer<2,1> flagged_bcc(r_mesh, 1.0);
        
        // Assembler for fine mesh flagged region. use bbc from before
        EllipticFlaggedMeshAssembler<2> elliptic_assembler(&r_mesh, mpPde, &flagged_bcc);
        
        // Solve the nutrient PDE
        Vec result_elliptic_restricted = elliptic_assembler.Solve();
        ReplicatableVector result_elliptic_repl(result_elliptic_restricted);

        std::map<unsigned, unsigned>& map = r_mesh.rGetSmasrmMap();
        std::map<unsigned, unsigned>::iterator map_iter = map.begin();

		std::vector<double> global_index;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> u;
        
        while (map_iter!=map.end())
        {
            unsigned node_index = map_iter->first;
            unsigned smasrm_index = map_iter->second;
            global_index.push_back((double) node_index);
            x.push_back(r_mesh.GetNode(node_index)->rGetLocation()[0]);
            y.push_back(r_mesh.GetNode(node_index)->rGetLocation()[1]);
            u.push_back(result_elliptic_repl[smasrm_index]);

            double o2_conc = result_elliptic_repl[smasrm_index];

            CellwiseData<DIM>::Instance()->SetValue(o2_conc, r_mesh.GetNode(node_index));
            map_iter++;
        }

//todo - using SimpleDataWriter is inefficient
        std::vector<std::vector<double> > data;
        data.push_back(global_index);
        data.push_back(x);
        data.push_back(y);
        data.push_back(u);
        
        static unsigned counter = 0;
        std::stringstream string_stream;
        string_stream << "nutrients_" << counter << ".dat";
        SimpleDataWriter writer(this->mOutputDirectory+"/nutrients/", string_stream.str(), data, false);
        counter++;
        
        VecDestroy(result_elliptic_restricted); // for the time being, while this is completely decoupled
        
        // update cells' hypoxic durations using their current oxygen concentration        
        for( typename Tissue<2>::Iterator cell_iter = this->mrTissue.Begin();
            cell_iter != this->mrTissue.End();
            ++cell_iter)
        {               
        	double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(&(*cell_iter));
        	
            // the oxygen concentration had better not be negative
        	assert(oxygen_concentration >= 0.0);
        	
//          TODO: change this line to something like
        	// if ( oxygen_concentration < CancerParameters::Instance()->GetHypoxicConcentration() )    	
        	if ( oxygen_concentration < 0.6 )
        	{
        		// add timestep to the hypoxic duration, since PostSolve() is called at the end of every timestep 
        		double curr_hyp_dur = cell_iter->GetHypoxicDuration();        		
            	cell_iter->SetHypoxicDuration(curr_hyp_dur + SimulationTime::Instance()->GetTimeStep());
        	}  
        	else // reset mHypoxicDuration
        	{
            	cell_iter->SetHypoxicDuration(0.0);
        	}          	   
        }        
        
    }


public:
    TissueSimulationWithNutrients(Tissue<DIM>& rTissue, AbstractLinearEllipticPde<DIM>* pPde, bool deleteTissue=false)
        : TissueSimulation<DIM>(rTissue, deleteTissue),
          mpPde(pPde)
    {
    }
};
#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
