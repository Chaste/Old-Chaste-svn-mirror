#ifndef TISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TISSUESIMULATIONWITHNUTRIENTS_HPP_

#include "TissueSimulation.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDataWriter.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "SimpleNonlinearEllipticAssembler.hpp"
#include "CellwiseData.cpp"

template<unsigned DIM>
class TissueSimulationWithNutrients : public TissueSimulation<DIM>
{
private :
    Vec mOxygenSolution;

    AbstractNonlinearEllipticPde<DIM>* mpPde;

    void PostSolve()
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = this->mrTissue.rGetMesh();

        CellwiseData<DIM>::Instance()->ReallocateMemory();
        
        std::set<unsigned> ghost_node_indices = this->mrTissue.GetGhostNodeIndices();
        
        // we shouldn't have any ghost nodes
        assert(ghost_node_indices.size()==0);
                
        // set up boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition;
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = r_mesh.GetBoundaryNodeIteratorBegin();
        while (node_iter != r_mesh.GetBoundaryNodeIteratorEnd())
        {
            p_boundary_condition = new ConstBoundaryCondition<2>(1.0);
            bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
            node_iter++;            
        }
                
        // set up assembler
        SimpleNonlinearEllipticAssembler<2,2> assembler(&r_mesh, mpPde, &bcc);
                
        // we cannot use the exact previous solution as initial guess as the size may be different
        // (due to cell birth/death) - could we just add in the necessary number of 1.0's say? 
        Vec initial_guess;
        
//        if(mOxygenSolution)
//        {
//            VecDuplicate(mOxygenSolution, &initial_guess);
//            VecCopy(mOxygenSolution, initial_guess);
//        }
//        else
//        {
//            initial_guess = assembler.CreateConstantInitialGuess(1.0);
//        }
        
        initial_guess = assembler.CreateConstantInitialGuess(1.0);
        
        // solve the nutrient PDE
        mOxygenSolution = assembler.Solve(initial_guess);        
        ReplicatableVector result_repl(mOxygenSolution);
                
        // save results to file
		std::vector<double> global_indices;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> u;
        
        for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            // we don't need to save the global node indices any more, 
            // since there are no ghost nodes, but we'd need to change 
            // the visualizer before we can take this out
            global_indices.push_back((double) i);
            x.push_back(r_mesh.GetNode(i)->rGetLocation()[0]);
            y.push_back(r_mesh.GetNode(i)->rGetLocation()[1]);
            u.push_back(result_repl[i]);

            double oxygen_conc = result_repl[i];
            CellwiseData<DIM>::Instance()->SetValue(oxygen_conc, r_mesh.GetNode(i));
        }

        // TODO: using SimpleDataWriter is inefficient
        std::vector<std::vector<double> > data;
        data.push_back(global_indices);
        data.push_back(x);
        data.push_back(y);
        data.push_back(u);
        
        static unsigned counter = 0;
        std::stringstream string_stream;
        string_stream << "nutrients_" << counter << ".dat";
        SimpleDataWriter writer(this->mOutputDirectory+"/nutrients/", string_stream.str(), data, false);
        counter++;
        
        VecDestroy(initial_guess);
        
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
    TissueSimulationWithNutrients(Tissue<DIM>& rTissue, AbstractNonlinearEllipticPde<DIM>* pPde, bool deleteTissue=false)
        : TissueSimulation<DIM>(rTissue, deleteTissue),
          mOxygenSolution(NULL),
          mpPde(pPde)
    {
    }
    
    ~TissueSimulationWithNutrients()
    {
        if(mOxygenSolution)
        {
            VecDestroy(mOxygenSolution);
        }
    }
};
#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
