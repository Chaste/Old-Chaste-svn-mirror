#ifndef TISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TISSUESIMULATIONWITHNUTRIENTS_HPP_

#include "TissueSimulation.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDataWriter.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "SimpleNonlinearEllipticAssembler.hpp"
#include "CellwiseData.cpp"
#include "PetscTools.hpp"

template<unsigned DIM>
class TissueSimulationWithNutrients : public TissueSimulation<DIM>
{
    // allow tests to access private members, in order to 
    // test computation of private functions 
    friend class TestTissueSimulationWithNutrients;
    
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
        
        // If we have a previous solution, then use 
        // this as the basis for the initial guess 
        if (mOxygenSolution)
        {        	
            // get the size of the previous solution
            PetscInt isize;
            VecGetSize(mOxygenSolution, &isize);
            unsigned size_of_previous_solution = isize;  
            
            if (size_of_previous_solution != r_mesh.GetNumNodes() )
            {
                initial_guess = assembler.CreateConstantInitialGuess(1.0); 
            }  
            
//            // if it's too small, then add the necessary number of extra entries
//            if (size_of_previous_solution < r_mesh.GetNumNodes() )
//            {
//                std::cout << "Need to add entries to previous solution...\n" << std::flush;
//                std::cout << "Previous solution size =" << size_of_previous_solution << "\n" << std::flush;
//                std::cout << "Current solution size =" << r_mesh.GetNumNodes() << "\n" << std::flush;
//                
//                ReplicatableVector previous_solution(mOxygenSolution);
//                std::vector<double> correctly_sized_guess;
//                
//                // use the previous solution
//                for (unsigned i=0; i<size_of_previous_solution; i++)
//                {
//                    correctly_sized_guess.push_back(previous_solution[i]);
//                }
//                // and add the the necessary number of entries
//                for (unsigned i=size_of_previous_solution; i<r_mesh.GetNumNodes(); i++)
//                {
//                    correctly_sized_guess.push_back(1.0);
//                }
//                initial_guess = PetscTools::CreateVec(correctly_sized_guess);
//            }
//            else if (size_of_previous_solution > r_mesh.GetNumNodes() )
//            {
//                std::cout << "Need to remove entries from previous solution...\n" << std::flush;
//                std::cout << "Previous solution size =" << size_of_previous_solution << "\n" << std::flush;
//                std::cout << "Current solution size =" << r_mesh.GetNumNodes() << "\n" << std::flush;
//                
//                ReplicatableVector previous_solution(mOxygenSolution);
//                std::vector<double> correctly_sized_guess;
//                
//                // use the previous solution
//                for (unsigned i=0; i<size_of_previous_solution; i++)
//                {
//                    correctly_sized_guess.push_back(previous_solution[i]);
//                }
//                // and remove the the necessary number of entries
//                for (unsigned i=size_of_previous_solution; i<r_mesh.GetNumNodes(); i++)
//                {
//                    correctly_sized_guess.pop_back();
//                }
//                initial_guess = PetscTools::CreateVec(correctly_sized_guess);
//            }
            else
            {
                VecDuplicate(mOxygenSolution, &initial_guess);
                VecCopy(mOxygenSolution, initial_guess);
            }
        }
        else
        {
            initial_guess = assembler.CreateConstantInitialGuess(1.0);        
        }
                
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
//        	assert(oxygen_concentration >= -1e-8);
        	
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
    
    /*
     * Records the final size of the mesh, for use in the visualizer
     */    
    void WriteFinalMeshSizeForVisualizer()
    {
        double time_now = SimulationTime::Instance()->GetDimensionalisedTime();
        std::ostringstream time_string;
        time_string << time_now;
            
        std::string results_directory = (this)->mOutputDirectory +"/results_from_time_" + time_string.str();
        
        OutputFileHandler output_file_handler(results_directory+"/vis_results/",false);
        out_stream p_setup_file = output_file_handler.OpenOutputFile("results.vizsetup");
        
        *p_setup_file << "FinalMeshSize\t" << std::max((this)->mrTissue.rGetMesh().GetWidth(0u),(this)->mrTissue.rGetMesh().GetWidth(1u));
        p_setup_file->close();
    }
    
    
    
};
#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
