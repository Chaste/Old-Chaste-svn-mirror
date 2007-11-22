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
    // allow tests to access private members, in order to test computation of private functions 
    friend class TestTissueSimulationWithNutrients;
    
private :
    
    Vec mOxygenSolution;

    AbstractNonlinearEllipticPde<DIM>* mpPde;  
    
    /** The file that the nutrient values are written out to. */ 
    out_stream mpNutrientResultsFile; 
    
    void SetupWriteNutrient() 
    { 
        OutputFileHandler output_file_handler2(this->mSimulationOutputDirectory+"/vis_results/",false); 
        mpNutrientResultsFile = output_file_handler2.OpenOutputFile("results.viznutrient");
        *this->mpSetupFile << "Nutrient \n" ;         
    } 
    
    void WriteNutrient()
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        double time = p_simulation_time->GetDimensionalisedTime();
        
        (*mpNutrientResultsFile) <<  time << "\t";
        
        double global_index;
        double x;
        double y;
        double nutrient;

        for (Tissue<2>::Iterator cell_iter = this->mrTissue.Begin();
               cell_iter != this->mrTissue.End();
               ++cell_iter)
        {
            // \todo: we don't need this anymore since there are no ghost nodes,
            // but we'd need to change the visualizer before we take this out
            global_index = (double) cell_iter.GetNode()->GetIndex();
            x = cell_iter.rGetLocation()[0];
            y = cell_iter.rGetLocation()[1];
                        
            nutrient = CellwiseData<2>::Instance()->GetValue(&(*cell_iter));
            
            (*mpNutrientResultsFile) << global_index << " " << x << " " << y << " " << nutrient << " ";
        }
        
        (*mpNutrientResultsFile) << "\n";
    }
    
    void SetupSolve()
    {
        if ( this->mrTissue.Begin() != this->mrTissue.End() )  // there are any cells
        {
            SetupWriteNutrient();
            WriteNutrient();
        }
    }

    void AfterSolve()
    {
        if ( this->mrTissue.Begin() != this->mrTissue.End() )  // there are any cells
        {
            mpNutrientResultsFile->close();
        }
    }
        
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
                
        // If we have a previous solution, then use this as the basis for the initial guess 
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
        
        // update cellwise data        
        for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            double oxygen_conc = result_repl[i];
            CellwiseData<DIM>::Instance()->SetValue(oxygen_conc, r_mesh.GetNode(i));
        }
        
        // save results to file
        WriteNutrient(); 
        
        // update each cell's hypoxic duration according to its current oxygen concentration        
        for( typename Tissue<2>::Iterator cell_iter = this->mrTissue.Begin();
            cell_iter != this->mrTissue.End();
            ++cell_iter)
        {               
        	double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(&(*cell_iter));
        	
            // the oxygen concentration had better not be negative
            // assert(oxygen_concentration >= -1e-8);
        	  	
        	if ( oxygen_concentration < CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration() )
        	{
        		// add timestep to the hypoxic duration, since PostSolve() is called at the end of every timestep 
        		double curr_hyp_dur = cell_iter->GetHypoxicDuration();        		
            	cell_iter->SetHypoxicDuration(curr_hyp_dur + SimulationTime::Instance()->GetTimeStep());
        	}  
        	else 
        	{
        		// reset the cell's hypoxic duration
            	cell_iter->SetHypoxicDuration(0.0);
        	}          	   
        }  
        
        VecDestroy(initial_guess);          
    }


public:
    TissueSimulationWithNutrients(Tissue<DIM>& rTissue, AbstractDiscreteTissueMechanicsSystem<DIM>* pMechanicsSystem=NULL, AbstractNonlinearEllipticPde<DIM>* pPde=NULL, bool deleteTissue=false) 
        : TissueSimulation<DIM>(rTissue, pMechanicsSystem, deleteTissue), 
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
        this->mpSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");
        
        *this->mpSetupFile << "FinalMeshSize\t" << std::max((this)->mrTissue.rGetMesh().GetWidth(0u),(this)->mrTissue.rGetMesh().GetWidth(1u));
        this->mpSetupFile->close();
    }
    
    
};

#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
