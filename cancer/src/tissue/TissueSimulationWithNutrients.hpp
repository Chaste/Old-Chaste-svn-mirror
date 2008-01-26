#ifndef TISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TISSUESIMULATIONWITHNUTRIENTS_HPP_

#include "TissueSimulation.cpp"
#include "SimpleDataWriter.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "TissueSimulationWithNutrientsAssembler.hpp"
#include "CellwiseData.cpp"
#include "PetscTools.hpp"

template<unsigned DIM>
class TissueSimulationWithNutrients : public TissueSimulation<DIM>
{
    // Allow tests to access private members, in order to test computation of private functions 
    friend class TestTissueSimulationWithNutrients;
    
private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {   
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>      
        archive & boost::serialization::base_object<TissueSimulation<DIM> >(*this);
    }
        
    Vec mOxygenSolution;

    AbstractLinearEllipticPde<DIM>* mpPde;  
    
    /** The file that the nutrient values are written out to. */ 
    out_stream mpNutrientResultsFile; 
    
    void SetupSolve()
    {
        if (this->mrTissue.Begin() != this->mrTissue.End())
        {
            SetupWriteNutrient();           
            double current_time = SimulationTime::Instance()->GetDimensionalisedTime();            
            WriteNutrient(current_time);
        }
    }

    void SetupWriteNutrient()
    { 
        OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/vis_results/",false);
    	if (output_file_handler.IsMaster())
    	{
    	    mpNutrientResultsFile = output_file_handler.OpenOutputFile("results.viznutrient");
    	    *this->mpSetupFile << "Nutrient \n" ;
    	}
    } 
    
    void WriteNutrient(double time)
    {
    	if (PetscTools::AmMaster())
    	{
            // Since there are no ghost nodes, the number of nodes must equal the number of real cells 
            assert(this->mrTissue.rGetMesh().GetNumNodes()==this->mrTissue.GetNumRealCells());
            
            (*mpNutrientResultsFile) << time << "\t";
            
    	    unsigned global_index; 
    	    double x;
    	    double y;
    	    double nutrient;
    
    	    for (typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
    		     cell_iter != this->mrTissue.End();
    		     ++cell_iter)
    	    {
        		global_index = cell_iter.GetNode()->GetIndex();
        		x = cell_iter.rGetLocation()[0];
        		y = cell_iter.rGetLocation()[1];                    
        		nutrient = CellwiseData<DIM>::Instance()->GetValue(&(*cell_iter));                    
        		
        		(*mpNutrientResultsFile) << global_index << " " << x << " " << y << " " << nutrient << " ";
    	    }            
    	    (*mpNutrientResultsFile) << "\n";
        }
    }
    
    void SolveNutrientPde()
    {
        assert(mpPde);
        
        ConformingTetrahedralMesh<DIM,DIM>& r_mesh = this->mrTissue.rGetMesh();
        CellwiseData<DIM>::Instance()->ReallocateMemory();
        std::set<unsigned> ghost_node_indices = this->mrTissue.GetGhostNodeIndices();
        
        // We shouldn't have any ghost nodes
        assert(ghost_node_indices.size()==0);
                
        // Set up boundary conditions
        BoundaryConditionsContainer<DIM,DIM,1> bcc;
        ConstBoundaryCondition<DIM>* p_boundary_condition = new ConstBoundaryCondition<DIM>(1.0);
        for (typename ConformingTetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = r_mesh.GetBoundaryNodeIteratorBegin();
             node_iter != r_mesh.GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
        }
        
        // Set up assembler - note this is a purpose made elliptic assembler
        // that interpolates the source terms from node onto gauss points,
        // as for a nutrients simulation the source will only be known at the
        // cells (nodes), not the gauss points
        TissueSimulationWithNutrientsAssembler<DIM> assembler(&r_mesh,mpPde,&bcc);
        
        PetscInt size_of_soln_previous_step = 0;
        if(mOxygenSolution)
        {
            VecGetSize(mOxygenSolution, &size_of_soln_previous_step);
        }
        
        if (size_of_soln_previous_step == (int)r_mesh.GetNumNodes())
        {
            // We make an initial guess which gets copied by the Solve method of
            // SimpleLinearSolver, so we need to delete it too.
            Vec initial_guess;
            VecDuplicate(mOxygenSolution, &initial_guess);
            VecCopy(mOxygenSolution, initial_guess);
            
            // Use current solution as the initial guess
            VecDestroy(mOxygenSolution);    // Solve method makes its own mOxygenSolution
            mOxygenSolution = assembler.Solve(initial_guess);
            VecDestroy(initial_guess);
        }
        else
        {
            if (mOxygenSolution)
            {
                VecDestroy(mOxygenSolution);
            }
            mOxygenSolution = assembler.Solve();
        }            

        ReplicatableVector result_repl(mOxygenSolution);

////  Uncomment this for non-linear pdes and find&replace Linear for NonLinear
//
//        SimpleNonlinearEllipticAssembler<DIM,DIM> assembler(&r_mesh, mpPde, &bcc);
//        
//        // We cannot use the exact previous solution as initial guess 
//        // as the size may be different (due to cell birth/death)
//        Vec initial_guess;
//        
//        // If we have a previous solution, then use this as the basis 
//        // for the initial guess
//        if (mOxygenSolution)
//        {
//            // Get the size of the previous solution
//            PetscInt isize;
//            VecGetSize(mOxygenSolution, &isize);
//            unsigned size_of_previous_solution = isize;
//            
//            if (size_of_previous_solution != r_mesh.GetNumNodes() )
//            {
//                initial_guess = assembler.CreateConstantInitialGuess(1.0);
//            }
//            else
//            {
//                VecDuplicate(mOxygenSolution, &initial_guess);
//                VecCopy(mOxygenSolution, initial_guess);
//            }
//            // Free memory
//            VecDestroy(mOxygenSolution);
//        }
//        else
//        {
//            initial_guess = assembler.CreateConstantInitialGuess(1.0);
//        }
//        
//        // Solve the nutrient PDE
//        mOxygenSolution = assembler.Solve(initial_guess);
//        VecDestroy(initial_guess);
//        ReplicatableVector result_repl(mOxygenSolution);

        // Update cellwise data
        for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            double oxygen_conc = result_repl[i];
            CellwiseData<DIM>::Instance()->SetValue(oxygen_conc, r_mesh.GetNode(i));
        }
    }
    

    void PostSolve()
    {
        SolveNutrientPde();
        // Save results to file
        SimulationTime* p_time = SimulationTime::Instance();
        if ((p_time->GetTimeStepsElapsed()+1)%this->mSamplingTimestepMultiple==0)
        {
            double time_next_step = p_time->GetDimensionalisedTime() + p_time->GetTimeStep();
            WriteNutrient(time_next_step);
        }
    }
    
    
    void AfterSolve()
    {
        if (this->mrTissue.Begin() != this->mrTissue.End() // if there are any cells
	    && PetscTools::AmMaster())
        {
            mpNutrientResultsFile->close();
        }
    }
    
    
public:

    /** 
     *  Constructor
     * 
     *  @param rTissue A tissue facade class (contains a mesh and cells)
     *  @param deleteTissue whether to delete the tissue on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
    TissueSimulationWithNutrients(MeshBasedTissue<DIM>& rTissue,
                                  AbstractDiscreteTissueMechanicsSystem<DIM>* pMechanicsSystem=NULL,
                                  AbstractLinearEllipticPde<DIM>* pPde=NULL,
                                  bool deleteTissue=false,
                                  bool initialiseCells=true) 
        : TissueSimulation<DIM>(rTissue, pMechanicsSystem, deleteTissue, initialiseCells), 
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
    
    
    /**
     * A small hack until we fully archive this class - 
     * needed to set the PDE after loading a simulation 
     * from an archive.
     */
    void SetPde(AbstractLinearEllipticPde<DIM>* pPde)
    {
        mpPde = pPde;
    }
    
    /**
     * Saves the whole tissue simulation for restarting later.
     *
     * Puts it in the folder mOutputDirectory/archive/
     * and the file "tissue_sim_at_time_<SIMULATION TIME>.arch"
     *
     * First archives simulation time then the simulation itself.
     * 
     * Note that this method has to be implemented in this class,
     * so you save the right sort of simulation to the archive.
     * Not really sure why this is needed, but...
     */
    void Save()
    {
        CommonSave(this);
    }

     /**
     * Loads a saved tissue simulation to run further.
     *
     * @param rArchiveDirectory the name of the simulation to load
     * (specified originally by simulation.SetOutputDirectory("wherever"); )
     * @param rTimeStamp the time at which to load the simulation (this must
     * be one of the times at which simulation.Save() was called) 
     * 
     * Note that this method has to be implemented in this class, since it's 
     * a static method.
     */
    static TissueSimulationWithNutrients<DIM>* Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
    {
        std::string archive_filename =
            TissueSimulation<DIM>::GetArchivePathname(rArchiveDirectory, rTimeStamp);

        // Create an input archive
        std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
        boost::archive::text_iarchive input_arch(ifs);

        TissueSimulation<DIM>::CommonLoad(input_arch);

        TissueSimulationWithNutrients<DIM>* p_sim; 
        input_arch >> p_sim;
         
        if (p_sim->rGetTissue().rGetMesh().GetNumNodes()!=p_sim->rGetTissue().rGetCells().size()) 
        { 
            #define COVERAGE_IGNORE 
            std::stringstream string_stream; 
            string_stream << "Error in Load(), number of nodes (" << p_sim->rGetTissue().rGetMesh().GetNumNodes() 
                          << ") is not equal to the number of cells (" << p_sim->rGetTissue().rGetCells().size()  
                          << ")"; 
            EXCEPTION(string_stream.str()); 
            #undef COVERAGE_IGNORE 
        } 
          
        return p_sim;         
    }       
    
};


#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TissueSimulationWithNutrients)


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulationWithNutrients.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TissueSimulationWithNutrients<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const MeshBasedTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    
    const AbstractDiscreteTissueMechanicsSystem<DIM> * p_spring_system = &(t->rGetMechanicsSystem());
    ar & p_spring_system;
}

/**
 * De-serialize constructor parameters and initialise tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TissueSimulationWithNutrients<DIM> * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    MeshBasedTissue<DIM>* p_tissue;
    ar >> p_tissue;
    
    AbstractDiscreteTissueMechanicsSystem<DIM>* p_spring_system;
    ar >> p_spring_system;
    
    // invoke inplace constructor to initialize instance
    ::new(t)TissueSimulationWithNutrients<DIM>(*p_tissue, p_spring_system, NULL, true, false);
}
}
} // namespace ...



#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
