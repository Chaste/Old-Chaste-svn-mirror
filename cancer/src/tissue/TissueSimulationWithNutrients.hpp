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

    /** 
     *  The current nutrient concentration, for use as an initial guess 
     *  when solving the nutrient PDE
     */  
    Vec mNutrientSolution;

    /** Pointer to the PDE satisfied by the nutrient. */ 
    AbstractLinearEllipticPde<DIM>* mpPde;  
    
    /** The file that the nutrient values are written out to. */ 
    out_stream mpNutrientResultsFile; 
    
    /** The file that the average radial nutrient distribution is written out to. */ 
    out_stream mpAverageRadialNutrientResultsFile;
    
    /** Whether to write to file the average radial nutrient distribution. */
    bool mWriteAverageRadialNutrientResults; 
    
    /** Whether to write the average radial nutrient distribution DAILY. */
    bool mWriteDailyAverageRadialNutrientResults;
    
    /** The number of radial 'bins' used to calculate the average radial nutrient distribution. */
    unsigned mNumRadialIntervals;
            
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
            if (mWriteAverageRadialNutrientResults)
            {
                mpAverageRadialNutrientResultsFile = output_file_handler.OpenOutputFile("radial_dist.dat");                     
            }
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
    
    void WriteAverageRadialNutrientDistribution(double time, unsigned num_intervals)
    {
        (*mpAverageRadialNutrientResultsFile) << time << " "; 
        
        // Get reference to the mesh and its size
        ConformingTetrahedralMesh<DIM,DIM>& r_mesh = this->mrTissue.rGetMesh();
        unsigned num_nodes = r_mesh.GetNumNodes();
        
        // Calculate the centre of the tissue
        c_vector<double,DIM> centre = zero_vector<double>(DIM);
        for (unsigned i=0; i< num_nodes; i++)
        {
            centre += r_mesh.GetNode(i)->rGetLocation();    
        }
        centre /= (double) num_nodes;
        
        // Calculate the distance between each node and the centre 
        // of the tissue, as well as the maximum of these
        std::map<double, TissueCell*> distance_cell_map;
        
        double max_distance_from_centre = 0.0;
                
        for (unsigned i=0; i<this->mrTissue.GetNumRealCells(); i++)
        {
            double distance = norm_2(r_mesh.GetNode(i)->rGetLocation()-centre);            
            distance_cell_map[distance] = &(this->mrTissue.rGetCellAtNodeIndex(i));
            
            if (distance > max_distance_from_centre)
            {
                max_distance_from_centre = distance;
            }
        }

        // Create vector of radius intervals                                  
        std::vector<double> radius_intervals;
        for (unsigned i=0; i<num_intervals; i++)
        {
            double upper_radius = max_distance_from_centre*((double) i+1)/((double) num_intervals);
            radius_intervals.push_back(upper_radius);       
        }
        
        // Calculate nutrient concentration in each radial interval
        double lower_radius = 0.0;
        for (unsigned i=0; i<num_intervals; i++)
        {   
            unsigned counter = 0;
            double average_conc = 0.0;
            
            for (std::map<double, TissueCell*>::iterator iter=distance_cell_map.begin();
                 iter != distance_cell_map.end();
                 ++iter)
            {
                if ((*iter).first > lower_radius && (*iter).first <= radius_intervals[i])
                {
                    average_conc += CellwiseData<DIM>::Instance()->GetValue((*iter).second); 
                    counter++;
                }
            }
            if (counter > 0)
            {
                average_conc /= (double) counter;
            }
            
            // Write results to file
            (*mpAverageRadialNutrientResultsFile) << radius_intervals[i] << " " << average_conc << " ";
            lower_radius = radius_intervals[i];
        }
        (*mpAverageRadialNutrientResultsFile) << "\n";
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
        
        // Set up assembler - note this is a purpose-made elliptic assembler
        // that interpolates the source terms from node onto gauss points,
        // as for a nutrients simulation the source will only be known at the
        // cells (nodes), not the gauss points
        TissueSimulationWithNutrientsAssembler<DIM> assembler(&r_mesh,mpPde,&bcc);
        
        PetscInt size_of_soln_previous_step = 0;
        
        if(mNutrientSolution)
        {
            VecGetSize(mNutrientSolution, &size_of_soln_previous_step);
        }
        
        if (size_of_soln_previous_step == (int)r_mesh.GetNumNodes())
        {
            // We make an initial guess which gets copied by the Solve method of
            // SimpleLinearSolver, so we need to delete it too.
            Vec initial_guess;
            VecDuplicate(mNutrientSolution, &initial_guess);
            VecCopy(mNutrientSolution, initial_guess);
            
            // Use current solution as the initial guess
            VecDestroy(mNutrientSolution);    // Solve method makes its own mNutrientSolution
            mNutrientSolution = assembler.Solve(initial_guess);
            VecDestroy(initial_guess);
        }
        else
        {
            if (mNutrientSolution)
            {
                VecDestroy(mNutrientSolution);
            }
            mNutrientSolution = assembler.Solve();
        }            

        ReplicatableVector result_repl(mNutrientSolution);

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
//        if (mNutrientSolution)
//        {
//            // Get the size of the previous solution
//            PetscInt isize;
//            VecGetSize(mNutrientSolution, &isize);
//            unsigned size_of_previous_solution = isize;
//            
//            if (size_of_previous_solution != r_mesh.GetNumNodes() )
//            {
//                initial_guess = assembler.CreateConstantInitialGuess(1.0);
//            }
//            else
//            {
//                VecDuplicate(mNutrientSolution, &initial_guess);
//                VecCopy(mNutrientSolution, initial_guess);
//            }
//            // Free memory
//            VecDestroy(mNutrientSolution);
//        }
//        else
//        {
//            initial_guess = assembler.CreateConstantInitialGuess(1.0);
//        }
//        
//        // Solve the nutrient PDE
//        mNutrientSolution = assembler.Solve(initial_guess);
//        VecDestroy(initial_guess);
//        ReplicatableVector result_repl(mNutrientSolution);

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
        
        double time_next_step = p_time->GetDimensionalisedTime() + p_time->GetTimeStep();
        
        if ((p_time->GetTimeStepsElapsed()+1)%this->mSamplingTimestepMultiple==0)
        {
            WriteNutrient(time_next_step);
        }
        
#define COVERAGE_IGNORE
        // Note: The number of timesteps per day is equal to 2880=24*120
        if ( mWriteDailyAverageRadialNutrientResults &&
             (p_time->GetTimeStepsElapsed()+1)%2880==0 )
        {
            WriteAverageRadialNutrientDistribution(time_next_step, mNumRadialIntervals);
        }
#undef COVERAGE_IGNORE

    }    
    
    void AfterSolve()
    {
        if (this->mrTissue.Begin() != this->mrTissue.End() // if there are any cells
	    && PetscTools::AmMaster())
        {
            mpNutrientResultsFile->close();
            
            if (mWriteAverageRadialNutrientResults)
            {
                mpAverageRadialNutrientResultsFile->close();
            }
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
          mNutrientSolution(NULL),
          mpPde(pPde),
          mWriteAverageRadialNutrientResults(false),
          mWriteDailyAverageRadialNutrientResults(false)
    {
    }    
    
    ~TissueSimulationWithNutrients()
    {
        if (mNutrientSolution)
        {
            VecDestroy(mNutrientSolution);
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
    
    void SetWriteAverageRadialNutrientResults(unsigned numRadialIntervals=10, bool writeDailyResults=false)
    {
        mWriteAverageRadialNutrientResults = true;
        mNumRadialIntervals = numRadialIntervals;
        mWriteDailyAverageRadialNutrientResults = writeDailyResults;
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
