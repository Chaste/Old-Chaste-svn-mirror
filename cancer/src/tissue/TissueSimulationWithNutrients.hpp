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
        
        archive & mWriteSpheroidStatistics;
    }
    
    bool mWriteSpheroidStatistics;
    
    Vec mOxygenSolution;

    AbstractNonlinearEllipticPde<DIM>* mpPde;  
    
    /** The file that the nutrient values are written out to. */ 
    out_stream mpNutrientResultsFile; 
    
    /** The file that the spheroid statistics are written out to. */ 
    out_stream mpSpheroidStatisticsFile; 
    
    void SetWriteSpheroidStatistics()
    {
        mWriteSpheroidStatistics = true; 
        this->mrTissue.CreateVoronoiTessellation();
    }
    
    void SetupWriteNutrient() 
    { 
        OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/vis_results/",false); 
        mpNutrientResultsFile = output_file_handler.OpenOutputFile("results.viznutrient");
        *this->mpSetupFile << "Nutrient \n" ;         
    } 
    
    
    void SetupWriteSpheroidStatistics() 
    { 
        OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/vis_results/",false); 
        mpSpheroidStatisticsFile = output_file_handler.OpenOutputFile("results.vizstatistics");        
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

        for (typename Tissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
               cell_iter != this->mrTissue.End();
               ++cell_iter)
        {
            // \todo: we don't need this anymore since there are no ghost nodes,
            // but we'd need to change the visualizer before we take this out
            global_index = (double) cell_iter.GetNode()->GetIndex();
            x = cell_iter.rGetLocation()[0];
            y = cell_iter.rGetLocation()[1];
                        
            nutrient = CellwiseData<DIM>::Instance()->GetValue(&(*cell_iter));
            
            (*mpNutrientResultsFile) << global_index << " " << x << " " << y << " " << nutrient << " ";
        }
        
        (*mpNutrientResultsFile) << "\n";
    }
    
    
    void WriteSpheroidStatistics()
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        double time = p_simulation_time->GetDimensionalisedTime();
        
        (*mpSpheroidStatisticsFile) <<  time << "\t";
        
        double spheroid_radius = GetSpheroidStatistics()[0];
        double necrotic_radius = GetSpheroidStatistics()[1];

        (*mpSpheroidStatisticsFile) << time << " " << spheroid_radius << " " << necrotic_radius << "\n";                
    }
    
    
    /**     
     * Calculates the volume of the spheroid and of the necrotic core.
     * 
     * Note that this only works for DIM=2 as this is asserted in the 
     * method VoronoiTessellation<DIM>::GetFaceArea().  
     * 
     * @return a c_vector<double,2> whose first entry is the spheroid 
     *         radius and whose second entry is the necrotic radius     
     */ 
    c_vector<double,2> GetSpheroidStatistics()
    {
#define COVERAGE_IGNORE
        assert(DIM==2);
#undef COVERAGE_IGNORE
        // First get references to the Voronoi tessellation and mesh
        VoronoiTessellation<DIM>& r_tessellation = this->mrTissue.rGetVoronoiTessellation();
        ConformingTetrahedralMesh<DIM,DIM>& r_mesh = this->mrTissue.rGetMesh();
                
        double spheroid_area = 0.0;
        double necrotic_area = 0.0;
        
        // Iterate over nodes (ghost nodes are not used in this class)        
        for (unsigned i=0; i<r_mesh.GetNumAllNodes(); i++)
        { 
            // Don't use contributions from boundary nodes
            if (r_mesh.GetNode(i)->IsBoundaryNode()==false)
            {
                double cell_area = r_tessellation.GetFace(i)->GetArea();
                spheroid_area += cell_area;
                
                if (this->mrTissue.rGetCellAtNodeIndex(i).GetCellType()==NECROTIC)
                {
                    necrotic_area += cell_area;                
                }
            }            
        }     
        
        c_vector<double,2> radii;        
        radii[0] = sqrt(spheroid_area/M_PI);
        radii[1] = sqrt(necrotic_area/M_PI);        
        
        return radii;
    }
    
    
    void SetupSolve()
    {
        if ( this->mrTissue.Begin() != this->mrTissue.End() )  // there are any cells
        {
            SetupWriteNutrient();
            WriteNutrient();
            
            if (mWriteSpheroidStatistics)
            {
                SetupWriteSpheroidStatistics();
                WriteSpheroidStatistics();
            }
        }
    }

        
    void PostSolve()
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
        
        // Set up assembler
        SimpleNonlinearEllipticAssembler<DIM,DIM> assembler(&r_mesh, mpPde, &bcc);
        
        // We cannot use the exact previous solution as initial guess 
        // as the size may be different (due to cell birth/death)
        Vec initial_guess;
        
        // If we have a previous solution, then use this as the basis 
        // for the initial guess
        if (mOxygenSolution)
        {
            // Get the size of the previous solution
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
            // Free memory
            VecDestroy(mOxygenSolution);
        }
        else
        {
            initial_guess = assembler.CreateConstantInitialGuess(1.0);
        }
        
        // Solve the nutrient PDE
        mOxygenSolution = assembler.Solve(initial_guess);
        VecDestroy(initial_guess);
        ReplicatableVector result_repl(mOxygenSolution);
        
        // Update cellwise data
        for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            double oxygen_conc = result_repl[i];
            CellwiseData<DIM>::Instance()->SetValue(oxygen_conc, r_mesh.GetNode(i));
        }
        
        // Save results to file
        WriteNutrient();
        
        if (mWriteSpheroidStatistics)
        {
            WriteSpheroidStatistics();
        }
    }
    
    
    void AfterSolve()
    {
        if ( this->mrTissue.Begin() != this->mrTissue.End() )  // if there are any cells
        {
            mpNutrientResultsFile->close();
            
            if (mWriteSpheroidStatistics)
            {
                mpSpheroidStatisticsFile->close();
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
    TissueSimulationWithNutrients(Tissue<DIM>& rTissue,
                                  AbstractDiscreteTissueMechanicsSystem<DIM>* pMechanicsSystem=NULL,
                                  AbstractNonlinearEllipticPde<DIM>* pPde=NULL,
                                  bool deleteTissue=false,
                                  bool initialiseCells=true) 
        : TissueSimulation<DIM>(rTissue, pMechanicsSystem, deleteTissue, initialiseCells), 
          mWriteSpheroidStatistics(false),
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
    void SetPde(AbstractNonlinearEllipticPde<DIM>* pPde)
    {
        mpPde = pPde;
    }
    
    
    /**
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
    const Tissue<DIM> * p_tissue = &(t->rGetTissue());
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
    Tissue<DIM>* p_tissue;
    ar >> p_tissue;
    
    AbstractDiscreteTissueMechanicsSystem<DIM>* p_spring_system;
    ar >> p_spring_system;
    
    // invoke inplace constructor to initialize instance
    ::new(t)TissueSimulationWithNutrients<DIM>(*p_tissue, p_spring_system, NULL, true, false);
}
}
} // namespace ...



#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
