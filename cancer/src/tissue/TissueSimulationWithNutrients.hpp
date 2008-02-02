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
     * The current nutrient concentration, for use as an initial guess 
     * when solving the nutrient PDE.
     */  
    Vec mNutrientSolution;

    /** 
     * Pointer to the PDE satisfied by the nutrient. 
     */ 
    AbstractLinearEllipticPde<DIM>* mpPde;  

    /** 
     * The file that the nutrient values are written out to. 
     */ 
    out_stream mpNutrientResultsFile; 

    /**
     * The file that the average radial nutrient distribution is written out to. 
     */ 
    out_stream mpAverageRadialNutrientResultsFile;

    /** 
     * Whether to write to file the average radial nutrient distribution. 
     */
    bool mWriteAverageRadialNutrientResults; 

    /** 
     * Whether to write the average radial nutrient distribution DAILY. 
     */
    bool mWriteDailyAverageRadialNutrientResults;

    /** 
     * The number of radial 'bins' used to calculate the average 
     * radial nutrient distribution. 
     */
    unsigned mNumRadialIntervals;

    /**
     * Overridden SetupSolve() method. 
     */ 
    void SetupSolve();
    
    /**
     * Set up the nutrient writer.
     */ 
    void SetupWriteNutrient();
    
    /**
     * Write the nutrient distribution to file at a specified time.
     * 
     * @param time The time at which to record the nutrient distribution
     */
    void WriteNutrient(double time);

    /**
     * Write the average radial nutrient distribution to file at a specified time.
     * 
     * @param time The time at which to record the average radial nutrient distribution
     * @param numIntervals  The number of radial intervals in which the average nutrient concentration is calculated 
     */
    void WriteAverageRadialNutrientDistribution(double time, unsigned numIntervals);
    
    /**
     * Solve the nutrient PDE. 
     */
    void SolveNutrientPde();

    /**
     * Overridden PostSolve() method. 
     */
    void PostSolve();
    
    /**
     * Overridden AfterSolve() method. 
     */
    void AfterSolve();

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
    
    /**
     * Destructor
     * 
     * Free any memory allocated by the constructor.
     * This frees the current nutrient distribution, if it exists.
     */
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
    
    /**
     * Write the final (and optionally also the daily) average 
     * radial nutrient distribution to file.
     *
     *  @param numRadialIntervals The number of radial intervals in which the average nutrient concentration is calculated
     *  @param writeDailyResults Whether to record the average radial nutrient distribution at the end of each day of the simulation 
     */ 
    
    void SetWriteAverageRadialNutrientResults(unsigned numRadialIntervals=10, 
                                              bool writeDailyResults=false);
    
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
    void Save();
    
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
    static TissueSimulationWithNutrients<DIM>* Load(const std::string& rArchiveDirectory, 
                                                    const double& rTimeStamp);
        
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
    // Save data required to construct instance
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
    // Retrieve data from archive required to construct new instance
    MeshBasedTissue<DIM>* p_tissue;
    ar >> p_tissue;
    
    AbstractDiscreteTissueMechanicsSystem<DIM>* p_spring_system;
    ar >> p_spring_system;
    
    // Invoke inplace constructor to initialize instance
    ::new(t)TissueSimulationWithNutrients<DIM>(*p_tissue, p_spring_system, NULL, true, false);
}
}
} // namespace ...


#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
