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

/** A pde which calculates the source term by adding the number of cells
 *  in the element containing that point and scaling by the element area
 */
template<unsigned DIM>
class AveragedSinksPde : public AbstractLinearEllipticPde<DIM>
{
private:
    MeshBasedTissue<DIM>& mrTissue;
    double mCoefficient;
    std::vector<double> mCellDensityOnCoarseElements;

public:
    AveragedSinksPde(MeshBasedTissue<DIM>& rTissue, double coefficient)
        : mrTissue(rTissue),
          mCoefficient(coefficient)
    {
    }

    void SetupSourceTerms(ConformingTetrahedralMesh<DIM,DIM>& rCoarseMesh) // must be called before solve
    {
        // allocate memory
        mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumNodes());
        for(unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
        {
            mCellDensityOnCoarseElements[elem_index]=0.0;
        } 
        //loop over cells, find which coarse element it is in, and add 1 to the mSourceTermOnCoarseElements[elem_index];
        for(typename MeshBasedTissue<DIM>::Iterator cell_iter = mrTissue.Begin();
            cell_iter != mrTissue.End();
            ++cell_iter)
        {
            const ChastePoint<DIM>& r_position_of_cell = cell_iter.rGetLocation();
            unsigned elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
            mCellDensityOnCoarseElements[elem_index] += 1.0;
        }    
        
        // then divide each entry of mSourceTermOnCoarseElements by the element's area
        for(unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
        {
            mCellDensityOnCoarseElements[elem_index]/= rCoarseMesh.GetElement(elem_index)->GetVolume();
        }
    
    }

    double ComputeConstantInUSourceTerm(const ChastePoint<DIM>& x)
    {
        return 0.0;
    }
    
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& x) //, const Element<DIM>& rElement) // now takes in element
    {
        return mCoefficient*mCellDensityOnCoarseElements[0] ;//rElement.GetIndex()];
    }
   
    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& )
    {
        return identity_matrix<double>(DIM);
    }   
};

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
        archive & mWriteAverageRadialNutrientResults;
        archive & mWriteDailyAverageRadialNutrientResults;
        archive & mNumRadialIntervals;
    }

    /** 
     * Current nutrient concentration, for use as an initial guess 
     * when solving the nutrient PDE.
     */  
    Vec mNutrientSolution;

    /** 
     * Pointer to the PDE satisfied by the nutrient. 
     */ 
    AbstractLinearEllipticPde<DIM>* mpPde;
    
    /** 
     * Pointer to the averaged sink PDE satisfied by the nutrient. 
     */ 
    AveragedSinksPde<DIM>* mpAveragedSinksPde;  

    /** 
     * File that the nutrient values are written out to. 
     */ 
    out_stream mpNutrientResultsFile; 

    /**
     * File that the average radial nutrient distribution is written out to. 
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
     *  Number of radial 'bins' used to calculate the average 
     * radial nutrient distribution. 
     */
    unsigned mNumRadialIntervals;
  
    /**
     * Coarse nutrient mesh on which to solve the nutrient PDE.
     */
    ConformingTetrahedralMesh<DIM,DIM>* mpCoarseNutrientMesh;
    
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
    
    void SolveNutrientPdeUsingCoarseMesh();

    /**
     * Overridden PostSolve() method. 
     */
    void PostSolve();
    
    /**
     * Overridden AfterSolve() method. 
     */
    void AfterSolve();
    
    /**
     * Create a coarse mesh on which to solve the nutrient PDE.
     */ 
    void CreateCoarseNutrientMesh(double coarseGrainScaleFactor);

public:

    /** 
     * Constructor
     * 
     * @param rTissue A tissue facade class (contains a mesh and cells)
     * @
     * @param deleteTissue whether to delete the tissue on destruction to free up memory
     * @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     * 
     */
     TissueSimulationWithNutrients(MeshBasedTissue<DIM>& rTissue,
                                   AbstractDiscreteTissueMechanicsSystem<DIM>* pMechanicsSystem=NULL,
                                   AbstractLinearEllipticPde<DIM>* pPde=NULL,
                                   AveragedSinksPde<DIM>* pAveragedSinksPde = NULL,
                                   bool deleteTissue=false,
                                   bool initialiseCells=true)
        : TissueSimulation<DIM>(rTissue, pMechanicsSystem, deleteTissue, initialiseCells),
          mNutrientSolution(NULL),
          mpPde(pPde),
          mpAveragedSinksPde(pAveragedSinksPde),
          mWriteAverageRadialNutrientResults(false),
          mWriteDailyAverageRadialNutrientResults(false),
          mpCoarseNutrientMesh(NULL)
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
        if (mpCoarseNutrientMesh)
        {
            delete mpCoarseNutrientMesh;
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
     * A small hack until we fully archive this class - 
     * needed to set the PDE after loading a simulation 
     * from an archive.
     */
    void SetAveragedSinksPde(AveragedSinksPde<DIM>* pAveragedSinksPde)
    {
        mpAveragedSinksPde = pAveragedSinksPde;
    }
    
    /**
     * Write the final (and optionally also the daily) average 
     * radial nutrient distribution to file.
     *
     * @param numRadialIntervals The number of radial intervals in which the average nutrient concentration is calculated
     * @param writeDailyResults Whether to record the average radial nutrient distribution at the end of each day of the simulation 
     */ 
    
    void SetWriteAverageRadialNutrientResults(unsigned numRadialIntervals=10, 
                                              bool writeDailyResults=false);

    /**
     * Solve the nutrient PDE on a coarse mesh.
     * 
     * @param coarseGrainScaleFactor The ratio of the width of the coarse nutrient mesh to the initial width of the tissue
     */ 
    void UseCoarseNutrientMesh(double coarseGrainScaleFactor=10.0);
    
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
    ::new(t)TissueSimulationWithNutrients<DIM>(*p_tissue, p_spring_system, NULL,NULL, true, false);
}
}
} // namespace ...


#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
