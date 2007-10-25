#ifndef CRYPTSIMULATION2D_HPP_
#define CRYPTSIMULATION2D_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "TissueSimulation.cpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

class CryptSimulation2d : public TissueSimulation<2>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2d;    

private :
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {   
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>      
        archive & boost::serialization::base_object<TissueSimulation<2> >(*this);
        archive & mUseJiggledBottomCells;   
    }
    
    /** Whether to use a flat bottom surface or to jiggle the cells on the bottom surface */
    bool mUseJiggledBottomCells;
    
    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent 
     * and daughter in opposing directions on this axis.
     * 
     * @param node_index The parent node index
     * 
     * @return daughter_coords The coordinates for the daughter cell.
     * 
     */
    c_vector<double, 2> CalculateDividingCellCentreLocations(Tissue<2>::Iterator parentCell)
    {     
        double separation = CancerParameters::Instance()->GetDivisionSeparation();
        c_vector<double, 2> parent_coords = parentCell.rGetLocation();
        c_vector<double, 2> daughter_coords;
            
        // pick a random direction and move the parent cell backwards by 0.5*sep in that
        // direction and return the position of the daughter cell (0.5*sep forwards in the
        // random vector direction
    
        // Make a random direction vector of the required length
        c_vector<double, 2> random_vector;        
        
        double random_angle = RandomNumberGenerator::Instance()->ranf();
        random_angle *= 2.0*M_PI;
        
        random_vector(0) = 0.5*separation*cos(random_angle);
        random_vector(1) = 0.5*separation*sin(random_angle);
        
        c_vector<double, 2> proposed_new_parent_coords = parent_coords-random_vector;
        c_vector<double, 2> proposed_new_daughter_coords = parent_coords+random_vector;
        
        if (   (proposed_new_parent_coords(1) >= 0.0)
            && (proposed_new_daughter_coords(1) >= 0.0))
        {
             // We are not too close to the bottom of the tissue
            // move parent
            parent_coords = proposed_new_parent_coords;
            daughter_coords = proposed_new_daughter_coords;
        }
        else
        {   
            proposed_new_daughter_coords = parent_coords+2.0*random_vector;
            while (proposed_new_daughter_coords(1) < 0.0)
            {
                random_angle = RandomNumberGenerator::Instance()->ranf();
                random_angle *= 2.0*M_PI;
                
                random_vector(0) = separation*cos(random_angle);
                random_vector(1) = separation*sin(random_angle);
                proposed_new_daughter_coords = parent_coords+random_vector;
            }
            daughter_coords = proposed_new_daughter_coords;
        }
        
        assert(daughter_coords(1)>=0.0); // to make sure dividing cells stay in the tissue
        assert(parent_coords(1)>=0.0);   // to make sure dividing cells stay in the tissue
                
        // set the parent to use this location
        ChastePoint<2> parent_coords_point(parent_coords);
        mrTissue.MoveCell(parentCell, parent_coords_point);
        return daughter_coords;           
    }
    
    
    /**
     * Moves each node to a new position for this timestep
     *
     * @param rDrDt the x and y force components on each node.
     */
    void UpdateNodePositions(const std::vector< c_vector<double, 2> >& rDrDt)
    {
        // update ghost positions first because they do not affect the real cells
        mrTissue.UpdateGhostPositions(mDt);
        // Iterate over all cells to update their positions.
        for (Tissue<2>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            TissueCell& cell = *cell_iter;
            unsigned index = cell.GetNodeIndex();
            
            ChastePoint<2> new_point(mrTissue.rGetMesh().GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
                        
            bool is_wnt_included = WntGradient::Instance()->IsGradientSetUp();
            if (!is_wnt_included) WntGradient::Destroy();
            // stem cells are fixed if no wnt, so reset the x-value to the old x-value           
            if ((cell.GetCellType()==STEM) && (!is_wnt_included))
            {
                new_point.rGetLocation()[0] = mrTissue.rGetMesh().GetNode(index)->rGetLocation()[0];
                new_point.rGetLocation()[1] = mrTissue.rGetMesh().GetNode(index)->rGetLocation()[1];
            }
            
            // for all cells - move up if below the bottom surface
            if (new_point.rGetLocation()[1] < 0.0) 
            {
                new_point.rGetLocation()[1] = 0.0; 
                if (mUseJiggledBottomCells)
                {
                   /*  
                    * Here we give the cell a push upwards so that it doesn't  
                    * get stuck on y=0 for ever (ticket:422). 
                    *  
                    * Note that all stem cells may get moved to same height and  
                    * random numbers try to ensure we aren't left with the same  
                    * problem at a different height! 
                    */ 
                    new_point.rGetLocation()[1] = 0.05*mpRandomGenerator->ranf();
                } 
            } 
            
            // move the cell
            assert(new_point[1]>=0.0);
            mrTissue.MoveCell(cell_iter, new_point); 
                    
        }
    }       

public :            

    /** 
     *  Constructor
     * 
     *  @param rTissue A tissue facade class (contains a mesh and cells)
     *  @param deleteTissue whether to delete the tissue on destruction to free up memory.
     */
    CryptSimulation2d(Tissue<2>& rTissue, bool deleteTissue=false)
        : TissueSimulation<2>(rTissue, deleteTissue),
          mUseJiggledBottomCells(false)
    {
    }
    
    
    void UseJiggledBottomCells()
    {            
        mUseJiggledBottomCells = true;                
    }
    
    /**
     * Saves the whole tissue simulation for restarting later.
     *
     * Puts it in the folder mOutputDirectory/archive/
     * and the file "tissue_sim_at_time_<SIMULATION TIME>.arch"
     *
     * First archives simulation time then the simulation itself.
     * 
     * Note that this method has to be implemented in this class (not sure why.. )
     */
    void Save()
    {
        TissueSimulation<2>::Save();
    }
    
    /**
     * Loads a saved tissue simulation to run further.
     *
     * @param rArchiveDirectory the name of the simulation to load
     * (specified originally by simulator.SetOutputDirectory("wherever"); )
     * @param rTimeStamp the time at which to load the simulation (this must
     * be one of the times at which the simulation.Save() was called) 
     * 
     * Note that this method has to be implemented in this class (not sure why.. )
     */
    static CryptSimulation2d* Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
    {
        return (CryptSimulation2d*)TissueSimulation<2>::Load(rArchiveDirectory, rTimeStamp);
    }       
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(CryptSimulation2d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation2d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation2d * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const Tissue<2> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    
    bool archive_wnt;
    archive_wnt=WntGradient::Instance()->IsGradientSetUp();
    ar & archive_wnt;
    if (archive_wnt)
    {
        WntGradient* p_wnt_gradient = WntGradient::Instance();
        ar & *p_wnt_gradient;
        ar & p_wnt_gradient;
    }    
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation2d * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    Tissue<2>* p_tissue;

    ar >> p_tissue;
    bool archive_wnt;
    ar & archive_wnt;
    if (archive_wnt)
    {
        WntGradient* p_wnt_gradient = WntGradient::Instance();
        ar & *p_wnt_gradient;
        ar & p_wnt_gradient;
    }
    // invoke inplace constructor to initialize instance
    ::new(t)CryptSimulation2d(*p_tissue, true);
}
}
} // namespace ...

#endif /*CRYPTSIMULATION2D_HPP_*/

