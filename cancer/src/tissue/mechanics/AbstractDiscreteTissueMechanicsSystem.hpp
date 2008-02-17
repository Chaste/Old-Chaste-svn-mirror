#ifndef ABSTRACTDISCRETETISSUEMECHANICSSYSTEM_HPP_
#define ABSTRACTDISCRETETISSUEMECHANICSSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>

#include "AbstractTissue.cpp"

/**
 * An abstract discrete tissue mechanics system that contains basic 
 * information to all mechanics systems.
 */
template<unsigned DIM>
class AbstractDiscreteTissueMechanicsSystem
{
private :	
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //archive & mrTissue; // done in load_construct_data of subclasses
        archive & mUseCutoffPoint;
        archive & mCutoffPoint;
    }
    	
protected :
	
	AbstractTissue<DIM>* mpTissue;
    
    /** Whether to have zero force if the cells are far enough apart */
    bool mUseCutoffPoint;
    
    /** Have zero force if the cells are this distance apart (and mUseCutoffPoint==true) */
    double mCutoffPoint;

public : 
	
    AbstractDiscreteTissueMechanicsSystem()
        : mpTissue(NULL)
    {
        mUseCutoffPoint = false;
        mCutoffPoint = 1e10;
    }
    
    /**
     * Use a cutoff point, ie specify zero force if two cells are greater 
     * than the cutoff distance apart
     */
    void UseCutoffPoint(double cutoffPoint)
    {
        assert(cutoffPoint > 0.0);
        mUseCutoffPoint = true;
        mCutoffPoint = cutoffPoint;
    }
    
    /**
     *  Pure method, overloaded in the concrete classes, which returns a reference 
     *  to the node velocities. Note, the overloaded must store the velocities
     *  as a member variable (eg see the meineke implementation of this method), as 
     *  a reference to it is taken.
     * 
     *  @param drdt The velocities at each node, a std::vector of c_vector's. drdt[i](j)
     *  is the jth component of the velocity of node i.
     */
    virtual std::vector<c_vector<double, DIM> >& rCalculateVelocitiesOfEachNode()=0;
    
    virtual ~AbstractDiscreteTissueMechanicsSystem()
    {
    }

    #define COVERAGE_IGNORE
    virtual bool NeedsVoronoiTessellation()
    {
        return false;
    }
    #undef COVERAGE_IGNORE
        
    /**
     *  Get the tissue. Needed for archiving
     */
    const AbstractTissue<DIM>& rGetTissue() const
    {
        return *mpTissue;
    }
        
    AbstractTissue<DIM>* GetTissue()
    {
        return mpTissue;
    }
};

namespace boost 
{
namespace serialization 
{
template<unsigned DIM>
struct is_abstract<AbstractDiscreteTissueMechanicsSystem<DIM> > 
{
    typedef mpl::bool_<true> type;
        BOOST_STATIC_CONSTANT(bool, value = true);
};
}
}


#endif /*ABSTRACTDISCRETETISSUEMECHANICSSYSTEM_HPP_*/
