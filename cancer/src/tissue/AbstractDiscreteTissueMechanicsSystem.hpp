#ifndef ABSTRACTDISCRETETISSUEMECHANICSSYSTEM_HPP_
#define ABSTRACTDISCRETETISSUEMECHANICSSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>

#include "Tissue.cpp"

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
    }
    	
protected :
	
	Tissue<DIM>& mrTissue;

public : 
	
    AbstractDiscreteTissueMechanicsSystem(Tissue<DIM>& rTissue)
        : mrTissue(rTissue)
    {
    }
    
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
    const Tissue<DIM>& rGetTissue() const
    {
        return mrTissue;
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
