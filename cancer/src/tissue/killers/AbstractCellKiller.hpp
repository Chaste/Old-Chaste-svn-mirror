#ifndef ABSTRACTCELLKILLER_HPP_
#define ABSTRACTCELLKILLER_HPP_

#include "AbstractTissue.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>

template <unsigned SPACE_DIM>
class AbstractCellKiller
{
public:
    
    AbstractCellKiller(AbstractTissue<SPACE_DIM>* pTissue);
        
    virtual ~AbstractCellKiller() {};

    /**
     *  Pure method which should call StartApoptosis() on any cell
     *  which should be about to undergo programmed death, or Kill()
     *  on any cell which should die immediately
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath()=0;
        
    const AbstractTissue<SPACE_DIM>* GetTissue() const;
    
protected:
    AbstractTissue<SPACE_DIM>* mpTissue;
    
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //archive & mpTissue; // done in load_construct_data of subclasses
    }
    
};

template <unsigned SPACE_DIM>
AbstractCellKiller<SPACE_DIM>::AbstractCellKiller(AbstractTissue<SPACE_DIM>* pTissue)
        : mpTissue(pTissue) 
{
}

template <unsigned SPACE_DIM>
const AbstractTissue<SPACE_DIM>* AbstractCellKiller<SPACE_DIM>::GetTissue() const
{
    return mpTissue;
}
    
namespace boost {
namespace serialization {
template<unsigned DIM>
struct is_abstract<AbstractCellKiller<DIM> > {
    typedef mpl::bool_<true> type;
        BOOST_STATIC_CONSTANT(bool, value = true);
};
}}


#endif /*ABSTRACTCELLKILLER_HPP_*/
