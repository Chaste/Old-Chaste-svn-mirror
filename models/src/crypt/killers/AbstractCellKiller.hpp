#ifndef ABSTRACTCELLKILLER_HPP_
#define ABSTRACTCELLKILLER_HPP_

#include "MeinekeCryptCell.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Crypt.cpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>

template <unsigned SPACE_DIM>
class AbstractCellKiller
{
public:
    virtual ~AbstractCellKiller()
    {}
    
    AbstractCellKiller(Crypt<SPACE_DIM>* pCrypt)
        : mpCrypt(pCrypt)
    {
    }

    /**
     *  Pure method which should call StartApoptosis() on any cell
     *  which should be about to die
     */
    virtual void TestAndLabelCellsForApoptosis()=0;
    
    /**
     *  Remove cells labelled as dead
     *  Calls RemoveDeadCells on the crypt
     * 
     *  N.B. This now calls DeleteNodePriorToReMesh() and therefore a 
     *  ReMesh MUST be done before any element information is used.
     * 
     *  @return The number of cells removed
     */
    unsigned RemoveDeadCells()
    {
        return this->mpCrypt->RemoveDeadCells();
    }
    
    const Crypt<SPACE_DIM>* GetCrypt() const
    {
        return mpCrypt;
    }
    
protected:
    Crypt<SPACE_DIM>* mpCrypt;
    
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //archive & mpCrypt; // done in load_construct_data of subclasses
    }
    
};


namespace boost {
namespace serialization {
template<unsigned DIM>
struct is_abstract<AbstractCellKiller<DIM> > {
    typedef mpl::bool_<true> type;
        BOOST_STATIC_CONSTANT(bool, value = true);
};
}}


#endif /*ABSTRACTCELLKILLER_HPP_*/
