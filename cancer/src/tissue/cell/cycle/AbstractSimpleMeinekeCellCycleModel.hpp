#ifndef _ABSTRACTSIMPLEMEINEKECELLCYCLEMODEL_HPP_
#define _ABSTRACTSIMPLEMEINEKECELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractSimpleCellCycleModel.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * This class contains all the things common to simple meineke cell cycle models
 * 
 * i.e. models where the length of cell cycle phases are determined when 
 * the cell cycle model is created, 
 * rather than evaluated 'on the fly' by ODEs and suchlike.
 * And also models which consider the generation of cells 
 * 
 * N.B. Whether or not the cell should actually divide may depend on 
 * Wnt / Oxygen etc. in subclasses...
 */
class AbstractSimpleMeinekeCellCycleModel : public AbstractSimpleCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
    }

protected:
    /** 
     * Protected constructor for creating an identical daughter cell 
     * (with the same G_ONE duration...)
     * */
    AbstractSimpleMeinekeCellCycleModel(double g1Duration, unsigned generation)
        : AbstractSimpleCellCycleModel(g1Duration, generation)
    {}

public:
    /**
     * Default constructor - creates an AbstractSimpleCellCycleModel
     */
    AbstractSimpleMeinekeCellCycleModel()
    {}
        
    /**
     * Default destructor
     */
    virtual ~AbstractSimpleMeinekeCellCycleModel()
    {}
    
    /**
     * Returns the cell types of the next generation of cells in a vector
     * [0] is the new mother cell type, [1] is the new daughter cell type
     * overwritten as new daughter cell type can depend on mother cell type.
     */
    std::vector<CellType> GetNewCellTypes();
        
    void ResetModel();
    
    void InitialiseDaughterCell();
    
};

BOOST_IS_ABSTRACT(AbstractSimpleMeinekeCellCycleModel)

#endif //_ABSTRACTSIMPLEMEINEKECELLCYCLEMODEL_HPP_
