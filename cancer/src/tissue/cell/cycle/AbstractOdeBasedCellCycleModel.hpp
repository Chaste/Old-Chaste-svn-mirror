#ifndef _ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_
#define _ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractCellCycleModel.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

class AbstractOdeBasedCellCycleModel : public AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mLastTime;
        archive & mDivideTime;
        archive & mReadyToDivide;
    }
    
protected:
    /** The last time the cell cycle ODEs were evaluated.*/
    double mLastTime;
    /** The time at which the cell should divide - Set this to DBL_MAX in constructor.*/
    double mDivideTime;
    /** Whether the cell is ready to divide or not */
    bool mReadyToDivide;
    
public:
    /**
     * This overrides the AbstractCellCycleModel::SetBirthTime(double birthTime)
     * because an ODE based cell cycle model has more to reset...
     * 
     * @param birthTime the simulation time when the cell was born
     */
    void SetBirthTime(double birthTime)
    {
        mLastTime = birthTime;
        mBirthTime = birthTime;
        mDivideTime = birthTime;
    }

};

BOOST_IS_ABSTRACT(AbstractOdeBasedCellCycleModel)

#endif //_ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_
