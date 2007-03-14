#ifndef GROWTHBYCONSTANTMASSODESYSTEM_HPP_
#define GROWTHBYCONSTANTMASSODESYSTEM_HPP_

#include <vector>
#include <cmath>
#include "AbstractOdeSystem.hpp"
#include "AbstractGrowingTumourSourceModel.hpp"


/**
 */
template<unsigned DIM>
class GrowthByConstantMassOdeSystem : public AbstractOdeSystem
{
private:
    double mRho;    
    unsigned mSourceModelIndex;
    
    AbstractGrowingTumourSourceModel<DIM>* mpSourceModel;
    
public:
    // Constructor
    GrowthByConstantMassOdeSystem(double rho, 
                                  unsigned sourceModelIndex, 
                                  AbstractGrowingTumourSourceModel<DIM>* pSourceModel)
       : AbstractOdeSystem(1)
    {
        mSourceModelIndex = sourceModelIndex;
        assert(rho>0);
        mRho = rho;
        assert(pSourceModel!=NULL);
        mpSourceModel = pSourceModel;
        
        mInitialConditions.push_back(1.0);
        
        SetStateVariables(mInitialConditions);  
    }
    
    ~GrowthByConstantMassOdeSystem()
    {
    }


    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
        double s = mpSourceModel->GetSourceValue(mSourceModelIndex);
        rDY[0] = (1.0/DIM)*rY[0]*mRho*s;
    }
    
}; 
#endif /*GROWTHBYCONSTANTMASSODESYSTEM_HPP_*/
