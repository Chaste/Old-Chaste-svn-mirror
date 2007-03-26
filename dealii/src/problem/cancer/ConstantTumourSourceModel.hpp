#ifndef CONSTANTTUMOURSOURCEMODEL_HPP_
#define CONSTANTTUMOURSOURCEMODEL_HPP_

#include "AbstractGrowingTumourSourceModel.hpp"

/**
 * A simple tumour source model for which s = constant at every evaluation
 * point, where the constant is taken in in the constructor
 */
template<unsigned DIM>
class ConstantTumourSourceModel : public AbstractGrowingTumourSourceModel<DIM>
{
private :
    double mValue;
public :
    ConstantTumourSourceModel(double value)
       : AbstractGrowingTumourSourceModel<DIM>(),
         mValue(value)
    {
    }
    
    void Run(double tStart, double tEnd, FiniteElasticityAssembler<DIM>* pFiniteElasticity)
    {
        typename std::map<unsigned,EvaluationPointInfo<DIM> >::iterator iter
        = this->mEvaluationPoints.begin();
        while (iter!=this->mEvaluationPoints.end())
        {
            iter->second.SourceValue = mValue;
            iter++;
        }
    }
};

#endif /*CONSTANTTUMOURSOURCEMODEL_HPP_*/
