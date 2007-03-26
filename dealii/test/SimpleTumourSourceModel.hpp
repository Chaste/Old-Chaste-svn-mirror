#ifndef SIMPLETUMOURSOURCEMODEL_HPP_
#define SIMPLETUMOURSOURCEMODEL_HPP_

#include "AbstractGrowingTumourSourceModel.hpp"

/**
 * A simple tumour source model (for testing) that just returns
 * s = n, at an evaluation point, where n is the index of the evaluation
 * point
 */
template<unsigned DIM>
class SimpleTumourSourceModel : public AbstractGrowingTumourSourceModel<DIM>
{
public :
    SimpleTumourSourceModel()
       : AbstractGrowingTumourSourceModel<DIM>()
    {
    }

    void Run(double tStart, double tEnd, FiniteElasticityAssembler<DIM>* pFiniteElasticity)
    {
        typename std::map<unsigned,EvaluationPointInfo<DIM> >::iterator iter
        = this->mEvaluationPoints.begin();
        while (iter!=this->mEvaluationPoints.end())
        {
            unsigned index = iter->first;
            iter->second.SourceValue = index;
            iter++;
        }
    }
};

#endif /*SIMPLETUMOURSOURCEMODEL_HPP_*/
