#ifndef UNITYTUMOURSOURCEMODEL_HPP_
#define UNITYTUMOURSOURCEMODEL_HPP_

#include "AbstractGrowingTumourSourceModel.hpp"

/* A simple tumour source model (for testing) that just returns
 * s = n, at an evaluation point, where n is the index of the evaluation
 * point 
 */
template<int DIM>
class SimpleTumourSourceModel : public AbstractGrowingTumourSourceModel<DIM>
{
public :
    void Run(double tStart, double tEnd)
    {
        typename std::map<unsigned,EvaluationPointInfo<DIM> >::iterator iter 
           = this->mEvaluationPoints.begin();
        while(iter!=this->mEvaluationPoints.end())
        {
            unsigned index = iter->first;
            iter->second.SourceValue = index;
            iter++;
        }
    }
};

#endif /*UNITYTUMOURSOURCEMODEL_HPP_*/
