#ifndef ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_
#define ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_

#include "FiniteElasticityAssembler.cpp"


template<int DIM>
struct EvaluationPointInfo
{
    Point<DIM> OldPosition;
    Point<DIM> NewPosition;
    unsigned MeshGlobalIndex;
    double SourceValue;
}; 


template<int DIM>
class AbstractGrowingTumourSourceModel
{
protected :
    std::map< unsigned, EvaluationPointInfo<DIM> > mEvaluationPoints; 

    
    
public :
    virtual void Run(double tStart, double tEnd)=0;

    
    virtual double GetSourceValue(unsigned index)
    {
        // could just do:
        //  return mEvaluationPoints[index].SourceValue;
        // but this doesn't check for existence
                
        typename std::map< unsigned, EvaluationPointInfo<DIM> >::iterator iter 
           = mEvaluationPoints.find(index);
           
        if(iter == mEvaluationPoints.end())
        {
            std::stringstream ss;
            ss << "No evaluation point corresponding to index=" << index;
            EXCEPTION(ss.str());
        }
        
        return iter->second.SourceValue;
    }
    

    void AddEvaluationPoint(unsigned index, Point<DIM> initialPosition, unsigned meshGlobalIndex)
    {
        // check an evaluation point with this index doesn't already 
        // exist
        typename std::map< unsigned, EvaluationPointInfo<DIM> >::iterator iter 
           = mEvaluationPoints.find(index);
        if(iter != mEvaluationPoints.end())
        {
            std::stringstream ss;
            ss << "An evaluation point corresponding to index=" << index << " already exists";
            EXCEPTION(ss.str());
        }             
        
        EvaluationPointInfo<DIM> evaluation_point_info;
        mEvaluationPoints[index] = evaluation_point_info;
 
        mEvaluationPoints[index].OldPosition = initialPosition;
        mEvaluationPoints[index].MeshGlobalIndex = meshGlobalIndex;
        mEvaluationPoints[index].SourceValue = 0.0;
    }

    void UpdateEvaluationPointNewPosition(FiniteElasticityAssembler<DIM>* pFiniteElasticityAssembler)
    {
        // loop over tumour points and update new position
    }
    
    unsigned GetNumEvaluationPoints()
    {
        return mEvaluationPoints.size();
    }
};

#endif /*ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_*/
